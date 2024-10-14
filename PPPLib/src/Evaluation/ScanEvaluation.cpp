#include <PPPLib/Evaluation/ScanEvaluation.h>
#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include <PPPLib/Logging.h>

// we want to make some statistics on the processing
#include <PPPLib/PostProcessingStatistics.h>

#include <cstdint>
#include <sstream>

using namespace Evaluation;
using namespace novac;

CScanEvaluation::CScanEvaluation(const Configuration::CUserConfiguration& userSettings, novac::ILogger& log)
    : ScanEvaluationBase(log), m_userSettings(userSettings)
{
}

CScanEvaluation::~CScanEvaluation()
{
}

std::unique_ptr<CScanResult> CScanEvaluation::EvaluateScan(
    novac::LogContext context,
    novac::CScanFileHandler& scan,
    const CFitWindow& fitWindow,
    const novac::SpectrometerModel& spectrometerModel,
    const Configuration::CDarkSettings* darkSettings)
{
    ValidateSetup(context, fitWindow); // Verify that the setup of the fit window is ok. Throws exception if it isn't

    std::unique_ptr<CEvaluationBase> eval; // the evaluator
    CFitWindow adjustedFitWindow = fitWindow; // we may need to make some small adjustments to the fit-window. This is a modified copy

    // Adjust the fit-low and fit-high parameters according to the spectra
    m_fitLow = adjustedFitWindow.fitLow;
    m_fitHigh = adjustedFitWindow.fitHigh;

    // sometimes the length of the spectra is not what we expect, 
    // we need to be able to handle this.
    adjustedFitWindow.interlaceStep = scan.GetInterlaceSteps();
    adjustedFitWindow.specLength = scan.GetSpectrumLength() * adjustedFitWindow.interlaceStep;
    adjustedFitWindow.startChannel = scan.GetStartChannel();

    // Now choose what we should do before the real evaluation. Should we;
    // 1) find the shift & squeeze from the Fraunhofer spectrum
    // 2) find the optimal shift & squeeze from the spectrum with the highest column
    // 3) do none of the above

    if (adjustedFitWindow.fraunhoferRef.m_path.size() > 4)
    {
        m_log.Debug(context, "Determining shift from FraunhoferReference");
        this->m_lastErrorMessage.clear();

        int result = adjustedFitWindow.fraunhoferRef.ReadCrossSectionDataFromFile();
        if (result != 0)
        {
            std::stringstream msg;
            msg << "Failed to read Fraunhofer reference spectrum file: " << adjustedFitWindow.fraunhoferRef.m_path;
            throw InvalidReferenceException(msg.str());
        }

        // If we have a solar-spectrum that we can use to determine the shift
        // & squeeze then fit that first so that we know the wavelength calibration
        novac::CEvaluationBase* newEval = FindOptimumShiftAndSqueezeFromFraunhoferReference(context, adjustedFitWindow, spectrometerModel, *darkSettings, m_userSettings.sky, scan);

        if (newEval != nullptr)
        {
            eval.reset(newEval);
        }
    }
    else if (fitWindow.findOptimalShift)
    {
        // Find the optimal shift & squeeze from the spectrum with the highest column
        CFitWindow window2 = adjustedFitWindow;
        for (int k = 0; k < window2.nRef; ++k)
        {
            window2.ref[k].m_shiftOption = SHIFT_TYPE::SHIFT_FIX;
            window2.ref[k].m_squeezeOption = SHIFT_TYPE::SHIFT_FIX;
            window2.ref[k].m_shiftValue = 0.0;
            window2.ref[k].m_squeezeValue = 1.0;
        }
        eval.reset(new CEvaluationBase(window2, m_log));

        // evaluate the scan one time
        std::unique_ptr<CScanResult> result = EvaluateOpenedScan(context, scan, eval, spectrometerModel, darkSettings);
        if (result == nullptr)
        {
            return 0;
        }

        if (m_indexOfMostAbsorbingSpectrum < 0)
        {
            m_log.Information(context, "Could not determine optimal shift & squeeze. No good spectra in scan");
            return 0;
        }

        // Make sure that this spectrum was ok and that the column-value is high enough
        int specieNum = 0; // TODO: Is this the correct specie to check for?
        double columnError = result->GetColumnError(m_indexOfMostAbsorbingSpectrum, specieNum); // <-- the column error that corresponds to the highest column-value
        double highestColumn = result->GetColumn(m_indexOfMostAbsorbingSpectrum, specieNum);
        if (highestColumn < 2 * columnError)
        {
            CString message;
            message.Format("Could not determine optimal shift & squeeze. Maximum column is too low (%.1lf +- %.1lf", highestColumn, columnError);
            m_log.Information(context, message.std_str());
            return 0;
        }

        novac::CEvaluationBase* newEval = FindOptimumShiftAndSqueeze(context, adjustedFitWindow, m_indexOfMostAbsorbingSpectrum, scan, darkSettings);
        if (newEval != nullptr)
        {
            eval.reset(newEval);
        }
    }

    if (eval == nullptr)
    {
        // The options above didn't apply, use the default.
        eval.reset(new CEvaluationBase(adjustedFitWindow, m_log));
    }

    // Make the real evaluation of the scan
    auto result = EvaluateOpenedScan(context, scan, eval, spectrometerModel, darkSettings);

    return result;
}

std::unique_ptr<CScanResult> CScanEvaluation::EvaluateOpenedScan(
    novac::LogContext logContext,
    novac::CScanFileHandler& scan,
    std::unique_ptr<novac::CEvaluationBase>& eval,
    const novac::SpectrometerModel& spectrometer,
    const Configuration::CDarkSettings* darkSettings)
{
    novac::CString message; // used for ShowMessage messages
    int curSpectrumIndex = 0;  // keeping track of the index of the current spectrum into the .pak-file
    double highestColumnInScan = 0.0; // the highest column-value in the evaluation

    CSpectrum dark, current;

    // ----------- Get the sky spectrum --------------
    // Get the sky and dark spectra and divide them by the number of 
    //     co-added spectra in it
    CSpectrum sky;
    if (!GetSky(scan, m_userSettings.sky, sky))
    {
        m_log.Error(logContext, "Failed to get the sky spectrum. Scan evaluation failed.");
        return nullptr;
    }
    CSpectrum skySpecBeforeDarkCorrection = sky;

    if (m_userSettings.sky.skyOption != Configuration::SKY_OPTION::USER_SUPPLIED)
    {
        // Get the dark-spectrum and remove it from the sky
        if (!GetDark(scan, sky, dark, darkSettings))
        {
            m_log.Error(logContext, "Failed to get the dark spectrum. Scan evaluation failed.");
            return nullptr;
        }
        sky.Sub(dark);
    }

    const double skyMaximumSaturationRatioInFitRegion = novac::GetMaximumSaturationRatioOfSpectrum(sky, spectrometer, m_fitLow, m_fitHigh);
    if (skyMaximumSaturationRatioInFitRegion < m_userSettings.m_minimumSaturationInFitRegion)
    {
        const double fitIntensity = sky.MaxValue(m_fitLow, m_fitHigh);
        message.Format("Sky spectrum has maximum saturation %.3lf (%.0lf counts) in fit region (at least %.3lf required). Skipping scan.", skyMaximumSaturationRatioInFitRegion, fitIntensity, m_userSettings.m_minimumSaturationInFitRegion);
        m_log.Information(logContext, message.std_str());
        return nullptr;
    }
    if (skyMaximumSaturationRatioInFitRegion > 0.95)
    {
        const double fitIntensity = sky.MaxValue(m_fitLow, m_fitHigh);
        message.Format("Sky spectrum has maximum saturation %.3lf (%.0lf counts) in fit region and judged to be saturated. Skipping scan.", skyMaximumSaturationRatioInFitRegion, fitIntensity, m_userSettings.m_minimumSaturationInFitRegion);
        m_log.Information(logContext, message.std_str());
        return nullptr;
    }

    if (sky.NumSpectra() > 0 && !m_averagedSpectra)
    {
        sky.Div(sky.NumSpectra());
        skySpecBeforeDarkCorrection.Div(skySpecBeforeDarkCorrection.NumSpectra());
    }

    // tell the evaluator which sky-spectrum to use
    eval->SetSkySpectrum(sky);

    // Adjust the fit-low and fit-high parameters according to the spectra
    m_fitLow -= sky.m_info.m_startChannel;
    m_fitHigh -= sky.m_info.m_startChannel;

    curSpectrumIndex = -1; // we're at spectrum number 0 in the .pak-file
    m_indexOfMostAbsorbingSpectrum = -1; // as far as we know, there's no absorption in any spectrum...

    // the data structure to keep track of the evaluation results
    std::unique_ptr<CScanResult> result = std::make_unique<CScanResult>();
    result->SetSkySpecInfo(skySpecBeforeDarkCorrection.m_info);
    result->SetDarkSpecInfo(dark.m_info);

    // Make sure that we'll start with the first spectrum in the scan
    scan.ResetCounter();

    // Evaluate all the spectra in the scan.
    while (1)
    {
        // remember which spectrum we're at
        const int spectrumIndex = current.ScanIndex();

        // a. Read the next spectrum from the file
        int ret = scan.GetNextSpectrum(logContext, current);

        if (ret == 0)
        {
            // if something went wrong when reading the spectrum
            if (scan.m_lastError == novac::CSpectrumIO::ERROR_SPECTRUM_NOT_FOUND || scan.m_lastError == novac::CSpectrumIO::ERROR_EOF)
            {
                // at the end of the file, quit the 'while' loop
                break;
            }
            else
            {
                novac::CString errMsg = "Faulty spectrum found in pak file.";
                switch (scan.m_lastError)
                {
                case novac::CSpectrumIO::ERROR_CHECKSUM_MISMATCH:
                    errMsg.Append(", Checksum mismatch. Spectrum ignored");
                    break;
                case novac::CSpectrumIO::ERROR_DECOMPRESS:
                    errMsg.Append(", Decompression error. Spectrum ignored");
                    break;
                default:
                    errMsg.Append(", Unknown error. Spectrum ignored");
                }
                m_log.Error(logContext, errMsg.std_str());
                // remember that this spectrum is corrupted
                result->MarkAsCorrupted(curSpectrumIndex);
                continue;
            }
        }

        ++curSpectrumIndex; // we'have just read the next spectrum in the .pak-file

        // If the read spectrum is the sky or the dark spectrum, 
        // then don't evaluate it...
        if (current.ScanIndex() == sky.ScanIndex() || current.ScanIndex() == dark.ScanIndex())
        {
            continue;
        }

        // If the spectrum is read out in an interlaced way then interpolate it back to it's original state
        if (current.m_info.m_interlaceStep > 1)
        {
            current.InterpolateSpectrum();
        }

        // b. Get the dark spectrum for this measured spectrum
        if (!GetDark(scan, current, dark, darkSettings))
        {
            m_log.Error(logContext, "Failed to get the dark spectrum for spectrum in scan. Scan evaluation failed.");
            return nullptr;
        }

        // b. Calculate the intensities, before we divide by the number of spectra
        //  and before we subtract the dark
        current.m_info.m_peakIntensity = (float)current.MaxValue(0, current.m_length - 2);
        current.m_info.m_fitIntensity = (float)current.MaxValue(m_fitLow, m_fitHigh);

        // Check if this spectrum is worth evaluating
        const double spectrumMaximumSaturationRatioInFitRegion = novac::GetMaximumSaturationRatioOfSpectrum(current, spectrometer, m_fitLow, m_fitHigh);
        if (spectrumMaximumSaturationRatioInFitRegion < m_userSettings.m_minimumSaturationInFitRegion)
        {
            message.Format("ignoring spectrum %d with maximum saturation %.3lf (%.0lf counts) in fit region (at least %.3lf required)", curSpectrumIndex, spectrumMaximumSaturationRatioInFitRegion, current.m_info.m_fitIntensity, m_userSettings.m_minimumSaturationInFitRegion);
            m_log.Information(logContext, message.std_str());
            continue;
        }

        // c. Divide the measured spectrum with the number of co-added spectra
        //     The sky and dark spectra should already be divided before this loop.
        if (current.NumSpectra() > 0 && !m_averagedSpectra)
        {
            current.Div(current.NumSpectra());
        }

        if (dark.NumSpectra() > 0 && !m_averagedSpectra)
        {
            dark.Div(dark.NumSpectra());
        }

        // Remove the dark current spectrum
        current.Sub(dark);

        // e. Evaluate the spectrum
        if (eval->Evaluate(current))
        {
            message.Format("Failed to evaluate spectrum %d out of %d in scan.", current.ScanIndex(), current.SpectraPerScan());
            if (eval->m_lastError.size() > 0)
            {
                message.AppendFormat("(%s)", eval->m_lastError.c_str());
            }

            m_log.Information(logContext, message.std_str());
            continue;
        }

        // e. Save the evaluation result
        result->AppendResult(eval->GetEvaluationResult(), current.m_info);

        // f. Check if this was an ok data point (CScanResult)
        result->CheckGoodnessOfFit(current.m_info, &spectrometer);

        // g. If it is ok, then check if the value is higher than any of the previous ones
        if (result->IsOk(result->GetEvaluatedNum() - 1) && std::abs(result->GetColumn(result->GetEvaluatedNum() - 1, 0)) > highestColumnInScan)
        {
            highestColumnInScan = std::abs(result->GetColumn(result->GetEvaluatedNum() - 1, 0));
            m_indexOfMostAbsorbingSpectrum = spectrumIndex;
        }
    } // end while(1)

    return result;
}

bool CScanEvaluation::GetDark(novac::CScanFileHandler& scan, const CSpectrum& spec, CSpectrum& dark, const Configuration::CDarkSettings* darkSettings)
{
    m_lastErrorMessage = "";
    const bool successs = ScanEvaluationBase::GetDark(scan, spec, dark, darkSettings);

    if (m_lastErrorMessage.size() > 0)
    {
        m_log.Error(m_lastErrorMessage);
    }

    return successs;
}

bool CScanEvaluation::GetSky(novac::CScanFileHandler& scan, const Configuration::CSkySettings& settings, CSpectrum& sky)
{
    m_lastErrorMessage = "";
    const bool successs = ScanEvaluationBase::GetSky(scan, settings, sky);

    if (m_lastErrorMessage.size() > 0)
    {
        m_log.Error(m_lastErrorMessage);
    }

    return successs;
}

// Sets the first reference to 'shfit free' and links the shift of all the other refernces to the first.
static void SetupFitWindowFitShiftDetermination(CFitWindow& window)
{
    window.ref[0].m_shiftOption = SHIFT_TYPE::SHIFT_FREE;
    window.ref[0].m_squeezeOption = SHIFT_TYPE::SHIFT_FIX;
    window.ref[0].m_squeezeValue = 1.0;
    for (int k = 1; k < window.nRef; ++k)
    {
        if (novac::Equals(window.ref[k].m_specieName, "FraunhoferRef"))
        {
            continue;
        }

        window.ref[k].m_shiftOption = SHIFT_TYPE::SHIFT_LINK;
        window.ref[k].m_squeezeOption = SHIFT_TYPE::SHIFT_LINK;
        window.ref[k].m_shiftValue = 0.0;
        window.ref[k].m_squeezeValue = 0.0;
    }
}

CEvaluationBase* CScanEvaluation::FindOptimumShiftAndSqueeze(novac::LogContext context, const CFitWindow& fitWindow, int indexOfMostAbsorbingSpectrum, novac::CScanFileHandler& scan, const Configuration::CDarkSettings* darkSettings)
{
    novac::CString message;
    CSpectrum spec, dark;

    // Evaluate this spectrum again with free (and linked) shift
    CFitWindow fitWindow2 = fitWindow;
    SetupFitWindowFitShiftDetermination(fitWindow2);

    // Get the sky-spectrum
    CSpectrum sky;
    if (!GetSky(scan, m_userSettings.sky, sky))
    {
        m_log.Error(context, "Failed to get the sky spectrum. Finding optimium shift and squeeze failed.");
        return nullptr;
    }

    if (m_userSettings.sky.skyOption != Configuration::SKY_OPTION::USER_SUPPLIED)
    {
        // Get the dark-spectrum and remove it from the sky
        if (!GetDark(scan, sky, dark, darkSettings))
        {
            m_log.Error(context, "Failed to get the dark spectrum. Finding optimium shift and squeeze failed.");
            return nullptr;
        }
        sky.Sub(dark);
    }

    if (sky.NumSpectra() > 0 && !m_averagedSpectra)
    {
        sky.Div(sky.NumSpectra());
    }

    // create the new evaluator
    std::unique_ptr<CEvaluationBase> intermediateEvaluator = std::make_unique<CEvaluationBase>(fitWindow2, m_log);
    intermediateEvaluator->SetSkySpectrum(sky);

    // Get the measured spectrum
    int ret = scan.GetSpectrum(context, spec, 2 + indexOfMostAbsorbingSpectrum); // The two comes from the sky and the dark spectra in the beginning
    if (ret != 1)
    {
        message.Format("Failed to read spectrum %d in file. Finding optimium shift and squeeze failed.", indexOfMostAbsorbingSpectrum);
        m_log.Error(context, message.std_str());
        return nullptr;
    }

    // Tell the user
    message.Format("Re-evaluating spectrum number %d (scan angle %.1lf, started at %2d:%2d:%2d) to determine optimum shift and squeeze", indexOfMostAbsorbingSpectrum, spec.ScanAngle(), spec.m_info.m_startTime.hour, spec.m_info.m_startTime.minute, spec.m_info.m_startTime.second);
    m_log.Information(context, message.std_str());

    if (spec.m_info.m_interlaceStep > 1)
    {
        spec.InterpolateSpectrum();
    }

    // Get the dark-spectrum and remove it
    if (!GetDark(scan, spec, dark, darkSettings))
    {
        m_log.Error(context, "Failed to get the dark spectrum for spectrum in scan. Finding optimium shift and squeeze failed.");
        return nullptr;
    }

    if (spec.NumSpectra() > 0 && !m_averagedSpectra)
    {
        spec.Div(spec.NumSpectra());
    }

    spec.Sub(dark);

    // Evaluate
    if (intermediateEvaluator->Evaluate(spec, 5000))
    {
        message.Format("Failed to evaluate spectrum. Finding optimium shift and squeeze failed.");
        if (intermediateEvaluator->m_lastError.size() > 0)
        {
            message.AppendFormat("Error message: '%s'", intermediateEvaluator->m_lastError.c_str());
        }

        m_log.Information(context, message.std_str());
        return nullptr;
    }

    // 4. See what the optimum value for the shift turned out to be.
    CEvaluationResult newResult = intermediateEvaluator->GetEvaluationResult();
    double optimumShift = newResult.m_referenceResult[0].m_shift;
    double optimumSqueeze = newResult.m_referenceResult[0].m_squeeze;

    // 5. Set the shift for all references to this value
    for (int k = 0; k < fitWindow2.nRef; ++k)
    {
        if (novac::Equals(fitWindow2.ref[k].m_specieName, "FraunhoferRef"))
        {
            continue;
        }

        fitWindow2.ref[k].m_shiftOption = SHIFT_TYPE::SHIFT_FIX;
        fitWindow2.ref[k].m_squeezeOption = SHIFT_TYPE::SHIFT_FIX;
        fitWindow2.ref[k].m_shiftValue = optimumShift;
        fitWindow2.ref[k].m_squeezeValue = optimumSqueeze;
    }

    CEvaluationBase* newEvaluator = new CEvaluationBase(fitWindow2, m_log);
    newEvaluator->SetSkySpectrum(sky);

    // 6. We're done!
    message.Format("Optimum shift set to : %.2lf. Optimum squeeze set to: %.2lf ", optimumShift, optimumSqueeze);
    m_log.Information(context, message.std_str());

    return newEvaluator;
}

void CScanEvaluation::ValidateSetup(novac::LogContext context, const novac::CFitWindow& window)
{
    if (window.fitHigh <= window.fitLow)
    {
        throw std::invalid_argument("The given fit window has an empty (fitLow, fitHigh) range");
    }
    if (window.nRef == 0)
    {
        throw std::invalid_argument("The given fit window has no references defined");
    }

    std::vector<std::string> paths;
    for (int refIdx = 0; refIdx < window.nRef; ++refIdx)
    {
        if (window.ref[refIdx].m_data == nullptr)
        {
            throw std::invalid_argument("At least one of the references of the fit window has no data (not read from disk?).");
        }

        for (const std::string& path : paths)
        {
            if (novac::Equals(path, window.ref[refIdx].m_path))
            {
                throw std::invalid_argument("The given fit window has one reference defined multiple times (" + path + ")");
            }
        }

        window.ref[refIdx].VerifyReferenceValues(window.fitLow, window.fitHigh);

        paths.push_back(window.ref[refIdx].m_path);
    }
}
