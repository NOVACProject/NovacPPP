#include <PPPLib/Evaluation/PostEvaluationController.h>
#include <PPPLib/Evaluation/ScanEvaluation.h>
#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/Meteorology/WindField.h>
#include <PPPLib/Evaluation/EvaluationUtils.h>
#include <PPPLib/Evaluation/PostEvaluationIO.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/VectorUtils.h>

#include <chrono>
#include <sstream>

using namespace Evaluation;
using namespace novac;

bool CPostEvaluationController::EvaluateScan(
    const novac::CString& pakFileName,
    const novac::CString& fitWindowName,
    novac::CString* txtFileName,
    CPlumeInScanProperty* plumeProperties)
{
    novac::CString errorMessage, serialNumber;
    Meteorology::CWindField windField;
    CDateTime startTime;
    novac::CSpectrumIO reader;
    CSpectrum skySpectrum;

    // The CScanFileHandler is a structure for reading the 
    //  spectral information from the scan-file
    CScanFileHandler scan(m_log);

    /** ------------- The process to evaluate a scan --------------- */

    const std::string pakFileNameStr((const char*)pakFileName);

    novac::LogContext context(novac::LogContext::FileName, novac::GetFileName(pakFileNameStr));

    // ------------------ Read the scan file -----------------------
    // --- this to make sure that the spectra in the file are ok ---
    if (!scan.CheckScanFile(context, pakFileNameStr))
    {
        std::stringstream message;
        message << "Could not read received pak file: '" << pakFileNameStr << "'";
        if (scan.m_lastErrorMessage.size() > 0)
        {
            message << " (" << scan.m_lastErrorMessage << ")";
        }
        throw PPPLib::FileIoException(message.str());
    }

    // ---------- Get the information we need about the instrument ----------

    //  Find the serial number of the spectrometer
    // reader.ReadSpectrum(pakFileName, 0, spec); // TODO: check for errors!!

    // Find the information in the configuration about this instrument.
    // Notice that these throws NotFoundException if the instrument, or its configuration could not be found.
    auto instrLocation = m_setup.GetInstrumentLocation(scan.GetDeviceSerial(), scan.GetScanStartTime());
    auto fitWindow = m_setup.GetFitWindow(scan.GetDeviceSerial(), scan.m_channel, scan.GetScanStartTime(), &fitWindowName);
    auto darkSettings = m_setup.GetDarkCorrection(scan.m_device, scan.m_startTime);

    // TODO: Should the model name be required?
    const SpectrometerModel spectrometerModel = CSpectrometerDatabase::GetInstance().GetModel(instrLocation.m_spectrometerModel);
    context = context.With(novac::LogContext::DeviceModel, spectrometerModel.modelName);

    // Check if we have already evaluated this scan. Only if this is a re-run of
    // an old processing...
    if (m_userSettings.m_fIsContinuation)
    {
        if (this->m_continuation.IsPreviouslyIgnored(pakFileNameStr))
        {
            m_log.Information(context, "Scan has already been evaluated and was ignored. Will proceed to the next scan");
            return true;
        }
        else
        {
            novac::CString archivePakFileName, archiveTxtFileName;

            // loop through all possible measurement modes and see if the evaluation log file already exists
            MEASUREMENT_MODE modes[] = { MEASUREMENT_MODE::MODE_FLUX, MEASUREMENT_MODE::MODE_WINDSPEED, MEASUREMENT_MODE::MODE_STRATOSPHERE, MEASUREMENT_MODE::MODE_DIRECT_SUN,
                                        MEASUREMENT_MODE::MODE_COMPOSITION, MEASUREMENT_MODE::MODE_LUNAR, MEASUREMENT_MODE::MODE_TROPOSPHERE, MEASUREMENT_MODE::MODE_MAXDOAS };
            for (int k = 0; k < 8; ++k)
            {
                PostEvaluationIO::GetArchivingfileName(m_log, archivePakFileName, archiveTxtFileName, fitWindowName, pakFileName, m_userSettings.m_outputDirectory.std_str(), modes[k]);
                if (Filesystem::IsExistingFile(archiveTxtFileName))
                {
                    m_log.Information(context, "Scan has already been evaluated and was ignored. Will proceed to the next scan");

                    txtFileName->Format(archiveTxtFileName);

                    return true;
                }
            }
        }
    }

    // ------------- Check that the measurement is good enough to evaluate -----------
    ReasonForScanRejection rejectionReason;
    std::string errorMessageStr;
    if (!IsGoodEnoughToEvaluate(scan, fitWindow, spectrometerModel, instrLocation, m_userSettings, rejectionReason, errorMessageStr))
    {
        m_processingStats.InsertRejection(scan.m_device, rejectionReason);
        m_log.Information(context, "Scan quality not good enough to evaluate (" + errorMessageStr + "), skipping scan.");
        return false;
    }

    // 6. Evaluate the scan
    auto startEvaluation = std::chrono::steady_clock::now();
    CScanEvaluation ev{ m_userSettings, m_log };
    std::unique_ptr<CScanResult> lastResult = ev.EvaluateScan(context, scan, fitWindow, spectrometerModel, &darkSettings);
    auto stopEvaluation = std::chrono::steady_clock::now();

    // 7. Check the reasonability of the evaluation
    if (lastResult == nullptr || lastResult->GetEvaluatedNum() == 0)
    {
        m_log.Information(context, "Zero spectra evaluated in recieved pak-file. Evaluation failed.");
        return false;
    }

    const int specieIndex = lastResult->GetSpecieIndex(CMolecule(m_userSettings.m_molecule).name);

    {
        std::stringstream msg;
        std::vector<double> columns = novac::GetColumns(*lastResult, specieIndex);
        std::vector<double> columnErrors = novac::GetColumnErrors(*lastResult, specieIndex);
        auto minMax = MinMax(columns);
        msg << "Scan evaluated in " << std::chrono::duration_cast<std::chrono::milliseconds>(stopEvaluation - startEvaluation).count() << " ms.";
        msg << " " << lastResult->GetEvaluatedNum() << " spectra evaluated.";
        msg << " Min column : " << minMax.first << ", max column : " << minMax.second;
        msg << ". Median column error: " << Median(columnErrors);
        m_log.Information(context, msg.str());
    }

    // 9. Get the mode of the evaluation
    lastResult->CheckMeasurementMode();
    lastResult->GetStartTime(0, startTime);

    // 10. Append the results to the evaluation-summary log
    PostEvaluationIO::AppendToEvaluationSummaryFile(m_userSettings.m_outputDirectory.std_str(), lastResult, &scan, &instrLocation, &fitWindow, windField);
    PostEvaluationIO::AppendToPakFileSummaryFile(m_userSettings.m_outputDirectory.std_str(), lastResult, &scan, &instrLocation, &fitWindow, windField);

    // 10. Append the result to the log file of the corresponding scanningInstrument
    if (RETURN_CODE::SUCCESS != PostEvaluationIO::WriteEvaluationResult(m_log, context, m_userSettings.m_outputDirectory.std_str(), spectrometerModel, lastResult, &scan, &instrLocation, &fitWindow, windField, txtFileName))
    {
        errorMessage.Format("Failed to write evaluation log file %s. No result produced", txtFileName);
        m_log.Error(context, errorMessage.std_str());
    }

    // 11. If this was a flux-measurement then we need to see the plume for the measurement to be useful
    //  this check should only be performed on the main fit window.
    if (Equals(fitWindow.name, m_userSettings.m_fitWindowsToUse[m_userSettings.m_mainFitWindow]))
    {
        if (!IsGoodEnoughToCalculateFlux(context, lastResult))
        {
            return false;
        }
    }

    // 12. Return the properties of the scan
    if (plumeProperties != nullptr)
    {
        lastResult->GetCalculatedPlumeProperties(*plumeProperties);

        std::stringstream msg;
        msg << "Scan sees the plume at scan angle: " << plumeProperties->plumeCenter << " +-" << plumeProperties->plumeCenterError << " [deg]. Completeness: " << plumeProperties->completeness;
        m_log.Information(context, msg.str());
    }

    PostEvaluationIO::CreatePlumespectrumFile(m_log, context, m_userSettings.m_outputDirectory.std_str(), lastResult, fitWindowName, scan, spectrometerModel, plumeProperties, specieIndex);

    return true;
}

bool CPostEvaluationController::IsGoodEnoughToCalculateFlux(LogContext context, std::unique_ptr<CScanResult>& result) const
{
    if (MEASUREMENT_MODE::MODE_FLUX != result->GetMeasurementMode())
    {
        m_log.Information(context, "Scan is not a flux measurement, no flux will be calculated.");
        return false;
    }

    if (0 != result->CalculateOffset(CMolecule(m_userSettings.m_molecule)))
    {
        m_log.Information(context, "Failed to calculate the offset of the scan, no flux will be calculated.");
        return false;
    }

    std::string message;
    if (!result->CalculatePlumeCentre(CMolecule(m_userSettings.m_molecule), message))
    {
        m_log.Information(context, message + " Scan does not see the plume, no flux will be calculated.");
        return false;
    }

    return true;
}
