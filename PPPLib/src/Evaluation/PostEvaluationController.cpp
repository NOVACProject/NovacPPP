#include <cstring>
#include <PPPLib/Evaluation/PostEvaluationController.h>
#include <PPPLib/Evaluation/ScanEvaluation.h>
#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/Meteorology/WindField.h>
#include <PPPLib/Evaluation/EvaluationUtils.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

// This is the information we need to continue an old processing
#include <PPPLib/ContinuationOfProcessing.h>

// ... support for handling the evaluation-log files...
#include <PPPLib/File/EvaluationLogFileHandler.h>

#include <Poco/Path.h>
#include <cmath>
#include <sstream>

using namespace Evaluation;
using namespace FileHandler;
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
    CScanFileHandler scan;

    /** ------------- The process to evaluate a scan --------------- */

    // ------------------ Read the scan file -----------------------
    // --- this to make sure that the spectra in the file are ok ---
    const std::string pakFileNameStr((const char*)pakFileName);
    if (!scan.CheckScanFile(pakFileNameStr))
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

    // Check if we have already evaluated this scan. Only if this is a re-run of
    // an old processing...
    if (m_userSettings.m_fIsContinuation)
    {
        if (this->m_continuation.IsPreviouslyIgnored(pakFileNameStr))
        {
            errorMessage.Format(" Scan %s has already been evaluated and was ignored. Will proceed to the next scan", pakFileNameStr.c_str());
            m_log.Information(errorMessage.std_str());
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
                GetArchivingfileName(archivePakFileName, archiveTxtFileName, fitWindowName, pakFileName, m_userSettings.m_outputDirectory.std_str(), modes[k]);
                if (Filesystem::IsExistingFile(archiveTxtFileName))
                {
                    errorMessage.Format(" Scan %s has already been evaluated. Will proceed to the next scan", (const char*)pakFileName);
                    m_log.Information(errorMessage.std_str());

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
        m_log.Information("Scan quality not good enough to evaluate (" + errorMessageStr + "), skipping pak file : " + pakFileNameStr);
        return false;
    }

    // 6. Evaluate the scan
    CScanEvaluation ev{ m_userSettings, m_log };
    std::unique_ptr<CScanResult> lastResult = ev.EvaluateScan(scan, fitWindow, spectrometerModel, &darkSettings);

    // 7. Check the reasonability of the evaluation
    if (lastResult == nullptr || lastResult->GetEvaluatedNum() == 0)
    {
        errorMessage.Format("Zero spectra evaluated in recieved pak-file %s. Evaluation failed.", (const char*)pakFileName);
        ShowMessage(errorMessage);
        return false;
    }

    // TODO: Make use of this really useful index...
    const int specieIndex = lastResult->GetSpecieIndex(CMolecule(m_userSettings.m_molecule).m_name);

    // 9. Get the mode of the evaluation
    lastResult->CheckMeasurementMode();
    lastResult->GetStartTime(0, startTime);

    // 10. Append the results to the evaluation-summary log
    AppendToEvaluationSummaryFile(lastResult, &scan, &instrLocation, &fitWindow, windField);
    AppendToPakFileSummaryFile(lastResult, &scan, &instrLocation, &fitWindow, windField);

    // 10. Append the result to the log file of the corresponding scanningInstrument
    if (RETURN_CODE::SUCCESS != WriteEvaluationResult(lastResult, &scan, &instrLocation, &fitWindow, windField, txtFileName))
    {
        errorMessage.Format("Failed to write evaluation log file %s. No result produced", txtFileName);
        ShowMessage(errorMessage);
    }

    // 11. If this was a flux-measurement then we need to see the plume for the measurement to be useful
    //  this check should only be performed on the main fit window.
    if (Equals(fitWindow.name, m_userSettings.m_fitWindowsToUse[m_userSettings.m_mainFitWindow]))
    {
        if (0 == CheckQualityOfFluxMeasurement(lastResult, pakFileName))
        {
            errorMessage.Format("Could not calculate flux from pak-file %s.", (const char*)pakFileName);
            m_log.Information(errorMessage.std_str());
            return false;
        }
    }

    // 12. Return the properties of the scan
    if (plumeProperties != nullptr)
    {
        lastResult->GetCalculatedPlumeProperties(*plumeProperties);
    }

    CreatePlumespectrumFile(lastResult, fitWindowName, scan, spectrometerModel, plumeProperties, specieIndex);

    return true;
}

void CPostEvaluationController::CreatePlumespectrumFile(
    const std::unique_ptr<CScanResult>& result,
    const novac::CString& fitWindowName,
    novac::CScanFileHandler& scan,
    const novac::SpectrometerModel& spectrometerModel,
    novac::CPlumeInScanProperty* plumeProperties,
    int specieIndex)
{
    const std::string outputDirectoryStr =
        std::string(m_userSettings.m_outputDirectory) +
        "/" + std::string(fitWindowName) +
        "/PlumeSpectra/" +
        result->GetSerial().std_str();

    int ret = Filesystem::CreateDirectoryStructure(outputDirectoryStr);
    if (ret)
    {
        novac::CString userMessage;
        userMessage.Format("Could not create directory for archiving plume spectra: %s", outputDirectoryStr.c_str());
        ShowMessage(userMessage);
    }
    else
    {
        Configuration::RatioEvaluationSettings plumeCalculationSettings;

        PlumeSpectrumSelector spectrumSelector;
        spectrumSelector.CreatePlumeSpectrumFile(
            scan,
            *result,
            *plumeProperties,
            plumeCalculationSettings,
            spectrometerModel,
            specieIndex,
            outputDirectoryStr);
    }
}

RETURN_CODE CPostEvaluationController::WriteRatioResult(const std::vector<Ratio>& result, const novac::CScanFileHandler& scan, const novac::CFitWindow& window)
{
    CSpectrum skySpec;
    scan.GetSky(skySpec);

    novac::CString fileName;
    fileName.Format("%s%c%s%cRatioSummary_%s.txt",
        (const char*)m_userSettings.m_outputDirectory,
        Poco::Path::separator(),
        window.name.c_str(),
        Poco::Path::separator(),
        skySpec.m_info.m_device.c_str());

    const bool writeHeaderLine = !Filesystem::IsExistingFile(fileName);

    FILE* f = fopen(fileName, "a");
    if (f == nullptr)
    {
        return RETURN_CODE::FAIL;
    }

    if (writeHeaderLine)
    {
        fprintf(f, "StartTime\t");
        for (const Ratio& r : result)
        {
            fprintf(f, "Ratio\tError\t");
            fprintf(f, "Column(%s)\tError(%s)\t", r.minorSpecieName.c_str(), r.minorSpecieName.c_str());
            fprintf(f, "Column(%s)\tError(%s)\t", r.majorSpecieName.c_str(), r.majorSpecieName.c_str());
        }
        fprintf(f, "\n");
    }

    // The start-time
    fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t", scan.m_startTime.year, scan.m_startTime.month, scan.m_startTime.day, scan.m_startTime.hour, scan.m_startTime.minute, scan.m_startTime.second);

    // The major specie
    for (const Ratio& r : result)
    {
        fprintf(f, "%.5e\t%.5e\t", r.ratio, r.error);
        fprintf(f, "%.5e\t%.5e\t", r.minorResult, r.minorError);
        fprintf(f, "%.5e\t%.5e\t", r.majorResult, r.majorError);
    }
    fprintf(f, "\n");

    fclose(f);

    return RETURN_CODE::SUCCESS;
}

RETURN_CODE CPostEvaluationController::WriteEvaluationResult(const std::unique_ptr<CScanResult>& result, const novac::CScanFileHandler* scan, const Configuration::CInstrumentLocation* instrLocation, const novac::CFitWindow* window, Meteorology::CWindField& windField, novac::CString* txtFileName)
{
    novac::CString string, string1, string2, string3, string4;
    long itSpectrum, itSpecie; // iterators
    novac::CString pakFile, txtFile, evalSummaryLog;
    novac::CString wsSrc, wdSrc, phSrc;
    CDateTime dateTime;

    // get the file-name that we want to have 
    GetArchivingfileName(pakFile, txtFile, window->name, scan->GetFileName(), m_userSettings.m_outputDirectory.std_str(), result->GetMeasurementMode());
    if (txtFileName != nullptr)
    {
        txtFileName->Format(txtFile);
    }

    // the spectrometer 
    const SpectrometerModel spectrometerModel = CSpectrometerDatabase::GetInstance().GetModel(instrLocation->m_spectrometerModel);

    // The date of the measurement & the serial-number of the spectrometer
    result->GetSkyStartTime(dateTime);

    // 0. Create the additional scan-information
    string.Format("\n<scaninformation>\n");
    string.AppendFormat("\tdate=%02d.%02d.%04d\n", scan->m_startTime.day, scan->m_startTime.month, scan->m_startTime.year);
    string.AppendFormat("\tstarttime=%02d:%02d:%02d\n", scan->m_startTime.hour, scan->m_startTime.minute, scan->m_startTime.second);
    string.AppendFormat("\tcompass=%.1lf\n", instrLocation->m_compass);
    string.AppendFormat("\ttilt=%.1lf\n", instrLocation->m_tilt);
    string.AppendFormat("\tlat=%.6lf\n", instrLocation->m_latitude);
    string.AppendFormat("\tlong=%.6lf\n", instrLocation->m_longitude);
    string.AppendFormat("\talt=%d\n", instrLocation->m_altitude);

    string.AppendFormat("\tvolcano=%s\n", (const char*)instrLocation->m_volcano);
    string.AppendFormat("\tsite=%s\n", (const char*)instrLocation->m_locationName);
    // string.AppendFormat("\tobservatory=%s\n",       m_common.SimplifyString(spectrometer.m_scanner.observatory));

    string.AppendFormat("\tserial=%s\n", (const char*)result->GetSerial());
    string.AppendFormat("\tspectrometer=%s\n", instrLocation->m_spectrometerModel.c_str());
    string.AppendFormat("\tspectrometerMaxIntensity=%lf\n", spectrometerModel.maximumIntensityForSingleReadout);

    string.AppendFormat("\tchannel=%d\n", window->channel);
    string.AppendFormat("\tconeangle=%.1lf\n", instrLocation->m_coneangle);
    string.AppendFormat("\tinterlacesteps=%d\n", scan->GetInterlaceSteps());
    string.AppendFormat("\tstartchannel=%d\n", scan->GetStartChannel());
    string.AppendFormat("\tspectrumlength=%d\n", scan->GetSpectrumLength());
    string.AppendFormat("\tflux=%.2lf\n", result->GetFlux());
    string.AppendFormat("\tbattery=%.2f\n", result->GetBatteryVoltage());
    string.AppendFormat("\ttemperature=%.2f\n", result->GetTemperature());

    // The mode
    if (result->IsDirectSunMeasurement())
        string.Append("\tmode=direct_sun\n");
    else if (result->IsLunarMeasurement())
        string.Append("\tmode=lunar\n");
    else if (result->IsWindMeasurement())
        string.Append("\tmode=wind\n");
    else if (result->IsStratosphereMeasurement())
        string.Append("\tmode=stratospheric\n");
    else if (result->IsCompositionMeasurement())
        string.Append("\tmode=composition\n");
    else
        string.Append("\tmode=plume\n");

    // The type of instrument used...
    if (instrLocation->m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG)
    {
        string.Append("\tinstrumenttype=gothenburg\n");
    }
    else if (instrLocation->m_instrumentType == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
    {
        string.Append("\tinstrumenttype=heidelberg\n");
    }

    // Finally, the version of the file and the version of the program
    string.Append("\tversion=2.2\n");
    string.Append("\tsoftware=NovacPPP\n");
    string.AppendFormat("\tcompiledate=%s\n", __DATE__);

    string.Append("</scaninformation>\n");
    // 0a. Write the additional scan-information to the evaluation log
    FILE* f = fopen(txtFile, "w");
    if (f != nullptr)
    {
        fprintf(f, "%s", string.c_str());
        fprintf(f, "\n");
    }

    // 0.1 Create an flux-information part and write it to the same file
    windField.GetWindSpeedSource(wsSrc);
    windField.GetWindDirectionSource(wdSrc);
    // windField.GetPlumeHeightSource(phSrc);

        // Get the information on where the plume is seen
    double plumeEdge1, plumeEdge2;
    double plumeCompleteness = result->GetCalculatedPlumeCompleteness();
    double plumeCentre1 = result->GetCalculatedPlumeCentre(0);
    double plumeCentre2 = result->GetCalculatedPlumeCentre(1);
    result->GetCalculatedPlumeEdges(plumeEdge1, plumeEdge2);

    string.Format("<fluxinfo>\n");
    string.AppendFormat("\tflux=%.4lf\n", result->GetFlux()); // ton/day
    string.AppendFormat("\twindspeed=%.4lf\n", windField.GetWindSpeed());
    string.AppendFormat("\twinddirection=%.4lf\n", windField.GetWindDirection());
    // string.AppendFormat("\tplumeheight=%.2lf\n",  windField.GetPlumeHeight());
    string.AppendFormat("\twindspeedsource=%s\n", (const char*)wsSrc);
    string.AppendFormat("\twinddirectionsource=%s\n", (const char*)wdSrc);
    string.AppendFormat("\tplumeheightsource=%s\n", (const char*)phSrc);
    if (fabs(instrLocation->m_compass) > 360.0)
        string.Append("\tcompasssource=compassreading\n");
    else
        string.Append("\tcompasssource=user\n");

    string.AppendFormat("\tplumecompleteness=%.2lf\n", plumeCompleteness);
    string.AppendFormat("\tplumecentre=%.2lf\n", plumeCentre1);
    if (instrLocation->m_instrumentType == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
        string.AppendFormat("\tplumecentre_phi=%.2lf\n", plumeCentre2);
    string.AppendFormat("\tplumeedge1=%.2lf\n", plumeEdge1);
    string.AppendFormat("\tplumeedge2=%.2lf\n", plumeEdge2);

    string.Append("</fluxinfo>");

    // 0.1b Write the flux-information to the evaluation-log
    if (f != nullptr)
    {
        fprintf(f, "%s", string.c_str());
        fprintf(f, "\n");
    }


    // 1. write the header
    if (instrLocation->m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG)
    {
        string.Format("#scanangle\t");
    }
    else if (instrLocation->m_instrumentType == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
    {
        string.Format("#observationangle\tazimuth\t");
    }
    string.Append("starttime\tstoptime\tname\tspecsaturation\tfitsaturation\tcounts_ms\tdelta\tchisquare\texposuretime\tnumspec\t");

    for (itSpecie = 0; itSpecie < window->nRef; ++itSpecie)
    {
        string.AppendFormat("column(%s)\tcolumnerror(%s)\t", window->ref[itSpecie].m_specieName.c_str(), window->ref[itSpecie].m_specieName.c_str());
        string.AppendFormat("shift(%s)\tshifterror(%s)\t", window->ref[itSpecie].m_specieName.c_str(), window->ref[itSpecie].m_specieName.c_str());
        string.AppendFormat("squeeze(%s)\tsqueezeerror(%s)\t", window->ref[itSpecie].m_specieName.c_str(), window->ref[itSpecie].m_specieName.c_str());
    }
    string.Append("isgoodpoint\toffset\tflag");

    // 1a. Write the header to the log file
    if (f != nullptr)
    {
        fprintf(f, "%s", string.c_str());
        fprintf(f, "\n<spectraldata>\n");
    }

    // ----------------------------------------------------------------------------------------------
    // 2. ----------------- Write the parameters for the sky and the dark-spectra -------------------
    // ----------------------------------------------------------------------------------------------
    CSpectrum sky, dark, darkCurrent, offset;
    string1.Format(""); string2.Format(""); string3.Format(""); string4.Format("");
    scan->GetSky(sky);
    if (sky.m_info.m_interlaceStep > 1)
        sky.InterpolateSpectrum();
    if (sky.m_length > 0)
    {
        sky.m_info.m_fitIntensity = (float)(sky.MaxValue(window->fitLow, window->fitHigh));
        if (sky.NumSpectra() > 0)
            sky.Div(sky.NumSpectra());
        CEvaluationLogFileHandler::FormatEvaluationResult(&sky.m_info, nullptr, instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * sky.NumSpectra(), window->nRef, string1);
    }
    scan->GetDark(dark);
    if (dark.m_info.m_interlaceStep > 1)
        dark.InterpolateSpectrum();
    if (dark.m_length > 0)
    {
        dark.m_info.m_fitIntensity = (float)(dark.MaxValue(window->fitLow, window->fitHigh));
        if (dark.NumSpectra() > 0)
            dark.Div(dark.NumSpectra());
        CEvaluationLogFileHandler::FormatEvaluationResult(&dark.m_info, nullptr, instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * dark.NumSpectra(), window->nRef, string2);
    }
    scan->GetOffset(offset);
    if (offset.m_info.m_interlaceStep > 1)
        offset.InterpolateSpectrum();
    if (offset.m_length > 0)
    {
        offset.m_info.m_fitIntensity = (float)(offset.MaxValue(window->fitLow, window->fitHigh));
        offset.Div(offset.NumSpectra());
        CEvaluationLogFileHandler::FormatEvaluationResult(&offset.m_info, nullptr, instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * offset.NumSpectra(), window->nRef, string3);
    }
    scan->GetDarkCurrent(darkCurrent);
    if (darkCurrent.m_info.m_interlaceStep > 1)
        darkCurrent.InterpolateSpectrum();
    if (darkCurrent.m_length > 0)
    {
        darkCurrent.m_info.m_fitIntensity = (float)(darkCurrent.MaxValue(window->fitLow, window->fitHigh));
        darkCurrent.Div(darkCurrent.NumSpectra());
        CEvaluationLogFileHandler::FormatEvaluationResult(&darkCurrent.m_info, nullptr, instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * darkCurrent.NumSpectra(), window->nRef, string4);
    }

    // 2b. Write it all to the evaluation log file
    if (f != nullptr)
    {
        if (strlen(string1) > 0)
        {
            fprintf(f, "%s", string1.c_str()); fprintf(f, "\n");
        }
        if (strlen(string2) > 0)
        {
            fprintf(f, "%s", string2.c_str()); fprintf(f, "\n");
        }
        if (strlen(string3) > 0)
        {
            fprintf(f, "%s", string3.c_str()); fprintf(f, "\n");
        }
        if (strlen(string4) > 0)
        {
            fprintf(f, "%s", string4.c_str()); fprintf(f, "\n");
        }
    }


    // ----------------------------------------------------------------------------------------------
    // 3. ------------------- Then write the parameters for each spectrum ---------------------------
    // ----------------------------------------------------------------------------------------------
    for (itSpectrum = 0; itSpectrum < result->GetEvaluatedNum(); ++itSpectrum)
    {
        int nSpectra = result->GetSpectrumInfo(itSpectrum).m_numSpec;

        // 3a. Pretty print the result and the spectral info into a string
        CEvaluationLogFileHandler::FormatEvaluationResult(&result->GetSpectrumInfo(itSpectrum), result->GetResult(itSpectrum), instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * nSpectra, window->nRef, string);

        // 3b. Write it all to the evaluation log file
        if (f != nullptr)
        {
            fprintf(f, "%s", string.c_str());
            fprintf(f, "\n");
        }
    }

    if (f != nullptr)
    {
        fprintf(f, "</spectraldata>\n");
        fclose(f);
    }

    return RETURN_CODE::SUCCESS;
}

RETURN_CODE CPostEvaluationController::AppendToEvaluationSummaryFile(const std::unique_ptr<CScanResult>& result, const novac::CScanFileHandler* scan, const Configuration::CInstrumentLocation* /*instrLocation*/, const novac::CFitWindow* window, Meteorology::CWindField& /*windField*/)
{
    novac::CString evalSummaryLog;
    bool fWriteHeaderLine = false;

    // we can also write an evaluation-summary log file
    evalSummaryLog.Format("%s%c%s%cEvaluationSummary_%s.txt",
        (const char*)m_userSettings.m_outputDirectory,
        Poco::Path::separator(),
        window->name.c_str(),
        Poco::Path::separator(),
        (const char*)result->GetSerial());

    if (!Filesystem::IsExistingFile(evalSummaryLog))
    {
        fWriteHeaderLine = true;
    }

    FILE* f = fopen(evalSummaryLog, "a");
    if (f == nullptr)
        return RETURN_CODE::FAIL;

    if (fWriteHeaderLine)
    {
        fprintf(f, "StartTime\tExpTime\tAppliedShift\tTemperature\tCalculatedOffset\tCalculatedPlumeCentre\tCalculatedPlumeCompleteness\t#Spectra\n");
    }

    // the start-time
    fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t", scan->m_startTime.year, scan->m_startTime.month, scan->m_startTime.day, scan->m_startTime.hour, scan->m_startTime.minute, scan->m_startTime.second);

    // The exposure time
    fprintf(f, "%ld\t", result->GetSkySpectrumInfo().m_exposureTime);

    // the shift applied
    fprintf(f, "%.2lf\t", result->GetResult(0)->m_referenceResult[0].m_shift);

    // the temperature of the spectrometer
    fprintf(f, "%.2lf\t", result->GetTemperature());

    // the calculated plume parameters
    fprintf(f, "%.2lf\t", result->GetOffset());
    fprintf(f, "%.2lf\t", result->GetCalculatedPlumeCentre());
    fprintf(f, "%.2lf\t", result->GetCalculatedPlumeCompleteness());

    // the number of evaluated spectra
    fprintf(f, "%ld\t", result->GetEvaluatedNum());

    // make a new line
    fprintf(f, "\n");

    // close the file
    fclose(f);

    return RETURN_CODE::SUCCESS;
}

RETURN_CODE CPostEvaluationController::AppendToPakFileSummaryFile(const std::unique_ptr<CScanResult>& result, const novac::CScanFileHandler* scan, const Configuration::CInstrumentLocation* /*instrLocation*/, const novac::CFitWindow* /*window*/, Meteorology::CWindField& /*windField*/)
{
    novac::CString pakSummaryLog;
    bool fWriteHeaderLine = false;

    // we can also write an evaluation-summary log file
    pakSummaryLog.Format("%s%cPakfileSummary.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());

    if (!Filesystem::IsExistingFile(pakSummaryLog))
    {
        fWriteHeaderLine = true;
    }

    FILE* f = fopen(pakSummaryLog, "a");
    if (f == nullptr)
    {
        return RETURN_CODE::FAIL;
    }

    if (fWriteHeaderLine)
    {
        fprintf(f, "Serial\tStartTime\tLat\tLong\tAlt\tExpTime\tBatteryVoltage\tTemperature\tElectronicOffset\n");
    }

    // the serial of the instrument
    fprintf(f, "%s\t", (const char*)result->GetSerial());

    // the start-time
    fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t", scan->m_startTime.year, scan->m_startTime.month, scan->m_startTime.day, scan->m_startTime.hour, scan->m_startTime.minute, scan->m_startTime.second);

    // the location
    const CGPSData& gps = scan->GetGPS();
    fprintf(f, "%.5lf\t%.5lf\t%.5lf\t", gps.m_latitude, gps.m_longitude, gps.m_altitude);

    // The exposure time
    fprintf(f, "%ld\t", result->GetSkySpectrumInfo().m_exposureTime);

    // the input-voltage at the time of measurement
    fprintf(f, "%.2lf\t", result->GetBatteryVoltage());

    // the temperature of the spectrometer
    fprintf(f, "%.2lf\t", result->GetTemperature());

    // the offset of the AD converter
    fprintf(f, "%.2lf", result->GetElectronicOffset(0));

    // make a new line
    fprintf(f, "\n");

    // close the file
    fclose(f);

    return RETURN_CODE::SUCCESS;
}

int CPostEvaluationController::CheckQualityOfFluxMeasurement(std::unique_ptr<CScanResult>& result, const novac::CString& pakFileName) const
{
    novac::CString errorMessage;

    // Check if this is a flux measurement at all
    if (MEASUREMENT_MODE::MODE_FLUX != result->GetMeasurementMode())
    {
        return -1;
    }
    if (0 == result->CalculateOffset(CMolecule(m_userSettings.m_molecule)))
    {
        if (0 == result->CalculatePlumeCentre(CMolecule(m_userSettings.m_molecule)))
        {
            // no plume found!
            errorMessage.Format(" - Scan %s does not see the plume. Scan ignored.", (const char*)pakFileName);
            ShowMessage(errorMessage);
            return 0;
        }
    }
    return 1;
}




