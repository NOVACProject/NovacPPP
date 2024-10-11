#include <PPPLib/Evaluation/PostEvaluationIO.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/File/EvaluationLogFileHandler.h>

#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/Configuration/RatioEvaluationSettings.h>
#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>

#include <cstring>

#ifndef MAX_PATH
#define MAX_PATH 512
#endif MAX_PATH

bool Evaluation::PostEvaluationIO::GetArchivingfileName(
    novac::ILogger& log,
    novac::CString& pakFile,
    novac::CString& txtFile,
    const novac::CString& fitWindowName,
    const novac::CString& temporaryScanFile,
    std::string outputDirectory,
    MEASUREMENT_MODE mode)
{
    novac::CSpectrumIO reader;
    novac::CSpectrum tmpSpec;
    novac::CString serialNumber, dateStr, timeStr, dateStr2, modeStr, userMessage;

    const char pathSeparator = '/';

    // 0. Make an initial assumption of the file-names
    int i = 0;
    while (1)
    {
        pakFile.Format("%s%cUnknownScans%c%d.pak", outputDirectory.c_str(), pathSeparator, pathSeparator, ++i);
        if (!Filesystem::IsExistingFile(pakFile))
        {
            break;
        }
    }
    txtFile.Format("%s%cUnknownScans%c%d.txt", outputDirectory.c_str(), pathSeparator, pathSeparator, i);

    // 1. Read the first spectrum in the scan
    const std::string temporaryScanFileStr((const char*)temporaryScanFile);
    if (!reader.ReadSpectrum(temporaryScanFileStr, 0, tmpSpec))
    {
        return false;
    }
    novac::CSpectrumInfo& info = tmpSpec.m_info;
    int channel = info.m_channel;

    // 1a. If the GPS had no connection with the satelites when collecting the sky-spectrum,
    //   then try to find a spectrum in the file for which it had connection...
    i = 1;
    while (info.m_startTime.year == 2004 && info.m_startTime.month == 3 && info.m_startTime.second == 22)
    {
        if (!reader.ReadSpectrum(temporaryScanFileStr, i++, tmpSpec))
        {
            break;
        }
        info = tmpSpec.m_info;
    }

    // 2. Get the serialNumber of the spectrometer
    serialNumber.Format("%s", info.m_device.c_str());

    // 3. Get the time and date when the scan started
    dateStr.Format("%02d%02d%02d", info.m_startTime.year % 1000, info.m_startTime.month, info.m_startTime.day);
    dateStr2.Format("%04d.%02d.%02d", info.m_startTime.year, info.m_startTime.month, info.m_startTime.day);
    timeStr.Format("%02d%02d", info.m_startTime.hour, info.m_startTime.minute);


    // 4. Write the archiving name of the spectrum file

    // 4a. Write the folder name
    pakFile.Format("%s%s%c%s%c%s%c", outputDirectory.c_str(), (const char*)fitWindowName, pathSeparator,
        (const char*)dateStr2, pathSeparator, (const char*)serialNumber, pathSeparator);
    txtFile.Format("%s", (const char*)pakFile);

    // 4b. Make sure that the folder exists
    int ret = Filesystem::CreateDirectoryStructure(pakFile);
    if (ret)
    {
        userMessage.Format("Could not create directory for archiving .pak-file: %s", (const char*)pakFile);
        log.Error(userMessage.std_str());
        return false;
    }

    // 4c. Write the code for the measurement mode
    switch (mode)
    {
    case MEASUREMENT_MODE::MODE_FLUX:   modeStr.Format("flux"); break;
    case MEASUREMENT_MODE::MODE_WINDSPEED: modeStr.Format("wind"); break;
    case MEASUREMENT_MODE::MODE_STRATOSPHERE: modeStr.Format("stra"); break;
    case MEASUREMENT_MODE::MODE_DIRECT_SUN: modeStr.Format("dsun"); break;
    case MEASUREMENT_MODE::MODE_COMPOSITION:  modeStr.Format("comp"); break;
    case MEASUREMENT_MODE::MODE_LUNAR:  modeStr.Format("luna"); break;
    case MEASUREMENT_MODE::MODE_TROPOSPHERE:  modeStr.Format("trop"); break;
    case MEASUREMENT_MODE::MODE_MAXDOAS:  modeStr.Format("maxd"); break;
    default:    modeStr.Format("unkn"); break;
    }

    // 4c. Write the name of the archiving file itself
    if (channel < 128 && channel > MAX_CHANNEL_NUM)
    {
        channel = channel % 16;
    }

    pakFile.AppendFormat("%s_%s_%s_%1d_%4s.pak", (const char*)serialNumber, (const char*)dateStr, (const char*)timeStr, channel, (const char*)modeStr);
    txtFile.AppendFormat("%s_%s_%s_%1d_%4s.txt", (const char*)serialNumber, (const char*)dateStr, (const char*)timeStr, channel, (const char*)modeStr);

    if (strlen(pakFile) > MAX_PATH - 2)
    {
        return false;
    }

    return true;
}

RETURN_CODE Evaluation::PostEvaluationIO::WriteEvaluationResult(
    novac::ILogger& log,
    const std::string& outputDirectory,
    novac::SpectrometerModel spectrometerModel,
    const std::unique_ptr<CScanResult>& result,
    const novac::CScanFileHandler* scan,
    const Configuration::CInstrumentLocation* instrLocation,
    const novac::CFitWindow* window,
    Meteorology::CWindField& windField,
    novac::CString* txtFileName)
{
    novac::CString string, string1, string2, string3, string4;
    long itSpectrum, itSpecie; // iterators
    novac::CString pakFile, txtFile, evalSummaryLog;
    novac::CString wsSrc, wdSrc, phSrc;
    novac::CDateTime dateTime;

    // get the file-name that we want to have 
    GetArchivingfileName(log, pakFile, txtFile, window->name, scan->GetFileName(), outputDirectory, result->GetMeasurementMode());
    if (txtFileName != nullptr)
    {
        txtFileName->Format(txtFile);
    }

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
    if (std::abs(instrLocation->m_compass) > 360.0)
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
    novac::CSpectrum sky, dark, darkCurrent, offset;
    string1.Format(""); string2.Format(""); string3.Format(""); string4.Format("");
    scan->GetSky(sky);
    if (sky.m_info.m_interlaceStep > 1)
        sky.InterpolateSpectrum();
    if (sky.m_length > 0)
    {
        sky.m_info.m_fitIntensity = (float)(sky.MaxValue(window->fitLow, window->fitHigh));
        if (sky.NumSpectra() > 0)
            sky.Div(sky.NumSpectra());
        FileHandler::CEvaluationLogFileHandler::FormatEvaluationResult(&sky.m_info, nullptr, instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * sky.NumSpectra(), window->nRef, string1);
    }
    scan->GetDark(dark);
    if (dark.m_info.m_interlaceStep > 1)
        dark.InterpolateSpectrum();
    if (dark.m_length > 0)
    {
        dark.m_info.m_fitIntensity = (float)(dark.MaxValue(window->fitLow, window->fitHigh));
        if (dark.NumSpectra() > 0)
            dark.Div(dark.NumSpectra());
        FileHandler::CEvaluationLogFileHandler::FormatEvaluationResult(&dark.m_info, nullptr, instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * dark.NumSpectra(), window->nRef, string2);
    }
    scan->GetOffset(offset);
    if (offset.m_info.m_interlaceStep > 1)
        offset.InterpolateSpectrum();
    if (offset.m_length > 0)
    {
        offset.m_info.m_fitIntensity = (float)(offset.MaxValue(window->fitLow, window->fitHigh));
        offset.Div(offset.NumSpectra());
        FileHandler::CEvaluationLogFileHandler::FormatEvaluationResult(&offset.m_info, nullptr, instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * offset.NumSpectra(), window->nRef, string3);
    }
    scan->GetDarkCurrent(darkCurrent);
    if (darkCurrent.m_info.m_interlaceStep > 1)
        darkCurrent.InterpolateSpectrum();
    if (darkCurrent.m_length > 0)
    {
        darkCurrent.m_info.m_fitIntensity = (float)(darkCurrent.MaxValue(window->fitLow, window->fitHigh));
        darkCurrent.Div(darkCurrent.NumSpectra());
        FileHandler::CEvaluationLogFileHandler::FormatEvaluationResult(&darkCurrent.m_info, nullptr, instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * darkCurrent.NumSpectra(), window->nRef, string4);
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
        FileHandler::CEvaluationLogFileHandler::FormatEvaluationResult(&result->GetSpectrumInfo(itSpectrum), result->GetResult(itSpectrum), instrLocation->m_instrumentType, spectrometerModel.maximumIntensityForSingleReadout * nSpectra, window->nRef, string);

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

void Evaluation::PostEvaluationIO::CreatePlumespectrumFile(
    novac::ILogger& log,
    const std::string& outputDirectory,
    const std::unique_ptr<Evaluation::CScanResult>& result,
    const novac::CString& fitWindowName,
    novac::CScanFileHandler& scan,
    const novac::SpectrometerModel& spectrometerModel,
    novac::CPlumeInScanProperty* plumeProperties,
    int specieIndex)
{
    const std::string outputDirectoryStr =
        outputDirectory +
        "/" + std::string(fitWindowName) +
        "/PlumeSpectra/" +
        result->GetSerial().std_str();

    int ret = Filesystem::CreateDirectoryStructure(outputDirectoryStr);
    if (ret)
    {
        novac::CString userMessage;
        userMessage.Format("Could not create directory for archiving plume spectra: %s", outputDirectoryStr.c_str());
        log.Error(userMessage.std_str());
    }
    else
    {
        Configuration::RatioEvaluationSettings plumeCalculationSettings;

        novac::PlumeSpectrumSelector spectrumSelector;
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

RETURN_CODE Evaluation::PostEvaluationIO::AppendToEvaluationSummaryFile(
    const std::string& outputDirectory,
    const std::unique_ptr<CScanResult>& result,
    const novac::CScanFileHandler* scan,
    const Configuration::CInstrumentLocation* /*instrLocation*/,
    const novac::CFitWindow* window,
    Meteorology::CWindField& /*windField*/)
{
    novac::CString evalSummaryLog;
    bool fWriteHeaderLine = false;

    // we can also write an evaluation-summary log file
    evalSummaryLog.Format("%s/%s/EvaluationSummary_%s.txt",
        outputDirectory.c_str(),
        window->name.c_str(),
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

RETURN_CODE Evaluation::PostEvaluationIO::AppendToPakFileSummaryFile(
    const std::string& outputDirectory,
    const std::unique_ptr<Evaluation::CScanResult>& result,
    const novac::CScanFileHandler* scan,
    const Configuration::CInstrumentLocation* /*instrLocation*/,
    const novac::CFitWindow* /*window*/,
    Meteorology::CWindField& /*windField*/)
{
    novac::CString pakSummaryLog;
    bool fWriteHeaderLine = false;

    // we can also write an evaluation-summary log file
    pakSummaryLog.Format("%s/PakfileSummary.txt", outputDirectory.c_str());

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
    const novac::CGPSData& gps = scan->GetGPS();
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

