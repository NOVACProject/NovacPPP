#include <PPPLib/Evaluation/EvaluationUtils.h>
#include <PPPLib/Logging.h>
#include <PPPLib/File/Filesystem.h>

#include <SpectralEvaluation/File/SpectrumIO.h>

#include <sstream>

#ifndef MAX_PATH
#define MAX_PATH 512
#endif MAX_PATH

bool Evaluation::IsGoodEnoughToEvaluate(
    const novac::CScanFileHandler& scan,
    const novac::CFitWindow& fitWindow,
    const novac::SpectrometerModel& model,
    const Configuration::CInstrumentLocation& instrLocation,
    const Configuration::CUserConfiguration& userSettings,
    ReasonForScanRejection& reason,
    std::string& reasonMessage)
{
    novac::CSpectrum skySpectrum;

    // Check that the sky-spectrum is ok
    scan.GetSky(skySpectrum);
    if (skySpectrum.IsDark())
    {
        reasonMessage = "Sky spectrum is dark";
        reason = ReasonForScanRejection::SkySpectrumDark;
        return false;
    }

    if ((instrLocation.m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG && skySpectrum.ExposureTime() > userSettings.m_maxExposureTime_got) ||
        (instrLocation.m_instrumentType == INSTRUMENT_TYPE::INSTR_HEIDELBERG && skySpectrum.ExposureTime() > userSettings.m_maxExposureTime_hei))
    {
        std::stringstream msg;
        msg << "Sky spectrum has too long exposure time (" << skySpectrum.ExposureTime() << " ms)";
        reasonMessage = msg.str();
        reason = ReasonForScanRejection::SkySpectrumTooLongExposureTime;
        return false;
    }

        const double dynamicRange = skySpectrum.NumSpectra() * model.maximumIntensityForSingleReadout;

    if (skySpectrum.MaxValue(fitWindow.fitLow, fitWindow.fitHigh) >= dynamicRange)
    {
        reasonMessage = "Sky spectrum is saturated in fit region";
        reason = ReasonForScanRejection::SkySpectrumSaturated;
        return false;
    }

    return true;
}


bool GetArchivingfileName(
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
        ShowMessage(userMessage);
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
