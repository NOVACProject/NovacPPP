#include <PPPLib/Calibration/PostCalibration.h>
#include <PPPLib/Calibration/PostCalibrationStatistics.h>
#include <PPPLib/PPPLib.h>
#include <PPPLib/Logging.h>
#include <PPPLib/File/Filesystem.h>
#include <SpectralEvaluation/DialogControllers/NovacProgramWavelengthCalibrationController.h>
#include <SpectralEvaluation/DialogControllers/ReferenceCreationController.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Calibration/StandardCrossSectionSetup.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/File/EvaluationConfigurationParser.h>
#include <sstream>
#include <algorithm>
#include <PPPLib/MFC/CFileUtils.h>

using namespace novac;

std::string FormatDateAndTimeOfSpectrum(const novac::CSpectrumInfo& spectrumInformation)
{
    CString dateAndTime;

    dateAndTime.Format(
        "%02d%02d%02d_%02d%02d",
        spectrumInformation.m_startTime.year % 1000,
        spectrumInformation.m_startTime.month,
        spectrumInformation.m_startTime.day,
        spectrumInformation.m_startTime.hour,
        spectrumInformation.m_startTime.minute);

    return std::string(dateAndTime);
}

std::string CreateOutputDirectoryForCalibration(const CSpectrumInfo& calibratedSpectrum, const Configuration::CUserConfiguration& userSettings)
{
    CString dateStr;
    dateStr.Format(
        "%04d.%02d.%02d",
        calibratedSpectrum.m_startTime.year,
        calibratedSpectrum.m_startTime.month,
        calibratedSpectrum.m_startTime.day);

    std::string directoryName{ (const char*)userSettings.m_outputDirectory };

    if (directoryName.back() != '/' && directoryName.back() != '\\')
    {
        directoryName += '/';
    }
    directoryName += dateStr + "/" + calibratedSpectrum.m_device + "/";

    // 4b. Make sure that the folder exists
    int ret = Filesystem::CreateDirectoryStructure(directoryName);
    if (ret)
    {
        std::stringstream message;
        message << "Could not create directory for archiving instrument calibration: " << directoryName;
        ShowMessage(message.str());

        // TODO: Another type of exception!
        throw std::invalid_argument(message.str());
    }

    return directoryName;
}

std::string GetCalibrationFileName(const novac::CSpectrumInfo& spectrumInformation)
{
    return "Calibration_" + spectrumInformation.m_device + "_" + FormatDateAndTimeOfSpectrum(spectrumInformation) + ".std";
}

void RunCalibration(
    NovacProgramWavelengthCalibrationController& calibrationController,
    const std::string& scanFile,
    const Configuration::CUserConfiguration& userSettings,
    const Configuration::CInstrumentCalibrationConfiguration& calibrationSettings)
{
    calibrationController.m_inputSpectrumFile = scanFile;
    calibrationController.m_solarSpectrumFile = userSettings.m_highResolutionSolarSpectrumFile;
    calibrationController.m_initialCalibrationFile = calibrationSettings.m_initialCalibrationFile;
    calibrationController.m_initialLineShapeFile = calibrationSettings.m_instrumentLineshapeFile;
    calibrationController.m_instrumentLineShapeFitOption = (WavelengthCalibrationController::InstrumentLineShapeFitOption)userSettings.m_calibrationInstrumentLineShapeFitOption;
    calibrationController.m_instrumentLineShapeFitRegion = userSettings.m_calibrationInstrumentLineShapeFitRegion;

    // Does the actual calibration. Throws a std::exception if the calibration fails.
    calibrationController.RunCalibration();
}

std::vector<novac::CReferenceFile> CreateStandardReferences(
    const novac::CSpectrumInfo& spectrumInformation,
    const std::unique_ptr<novac::InstrumentCalibration>& calibration,
    const novac::StandardCrossSectionSetup standardCrossSections,
    const std::string& directoryName)
{
    ReferenceCreationController referenceController;
    std::vector<novac::CReferenceFile> referencesCreated;

    referenceController.m_highPassFilter = false; // In the NovacPPP the references are never filtered on disk (but may be filtered after having been read in)
    referenceController.m_unitSelection = 1; // Default to molecules/cm2 in NovacPPP

    // First the ordinary references
    for (size_t ii = 0; ii < standardCrossSections.NumberOfReferences(); ++ii)
    {
        referenceController.m_convertToAir = standardCrossSections.IsReferenceInVacuum(ii);
        referenceController.m_highResolutionCrossSection = standardCrossSections.ReferenceFileName(ii);
        referenceController.m_isPseudoAbsorber = standardCrossSections.IsAdditionalAbsorber(ii);
        referenceController.ConvolveReference(*calibration);

        // Save the result
        const std::string filteringStr = (referenceController.m_highPassFilter) ?
            "_HP500_PPMM" :
            "";
        const std::string dstFileName =
            directoryName +
            spectrumInformation.m_device +
            "_" +
            standardCrossSections.ReferenceSpecieName(ii) +
            filteringStr +
            "_" +
            FormatDateAndTimeOfSpectrum(spectrumInformation) +
            ".txt";
        novac::SaveCrossSectionFile(dstFileName, *(referenceController.m_resultingCrossSection));

        novac::CReferenceFile newReference;
        newReference.m_specieName = standardCrossSections.ReferenceSpecieName(ii);
        newReference.m_path = dstFileName;
        newReference.m_isFiltered = referenceController.m_highPassFilter;
        referencesCreated.push_back(newReference);
    }

    // Save the Fraunhofer reference as well
    if (IsExistingFile(standardCrossSections.FraunhoferReferenceFileName()))
    {
        // Do the convolution
        referenceController.m_convertToAir = false;
        referenceController.m_highResolutionCrossSection = standardCrossSections.FraunhoferReferenceFileName();
        referenceController.m_isPseudoAbsorber = true;
        referenceController.ConvolveReference(*calibration);

        // Save the result
        const std::string dstFileName =
            directoryName +
            spectrumInformation.m_device +
            "_Fraunhofer_" +
            FormatDateAndTimeOfSpectrum(spectrumInformation) +
            ".txt";

        novac::SaveCrossSectionFile(dstFileName, *(referenceController.m_resultingCrossSection));

        novac::CReferenceFile newReference;
        newReference.m_specieName = "Fraunhofer";
        newReference.m_path = dstFileName;
        newReference.m_isFiltered = referenceController.m_highPassFilter;
        referencesCreated.push_back(newReference);
    }

    return referencesCreated;
}

// ----------- PostCalibration class -----------

bool ScanIsMeasuredInConfiguredTimeOfDayForCalibration(const novac::CDateTime& scanStartTime, const Configuration::CUserConfiguration& userSettings)
{
    // Check if the scan lies within the configured interval. Slightly complex logic here since the interval may wrap around midnight UTC.
    // Notice that this logic is shared with the real-time calibration in the NovacProgram.
    const bool calibrationIntervalWrapsMidnight = userSettings.m_calibrationIntervalTimeOfDayLow > userSettings.m_calibrationIntervalTimeOfDayHigh;
    bool scanTimeLiesWithinCalibrationTimeInterval = false;
    if (calibrationIntervalWrapsMidnight)
    {
        scanTimeLiesWithinCalibrationTimeInterval =
            scanStartTime.SecondsSinceMidnight() >= userSettings.m_calibrationIntervalTimeOfDayLow ||
            scanStartTime.SecondsSinceMidnight() <= userSettings.m_calibrationIntervalTimeOfDayHigh;
    }
    else
    {
        scanTimeLiesWithinCalibrationTimeInterval =
            scanStartTime.SecondsSinceMidnight() >= userSettings.m_calibrationIntervalTimeOfDayLow &&
            scanStartTime.SecondsSinceMidnight() <= userSettings.m_calibrationIntervalTimeOfDayHigh;
    }

    if (!scanTimeLiesWithinCalibrationTimeInterval)
    {
        std::stringstream message;
        message << "Measurement time (" << scanStartTime.SecondsSinceMidnight() << ") is outside of configured interval [";
        message << userSettings.m_calibrationIntervalTimeOfDayLow << " to " << userSettings.m_calibrationIntervalTimeOfDayHigh << "]";
        ShowMessage(message.str());
        return false;
    }

    // no further objections found.
    return true;
}

std::map<SpectrometerId, std::vector<CPostCalibration::BasicScanInfo>> CPostCalibration::SortScanFilesByInstrument(novac::ILogger& log, const std::vector<std::string>& scanFileList)
{
    std::map<SpectrometerId, std::vector<CPostCalibration::BasicScanInfo>> result;

    for (const auto& scanFile : scanFileList)
    {
        novac::CDateTime startTime;
        CString serial;
        int channel;
        MeasurementMode mode;

        BasicScanInfo info;

        CString fileName = scanFile.c_str();
        if (!CFileUtils::GetInfoFromFileName(fileName, startTime, serial, channel, mode))
        {
            novac::LogContext context(novac::LogContext::FileName, scanFile);
            CScanFileHandler scan(log);
            if (!scan.CheckScanFile(context, scanFile))
            {
                log.Error(context, "Could not read pak file");
                continue;
            }

            CSpectrum skySpec;
            if (scan.GetSky(skySpec))
            {
                log.Error(context, "Could not read the sky spectrum from pak file");
                continue;
            }

            info.fullPath = scanFile;
            info.serial = skySpec.m_info.m_device;
            info.channel = skySpec.m_info.m_channel;
            info.startTime = skySpec.m_info.m_startTime;
        }
        else
        {
            info.fullPath = scanFile;
            info.serial = serial;
            info.channel = channel;
            info.startTime = startTime;
        }

        SpectrometerId id(info.serial, info.channel);

        // Insert the result
        auto pos = result.find(id);
        if (pos == result.end())
        {
            std::vector<CPostCalibration::BasicScanInfo> newCollection;
            newCollection.push_back(info);
            result[id] = newCollection;
        }
        else
        {
            pos->second.push_back(info);
        }
    }

    return result;
}

int CPostCalibration::RunInstrumentCalibration(const std::vector<std::string>& scanFileList, CPostCalibrationStatistics& statistics)
{
    auto sortedScanFileList = SortScanFilesByInstrument(m_log, scanFileList);
    {
        std::stringstream message;
        message << "Located pak files from " << sortedScanFileList.size() << " devices";
        ShowMessage(message.str());
    }

    int numberOfCalibrations = 0;

    for (auto& scanFileInfo : sortedScanFileList)
    {
        try
        {
            {
                std::stringstream message;
                message << "Performing calibrations for " << scanFileInfo.first.serial << " (channel: " << scanFileInfo.first.channel << ")";
                ShowMessage(message.str());
            }

            // Sort the files by increasing start time.
            std::sort(
                begin(scanFileInfo.second),
                end(scanFileInfo.second),
                [](const BasicScanInfo& first, const BasicScanInfo& second)
                {
                    return first.startTime < second.startTime;
                });

            novac::CDateTime timeOfLastCalibration;
            for (const auto& basicFileInfo : scanFileInfo.second)
            {
                {
                    std::stringstream message;
                    message << "Checking pak file: " << basicFileInfo.fullPath;
                    ShowMessage(message.str());
                }

                if (!ScanIsMeasuredInConfiguredTimeOfDayForCalibration(basicFileInfo.startTime, m_userSettings))
                {
                    continue;
                }

                const double secondsSinceLastCalibration = timeOfLastCalibration.year > 0 ? novac::CDateTime::Difference(basicFileInfo.startTime, timeOfLastCalibration) : 1e99;
                if (secondsSinceLastCalibration < 3600.0 * m_userSettings.m_calibrationIntervalHours)
                {
                    std::stringstream message;
                    message << "Interval since last performed calibration ( " << (secondsSinceLastCalibration / 3600.0) << " hours) is too small. Skipping scan.";
                    ShowMessage(message.str());
                    continue;
                }

                if (RunInstrumentCalibration(basicFileInfo.fullPath, statistics))
                {
                    ++numberOfCalibrations;
                    timeOfLastCalibration = basicFileInfo.startTime;
                }
            }

            // All calibrations for this particular spectrometer are now done.
            CreateEvaluationSettings(scanFileInfo.first, statistics);
        }
        catch (std::exception& e)
        {
            std::stringstream message;
            message << "Failed to create evaluation data for " << scanFileInfo.first.serial << " (channel: " << scanFileInfo.first.channel << "). ";
            message << "Exception: " << e.what();
            ShowMessage(message.str());
        }
    }

    return numberOfCalibrations;
}

bool CPostCalibration::RunInstrumentCalibration(const std::string& scanFile, CPostCalibrationStatistics& statistics)
{
    try
    {
        novac::ConsoleLog log;
        novac::LogContext context;
        CScanFileHandler scan(log);
        if (!scan.CheckScanFile(context, scanFile))
        {
            std::stringstream message;
            message << "Could not read recieved pak-file '" << scanFile << "' . Will not perform calibration.";
            ShowMessage(message.str());
            return false;
        }

        // Get the sky-spectrum. Read out serial-number and start-time from this
        CSpectrum skySpec;
        if (scan.GetSky(skySpec))
        {
            std::stringstream message;
            message << "Could not read a sky spectrum from pak-file '" << scanFile << "' . Will not perform calibration.";
            ShowMessage(message.str());
            return false;
        }

        const auto* instrument = m_setup.GetInstrument(skySpec.m_info.m_device);
        if (instrument == nullptr)
        {
            std::stringstream message;
            message << "Received pak file from not configured instrument: '" << skySpec.m_info.m_device << "' . Will not perform calibration.";
            ShowMessage(message.str());
            return false;
        }

        // Use the WavelengthCalibrationController, which is also used when the 
        //  user performs the instrument calibrations using the CCalibratePixelToWavelengthDialog.
        // This makes sure we get the same behavior in the dialog and here.
        NovacProgramWavelengthCalibrationController calibrationController(m_log);
        RunCalibration(calibrationController, scanFile, m_userSettings, instrument->m_instrumentCalibration);

        // Save new instrument calibration.
        std::string directoryName = CreateOutputDirectoryForCalibration(skySpec.m_info, m_userSettings);
        const std::string calibrationFileName = directoryName + GetCalibrationFileName(calibrationController.m_calibrationDebug.spectrumInfo);
        calibrationController.SaveResultAsStd(calibrationFileName);

        // Create the standard references.
        const SpectrometerId device{ skySpec.m_info.m_device, skySpec.m_info.m_channel };
        const auto finalCalibration = calibrationController.GetFinalCalibration();
        auto referencesCreated = CreateStandardReferences(
            calibrationController.m_calibrationDebug.spectrumInfo,
            finalCalibration,
            m_standardCrossSections,
            directoryName);

        statistics.RememberCalibrationPerformed(
            device,
            skySpec.m_info.m_startTime,
            referencesCreated);

        {
            std::stringstream message;
            message << "Peformed calibration on '" << scanFile << "'. Created " << referencesCreated.size() << " references.";
            ShowMessage(message.str());
        }

        return true;
    }
    catch (std::exception& e)
    {
        ShowError(e.what());
    }
    return false;
}

void CPostCalibration::CreateEvaluationSettings(const SpectrometerId& spectrometer, const CPostCalibrationStatistics& statistics)
{
    /* 1. Get the existing fit-windows (with time) for this spectrometer
    *  2. Create new fit-windows with the new references and a splitting up the validity time
    *  3. Save the new .exml file.
    */
    const Configuration::CInstrumentConfiguration* instrument = m_setup.GetInstrument(spectrometer.serial);
    if (instrument == nullptr)
    {
        std::stringstream message;
        message << "Failed to retrieve instrument settings for: " << spectrometer.serial << " no new evaluation settings could be generated";
        ShowMessage(message.str());
        return;
    }

    // Get the windows defined for this instrument
    std::vector<Configuration::FitWindowWithTime> originalWindowsForThisChannel;
    std::vector<Configuration::FitWindowWithTime> originalWindowsForOtherChannels;
    for (int idx = 0; idx < instrument->m_eval.NumberOfFitWindows(); ++idx)
    {
        Configuration::FitWindowWithTime window;
        if (!instrument->m_eval.GetFitWindow(idx, window.window, window.validFrom, window.validTo))
        {
            if (window.window.channel == spectrometer.channel)
            {
                originalWindowsForThisChannel.push_back(window);
            }
            else
            {
                originalWindowsForOtherChannels.push_back(window);
            }
        }
    }

    if (originalWindowsForThisChannel.size() == 0)
    {
        std::stringstream message;
        message << "Failed to retrieve any fit window defined for: " << spectrometer.serial << " and channel ";
        message << spectrometer.channel << " no new evaluation settings could be generated";
        ShowMessage(message.str());
        return;
    }

    // If there are still multiple windows, then warn the user that we can only (so far) process the first one.
    if (originalWindowsForThisChannel.size() > 1)
    {
        std::stringstream message;
        message << "Found " << originalWindowsForThisChannel.size() << " fit windows defined for: " << spectrometer.serial << " and channel " << spectrometer.channel;
        message << ". Will generate new evaluation with settings from the first, different configurations is not supported";
        ShowMessage(message.str());
    }

    // Get the newly defined statistics
    const int calibrationNum = statistics.GetNumberOfCalibrationsPerformedFor(spectrometer);

    // Use the new references and the window given above to create a new (timed) fit window.
    std::vector<Configuration::FitWindowWithTime> result;
    for (int idx = 0; idx < calibrationNum; ++idx)
    {
        Configuration::FitWindowWithTime evaluationWindow;
        std::vector<novac::CReferenceFile> references;
        statistics.GetCalibration(spectrometer, idx, evaluationWindow.validFrom, evaluationWindow.validTo, references);

        evaluationWindow.window = originalWindowsForThisChannel.front().window;
        evaluationWindow.window.nRef = 0;
        for (const auto& reference : references)
        {
            if (Equals(reference.m_specieName, "Fraunhofer") && evaluationWindow.window.fraunhoferRef.m_path.size() != 0)
            {
                evaluationWindow.window.fraunhoferRef = reference;
            }
            else
            {
                evaluationWindow.window.ref[evaluationWindow.window.nRef] = reference;
                evaluationWindow.window.ref[evaluationWindow.window.nRef].m_shiftOption = novac::SHIFT_TYPE::SHIFT_FIX;
                evaluationWindow.window.ref[evaluationWindow.window.nRef].m_shiftValue = 0.0;
                evaluationWindow.window.ref[evaluationWindow.window.nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
                evaluationWindow.window.ref[evaluationWindow.window.nRef].m_squeezeValue = 1.0;

                ++evaluationWindow.window.nRef;
            }
        }

        result.push_back(evaluationWindow);
    }

    // Write the result to file (together with the other previous settings).
    std::string directoryName{ (const char*)m_userSettings.m_outputDirectory };
    if (directoryName.back() != '/' && directoryName.back() != '\\')
    {
        directoryName += '/';
    }
    directoryName += "calibration/";
    if (Filesystem::CreateDirectoryStructure(directoryName))
    {
        std::stringstream message;
        message << "Could not create directory for saving evaluation configuration: " << directoryName;
        ShowMessage(message.str());

        // TODO: Another type of exception!
        throw std::invalid_argument(message.str());
    }
    std::string fileName = directoryName + spectrometer.serial + ".exml";

    Configuration::CEvaluationConfiguration newEvaluationSettings;
    newEvaluationSettings.m_serial = spectrometer.serial;

    // Insert the updated windows, as well as the other, not updated, fit windows.
    for (const auto& evaluationWindow : result)
    {
        newEvaluationSettings.InsertFitWindow(evaluationWindow.window, evaluationWindow.validFrom, evaluationWindow.validTo);
    }
    for (size_t otherWindowIdx = 1; otherWindowIdx < originalWindowsForThisChannel.size(); ++otherWindowIdx)
    {
        newEvaluationSettings.InsertFitWindow(originalWindowsForThisChannel[otherWindowIdx].window, originalWindowsForThisChannel[otherWindowIdx].validFrom, originalWindowsForThisChannel[otherWindowIdx].validTo);
    }
    for (const auto& evaluationWindow : originalWindowsForOtherChannels)
    {
        newEvaluationSettings.InsertFitWindow(evaluationWindow.window, evaluationWindow.validFrom, evaluationWindow.validTo);
    }

    FileHandler::CEvaluationConfigurationParser writer{ m_log };
    writer.WriteConfigurationFile(fileName, newEvaluationSettings, instrument->m_darkCurrentCorrection, instrument->m_instrumentCalibration);

}