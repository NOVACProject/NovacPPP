#include "PostCalibration.h"
#include "../Common/Common.h"
#include <SpectralEvaluation/DialogControllers/NovacProgramWavelengthCalibrationController.h>
#include <SpectralEvaluation/DialogControllers/ReferenceCreationController.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Calibration/StandardCrossSectionSetup.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include "../Configuration/NovacPPPConfiguration.h"
#include "../Configuration/UserConfiguration.h"
#include "PostCalibrationStatistics.h"
#include <sstream>
#include <algorithm>
#include <PPPLib/CFileUtils.h>

extern Configuration::CNovacPPPConfiguration g_setup; // <-- The setup of the instruments
extern Configuration::CUserConfiguration g_userSettings; // <-- The settings of the user for what to process

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

std::string CreateOutputDirectoryForCalibration(const CSpectrumInfo& calibratedSpectrum)
{
    CString dateStr;
    dateStr.Format(
        "%04d.%02d.%02d",
        calibratedSpectrum.m_startTime.year,
        calibratedSpectrum.m_startTime.month,
        calibratedSpectrum.m_startTime.day);

    std::string directoryName{ (const char*)g_userSettings.m_outputDirectory };
    directoryName += "/" + dateStr + "/" + calibratedSpectrum.m_device + "/";

    // 4b. Make sure that the folder exists
    int ret = CreateDirectoryStructure(directoryName);
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

void RunCalibration(NovacProgramWavelengthCalibrationController& calibrationController, const std::string& scanFile, const Configuration::CInstrumentCalibrationConfiguration& calibrationSettings)
{
    calibrationController.m_inputSpectrumFile = scanFile;
    calibrationController.m_solarSpectrumFile = g_userSettings.m_highResolutionSolarSpectrumFile;
    calibrationController.m_initialCalibrationFile = calibrationSettings.m_initialCalibrationFile;
    calibrationController.m_initialLineShapeFile = calibrationSettings.m_instrumentLineshapeFile;
    calibrationController.m_instrumentLineShapeFitOption = (WavelengthCalibrationController::InstrumentLineShapeFitOption)g_userSettings.m_calibrationInstrumentLineShapeFitOption;
    calibrationController.m_instrumentLineShapeFitRegion = g_userSettings.m_calibrationInstrumentLineShapeFitRegion;

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

    referenceController.m_highPassFilter = false; // TODO: Are the references sometimes filtered here?
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
    {
        // Do the convolution
        referenceController.m_highPassFilter = false;
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

bool ScanIsMeasuredInConfiguredTimeOfDayForCalibration(const novac::CDateTime& scanStartTime)
{
    // Check if the scan lies within the configured interval. Slightly complex logic here since the interval may wrap around midnight UTC.
    // Notice that this logic is shared with the real-time calibration in the NovacProgram.
    const bool calibrationIntervalWrapsMidnight = g_userSettings.m_calibrationIntervalTimeOfDayLow > g_userSettings.m_calibrationIntervalTimeOfDayHigh;
    bool scanTimeLiesWithinCalibrationTimeInterval = false;
    if (calibrationIntervalWrapsMidnight)
    {
        scanTimeLiesWithinCalibrationTimeInterval =
            scanStartTime.SecondsSinceMidnight() >= g_userSettings.m_calibrationIntervalTimeOfDayLow ||
            scanStartTime.SecondsSinceMidnight() <= g_userSettings.m_calibrationIntervalTimeOfDayHigh;
    }
    else
    {
        scanTimeLiesWithinCalibrationTimeInterval =
            scanStartTime.SecondsSinceMidnight() >= g_userSettings.m_calibrationIntervalTimeOfDayLow &&
            scanStartTime.SecondsSinceMidnight() <= g_userSettings.m_calibrationIntervalTimeOfDayHigh;
    }

    if (!scanTimeLiesWithinCalibrationTimeInterval)
    {
        std::stringstream message;
        message << "Measurement time (" << scanStartTime.SecondsSinceMidnight() << ") is outside of configured interval [";
        message << g_userSettings.m_calibrationIntervalTimeOfDayLow << " to " << g_userSettings.m_calibrationIntervalTimeOfDayHigh << "]";
        ShowMessage(message.str());
        return false;
    }

    // no further objections found.
    return true;
}

std::map<CPostCalibration::SpectrometerId, std::vector<CPostCalibration::BasicScanInfo>> CPostCalibration::SortScanFilesByInstrument(const std::vector<std::string>& scanFileList)
{
    std::map<CPostCalibration::SpectrometerId, std::vector<CPostCalibration::BasicScanInfo>> result;

    for (const auto& scanFile : scanFileList)
    {
        novac::CDateTime startTime;
        CString serial;
        int channel;
        MEASUREMENT_MODE mode;

        BasicScanInfo info;

        CString fileName = scanFile.c_str();
        if (!CFileUtils::GetInfoFromFileName(fileName, startTime, serial, channel, mode))
        {
            CScanFileHandler scan;
            if (SUCCESS != scan.CheckScanFile(scanFile))
            {
                std::stringstream message;
                message << "Could not read pak-file '" << scanFile << "'";
                ShowMessage(message.str());
                continue;
            }

            CSpectrum skySpec;
            if (scan.GetSky(skySpec))
            {
                std::stringstream message;
                message << "Could not read a sky spectrum from pak-file '" << scanFile << "'";
                ShowMessage(message.str());
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
    auto sortedScanFileList = SortScanFilesByInstrument(scanFileList);
    {
        std::stringstream message;
        message << "Located pak files from " << sortedScanFileList.size() << " devices";
        ShowMessage(message.str());
    }

    int numberOfCalibrations = 0;

    for (auto& scanFileInfo : sortedScanFileList)
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
            [](const BasicScanInfo& first, const BasicScanInfo& second) {
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

            if (!ScanIsMeasuredInConfiguredTimeOfDayForCalibration(basicFileInfo.startTime))
            {
                continue;
            }

            const double secondsSinceLastCalibration = timeOfLastCalibration.year > 0 ? novac::CDateTime::Difference(basicFileInfo.startTime, timeOfLastCalibration) : 1e99;
            if (secondsSinceLastCalibration < 3600.0 * g_userSettings.m_calibrationIntervalHours)
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
    }

    return numberOfCalibrations;
}

bool CPostCalibration::RunInstrumentCalibration(const std::string& scanFile, CPostCalibrationStatistics& statistics)
{
    try
    {
        CScanFileHandler scan;
        if (SUCCESS != scan.CheckScanFile(scanFile))
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

        const auto* instrument = g_setup.GetInstrument(skySpec.m_info.m_device);
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
        NovacProgramWavelengthCalibrationController calibrationController;
        RunCalibration(calibrationController, scanFile, instrument->m_instrumentCalibration);

        // Save new instrument calibration.
        std::string directoryName = CreateOutputDirectoryForCalibration(skySpec.m_info);
        const std::string calibrationFileName = directoryName + GetCalibrationFileName(calibrationController.m_calibrationDebug.spectrumInfo);
        calibrationController.SaveResultAsStd(calibrationFileName);

        // Create the standard references.
        const auto finalCalibration = calibrationController.GetFinalCalibration();
        auto referencesCreated = CreateStandardReferences(
            calibrationController.m_calibrationDebug.spectrumInfo,
            finalCalibration,
            m_standardCrossSections,
            directoryName);

        statistics.RememberCalibrationPerformed(
            skySpec.m_info.m_device,
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