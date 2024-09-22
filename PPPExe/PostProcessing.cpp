#include "stdafx.h"
#include "PostProcessing.h"
#include <SpectralEvaluation/File/File.h>

#undef min
#undef max

#include <algorithm>
#include <sstream>

// the PostEvaluationController takes care of the DOAS evaluations
#include "Evaluation/PostEvaluationController.h"

// The PostCalibration takes care of the instrument calibrations.
#include <PPPLib/Calibration/PostCalibration.h>
#include <PPPLib/Calibration/PostCalibrationStatistics.h>

// The FluxCalculator takes care of calculating the fluxes
#include "Flux/FluxCalculator.h"

// The Stratospherecalculator takes care of calculating Stratospheric VCD's
#include "Stratosphere/StratosphereCalculator.h"

// The flux CFluxStatistics takes care of the statistical part of the fluxes
#include "Flux/FluxStatistics.h"

// We also need to read the evaluation-log files
#include "Common/EvaluationLogFileHandler.h"

#include "WindMeasurement/WindSpeedCalculator.h"

#include <PPPLib/Meteorology/XMLWindFileReader.h>
#include <PPPLib/File/Filesystem.h>
#include "Common/EvaluationLogFileHandler.h"

#include <PPPLib/VolcanoInfo.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <PPPLib/ThreadUtils.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Calibration/StandardCrossSectionSetup.h>

// we need to be able to download data from the FTP-server
#include <PPPLib/Communication/FTPServerConnection.h>

#include <Poco/DirectoryIterator.h>
#include <Poco/Exception.h>

#undef min
#undef max

extern novac::CVolcanoInfo  g_volcanoes;   // <-- A list of all known volcanoes

using namespace novac;


// this is the working-thread that takes care of evaluating a portion of the scans
void EvaluateScansThread(const Configuration::CNovacPPPConfiguration& setup, const Configuration::CUserConfiguration& userSettings, CPostProcessingStatistics& processingStats);

// this takes care of adding the evaluated log-files to the list in an synchronized way
//  the parameter passed in a reference to an array of strings holding the names of the 
//  eval-log files generated
void AddResultToList(const novac::CString& pakFileName, const novac::CString(&evalLog)[MAX_FIT_WINDOWS], const CPlumeInScanProperty& scanProperties, const Configuration::CUserConfiguration& userSettings, CPostProcessingStatistics& processingStats);

CPostProcessing::CPostProcessing(ILogger& logger, Configuration::CNovacPPPConfiguration setup, Configuration::CUserConfiguration userSettings)
    : m_log(logger), m_setup(setup), m_userSettings(userSettings)
{
}

void CPostProcessing::DoPostProcessing_Flux()
{
    novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&> evalLogFiles;
    novac::CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*> geometryResults;
    novac::CString messageToUser, windFileName;

    ShowMessage("--- Prepairing to perform Flux Calculations --- ");
    // --------------- PREPARING FOR THE PROCESSING -----------

    // Checks that the evaluation is ok and that all settings makes sense.
    CheckProcessingSettings();

    // Prepare for the flux-calculation by reading in the wind-field
    ReadWindField();

    // Prepare for the flux-calculation by compiling a set of pausible plume heights
    if (PreparePlumeHeights())
    {
        ShowMessage("Exiting post processing");
        return;
    }

    // --------------- DOING THE EVALUATIONS -----------

    if (m_userSettings.m_doEvaluations)
    {
        // Prepare for the evaluation by reading in the reference files
        PrepareEvaluation();

        // 1. Find all .pak files in the directory.
        ShowMessage("--- Locating Pak Files --- ");
        std::vector<std::string> pakFileList;
        if (m_userSettings.m_LocalDirectory.GetLength() > 3)
        {
            messageToUser.Format("Searching for .pak - files in directory %s", (const char*)m_userSettings.m_LocalDirectory);
            ShowMessage(messageToUser);

            const bool includeSubDirs = (m_userSettings.m_includeSubDirectories_Local > 0);
            Filesystem::FileSearchCriterion limits;
            limits.startTime = m_userSettings.m_fromDate;
            limits.endTime = m_userSettings.m_toDate;
            limits.fileExtension = ".pak";
            Filesystem::SearchDirectoryForFiles(m_userSettings.m_LocalDirectory, includeSubDirs, pakFileList, &limits);
        }

        if (m_userSettings.m_FTPDirectory.GetLength() > 9)
        {
            CheckForSpectraOnFTPServer(pakFileList);
        }
        if (pakFileList.size() == 0)
        {
            ShowMessage("No spectrum files found. Exiting");
            return;
        }

        // Evaluate the scans. This at the same time generates a list of evaluation-log
        // files with the evaluated results
        ShowMessage("--- Running Evaluations --- ");
        EvaluateScans(pakFileList, evalLogFiles);
        messageToUser.Format("%d evaluation log files accepted", evalLogFiles.GetCount());
        ShowMessage(messageToUser);
    }
    else
    {
        // Don't evaluate the scans, just read the log-files and calculate fluxes from there.
        LocateEvaluationLogFiles(m_userSettings.m_outputDirectory, evalLogFiles);
    }

    // Sort the evaluation-logs in order of increasing start-time, this to make
    // the looking for matching files in 'CalculateGeometries' faster
    ShowMessage("Evaluation done. Sorting the evaluation log files");
    SortEvaluationLogs(evalLogFiles);
    ShowMessage("Sort done.");

    // 3. Loop through list with output text files from evaluation and calculate
    //      the geometries
    CalculateGeometries(evalLogFiles, geometryResults);

    // 4.1 write the calculations to file, for later checking or other uses...
    WriteCalculatedGeometriesToFile(geometryResults);

    // 4.2 Insert the calculated geometries into the plume height database
    InsertCalculatedGeometriesIntoDataBase(geometryResults);

    // 5. Calculate the wind-speeds from the wind-speed measurements
    //  the plume heights are taken from the database
    CalculateDualBeamWindSpeeds(evalLogFiles);

    // 6. Calculate flux from evaluation text files
    CalculateFluxes(evalLogFiles);

    // 7. Write the statistics
    novac::CString statFileName;
    statFileName.Format("%s%cProcessingStatistics.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    Common::ArchiveFile(statFileName);
    m_processingStats.WriteStatToFile(statFileName);

    // 8. Also write the wind field that we have created to file
    windFileName.Format("%s%cGeneratedWindField.wxml", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    Common::ArchiveFile(windFileName);
    m_windDataBase.WriteToFile(windFileName);

    // 9. Upload the results to the FTP-server
    if (m_userSettings.m_uploadResults)
    {
        UploadResultsToFTP();
    }

    // ------------ Clean up -----------
    while (geometryResults.GetSize() != 0)
    {
        auto p = geometryResults.GetTailPosition();
        Geometry::CGeometryResult* g = geometryResults.GetAt(p);
        delete g;
        geometryResults.RemoveTail();
    }
}

void CPostProcessing::DoPostProcessing_InstrumentCalibration()
{
    ShowMessage("--- Prepairing to perform Instrument Calibrations --- ");

    ShowMessage("--- Validating settings --- ");
    if (CheckInstrumentCalibrationSettings())
    {
        ShowMessage("Exiting post processing");
        return;
    }

    novac::StandardCrossSectionSetup standardCrossSections{ m_exePath };
    if (standardCrossSections.NumberOfReferences() == 0)
    {
        ShowMessage("The StandardReferences folder is either missing or contains no references. No references can be created.");
        return;
    }
    else if (!IsExistingFile(standardCrossSections.FraunhoferReferenceFileName()))
    {
        ShowMessage("Cannot locate the Fraunhofer reference in the StandardReferences folder. The file is either missing or path is invalid. No Fraunhofer reference can be created.");
    }


    // 1. Find all .pak files in the directory.
    ShowMessage("--- Locating Pak Files --- ");
    std::vector<std::string> pakFileList;
    if (m_userSettings.m_LocalDirectory.GetLength() > 3)
    {
        CString messageToUser;
        messageToUser.Format("Searching for .pak - files in directory %s", (const char*)m_userSettings.m_LocalDirectory);
        ShowMessage(messageToUser);

        const bool includeSubDirs = (m_userSettings.m_includeSubDirectories_Local > 0);
        Filesystem::FileSearchCriterion limits;
        limits.startTime = m_userSettings.m_fromDate;
        limits.endTime = m_userSettings.m_toDate;
        limits.fileExtension = ".pak";
        Filesystem::SearchDirectoryForFiles(m_userSettings.m_LocalDirectory, includeSubDirs, pakFileList, &limits);
    }

    if (m_userSettings.m_FTPDirectory.GetLength() > 9)
    {
        CheckForSpectraOnFTPServer(pakFileList);
    }
    if (pakFileList.size() == 0)
    {
        ShowMessage("No spectrum files found. Exiting");
        return;
    }

    ShowMessage("--- Running Calibrations --- ");

    novac::CPostCalibration calibrationController{ standardCrossSections, m_log };
    novac::CPostCalibrationStatistics calibrationStatistics;

    // Unlike other parts of the NovacPPP, this function is intended to be single threaded. The reason for this is that the
    //  classes which are called upon are strongly threaded themselves and further threading will not help the performance.
    int numberOfCalibrations = calibrationController.RunInstrumentCalibration(pakFileList, calibrationStatistics);

    {
        CString messageToUser;
        messageToUser.Format("%d instrument calibrations performed", numberOfCalibrations);
        ShowMessage(messageToUser);
    }
}

void CPostProcessing::DoPostProcessing_Strat()
{
    novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&> evalLogFiles;
    novac::CString messageToUser, statFileName;
    Stratosphere::CStratosphereCalculator strat;

    // --------------- PREPARING FOR THE PROCESSING -----------

    // Checks that the evaluation is ok and that all settings makes sense.
    CheckProcessingSettings();

    // Prepare for the evaluation by reading in the reference files
    PrepareEvaluation();

    // --------------- DOING THE PROCESSING -----------

    // 1. Find all .pak files in the directory.
    std::vector<std::string> pakFileList;
    if (m_userSettings.m_LocalDirectory.GetLength() > 3)
    {
        messageToUser.Format("Searching for .pak - files in directory %s", (const char*)m_userSettings.m_LocalDirectory);
        ShowMessage(messageToUser);

        const bool includeSubDirs = (m_userSettings.m_includeSubDirectories_Local > 0);
        Filesystem::FileSearchCriterion limits;
        limits.startTime = m_userSettings.m_fromDate;
        limits.endTime = m_userSettings.m_toDate;
        limits.fileExtension = ".pak";
        Filesystem::SearchDirectoryForFiles(m_userSettings.m_LocalDirectory, includeSubDirs, pakFileList, &limits);
    }
    if (m_userSettings.m_FTPDirectory.GetLength() > 9)
    {
        CheckForSpectraOnFTPServer(pakFileList);
    }
    if (pakFileList.size() == 0)
    {
        ShowMessage("No spectrum files found. Exiting");
        return;
    }

    // Evaluate the scans. This at the same time generates a list of evaluation-log
    // files with the evaluated results
    EvaluateScans(pakFileList, evalLogFiles);
    messageToUser.Format("%d evaluation log files accepted", evalLogFiles.GetCount());
    ShowMessage(messageToUser);

    // Sort the evaluation-logs in order of increasing start-time, this to make
    // the looking for matching files in 'CalculateGeometries' faster
    ShowMessage("Evaluation done. Sorting the evaluation log files");
    SortEvaluationLogs(evalLogFiles);
    ShowMessage("Sort done.");

    ShowMessage("Calculating stratospheric columns");
    // strat.CalculateVCDs(evalLogFiles);
    ShowMessage("Calculation Done");

    // 7. Write the statistics
    statFileName.Format("%s%cProcessingStatistics.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    Common::ArchiveFile(statFileName);
    m_processingStats.WriteStatToFile(statFileName);
}

void CPostProcessing::CheckForSpectraOnFTPServer(std::vector<std::string>& fileList)
{
    Communication::CFTPServerConnection serverDownload;

    int ret = serverDownload.DownloadDataFromFTP(
        m_userSettings.m_FTPDirectory,
        m_userSettings.m_FTPUsername,
        m_userSettings.m_FTPPassword,
        fileList);

    if (ret == 0)
    {
        ShowMessage("Successfully downloaded all data files.");
    }
    else
    {
        ShowMessage("Error happened when downloading data from FTP.");
    }
}

novac::GuardedList<std::string> s_pakFilesRemaining;
novac::GuardedList<Evaluation::CExtendedScanResult> s_evalLogs;

volatile unsigned long s_nFilesToProcess;

void CPostProcessing::EvaluateScans(
    const std::vector<std::string>& pakFileList,
    novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&>& evalLogFiles) const
{
    s_nFilesToProcess = (long)pakFileList.size();
    novac::CString messageToUser;

    // share the list of pak-files with the other functions around here
    for (const std::string& file : pakFileList)
    {
        s_pakFilesRemaining.AddItem(file);
    }

    // Keep the user informed about what we're doing
    messageToUser.Format("%ld spectrum files found. Begin evaluation using %d threads.", s_nFilesToProcess, m_userSettings.m_maxThreadNum);
    ShowMessage(messageToUser);

    // start the threads
    std::vector<std::thread> evalThreads(m_userSettings.m_maxThreadNum);
    for (unsigned int threadIdx = 0; threadIdx < m_userSettings.m_maxThreadNum; ++threadIdx)
    {
        std::thread t{ EvaluateScansThread, m_setup, m_userSettings, m_processingStats };
        evalThreads[threadIdx] = std::move(t);
    }

    // make sure that all threads have time to finish before we say that we're ready
    for (unsigned int threadIdx = 0; threadIdx < m_userSettings.m_maxThreadNum; ++threadIdx)
    {
        evalThreads[threadIdx].join();
    }

    // copy out the result
    s_evalLogs.CopyTo(evalLogFiles);

    messageToUser.Format("All %ld scans evaluated.", s_nFilesToProcess);
    ShowMessage(messageToUser);
}

void EvaluateScansThread(const Configuration::CNovacPPPConfiguration& setup, const Configuration::CUserConfiguration& userSettings, CPostProcessingStatistics& processingStats)
{
    std::string fileName;

    // create a new CPostEvaluationController
    Evaluation::CPostEvaluationController eval{ setup, userSettings, processingStats };

    // while there are more .pak-files
    while (s_pakFilesRemaining.PopFront(fileName))
    {
        novac::CString evalLog[MAX_FIT_WINDOWS];
        CPlumeInScanProperty scanProperties[MAX_FIT_WINDOWS];

        // evaluate the .pak-file in all the specified fit-windows and retrieve the name of the 
        // eval-logs. If any of the fit-windows fails then the scan is not inserted.
        bool evaluationSucceeded = true;
        for (int fitWindowIndex = 0; fitWindowIndex < userSettings.m_nFitWindowsToUse; ++fitWindowIndex)
        {
            if (0 != eval.EvaluateScan(fileName, userSettings.m_fitWindowsToUse[fitWindowIndex], &evalLog[fitWindowIndex], &scanProperties[fitWindowIndex]))
            {
                evaluationSucceeded = false;
                break;
            }
        }

        if (evaluationSucceeded)
        {
            // If we made it this far then the measurement is ok, insert it into the list!
            AddResultToList(fileName, evalLog, scanProperties[userSettings.m_mainFitWindow], userSettings, processingStats);

            // Tell the user what is happening
            novac::CString messageToUser;
            messageToUser.Format(" + Inserted scan %s into list of evaluation logs", fileName.c_str());
            ShowMessage(messageToUser);
        }
        else
        {
            novac::CString messageToUser;
            messageToUser.Format(" - No flux calculated for scan %s ", fileName.c_str());
            ShowMessage(messageToUser);
        }
    }
}

void AddResultToList(
    const novac::CString& pakFileName,
    const novac::CString(&evalLog)[MAX_FIT_WINDOWS],
    const CPlumeInScanProperty& scanProperties,
    const Configuration::CUserConfiguration& userSettings,
    CPostProcessingStatistics& processingStats)
{
    // these are not used...
    novac::CString serial;
    int channel;
    MEASUREMENT_MODE mode;

    // Create a new Extended scan result and add it to the end of the list
    Evaluation::CExtendedScanResult newResult;
    newResult.m_pakFile.Format(pakFileName);
    for (int fitWindowIndex = 0; fitWindowIndex < userSettings.m_nFitWindowsToUse; ++fitWindowIndex)
    {
        newResult.m_evalLogFile[fitWindowIndex].Format(evalLog[fitWindowIndex]);
        newResult.m_fitWindowName[fitWindowIndex].Format(userSettings.m_fitWindowsToUse[fitWindowIndex]);
    }
    novac::CFileUtils::GetInfoFromFileName(evalLog[0], newResult.m_startTime, serial, channel, mode);
    newResult.m_scanProperties = scanProperties;

    // store the name of the evaluation-log file generated
    s_evalLogs.AddItem(newResult);

    // update the statistics
    processingStats.InsertAcception(serial);
}

int CPostProcessing::CheckInstrumentCalibrationSettings() const
{
    novac::CString errorMessage;
    CDateTime now;
    int returnCode = 0;

    if (!novac::IsExistingFile(m_userSettings.m_highResolutionSolarSpectrumFile))
    {
        std::stringstream message;
        message << "The high resolution solar spectrum file could not be found. File given: '" << m_userSettings.m_highResolutionSolarSpectrumFile << "'" << std::endl;
        ShowError(message.str());
        returnCode = 1;
    }

    if (m_userSettings.m_calibrationIntervalTimeOfDayLow < 0 || m_userSettings.m_calibrationIntervalTimeOfDayLow > 86400 ||
        m_userSettings.m_calibrationIntervalTimeOfDayHigh < 0 || m_userSettings.m_calibrationIntervalTimeOfDayHigh > 86400)
    {
        std::stringstream message;
        message << "The calibration intervals lie outside of the allowed interval [0, 86400]" << std::endl;
        ShowError(message.str());
        returnCode = 1;
    }

    // Verify that each instrument has an initial instrument line shape file.
    for (int instrumentIdx = 0; instrumentIdx < m_setup.NumberOfInstruments(); ++instrumentIdx)
    {
        if (m_setup.m_instrument[instrumentIdx].m_instrumentCalibration.m_initialCalibrationFile.size() == 0)
        {
            std::stringstream message;
            message << "No initial calibration file defined for instrument " << (const char*)m_setup.m_instrument[instrumentIdx].m_serial << std::endl;
            ShowError(message.str());
            returnCode = 1;
        }
        if (!novac::IsExistingFile(m_setup.m_instrument[instrumentIdx].m_instrumentCalibration.m_initialCalibrationFile))
        {
            std::stringstream message;
            message << "Cannot locate the initial calibration file defined for instrument " << (const char*)m_setup.m_instrument[instrumentIdx].m_serial << std::endl;
            ShowError(message.str());
            returnCode = 1;
        }

        // The initial instrument line shape file is optional, but if it is defined then it must also exist.
        if (m_setup.m_instrument[instrumentIdx].m_instrumentCalibration.m_instrumentLineshapeFile.size() > 0 &&
            !novac::IsExistingFile(m_setup.m_instrument[instrumentIdx].m_instrumentCalibration.m_instrumentLineshapeFile))
        {
            std::stringstream message;
            message << "Cannot locate the initial instrument line shape file defined for instrument " << (const char*)m_setup.m_instrument[instrumentIdx].m_serial << std::endl;
            ShowError(message.str());
            returnCode = 1;
        }
    }
    return returnCode;
}

void CPostProcessing::CheckProcessingSettings() const
{
    novac::CString errorMessage;
    CDateTime now;

    ShowMessage("--- Validating settings --- ");

    // Check that no instrument is duplicated in the list of instruments...
    for (int j = 0; j < m_setup.NumberOfInstruments(); ++j)
    {
        for (int k = j + 1; k < m_setup.NumberOfInstruments(); ++k)
        {
            if (Equals(m_setup.m_instrument[j].m_serial, m_setup.m_instrument[k].m_serial))
            {
                errorMessage.Format("The instrument %s is defined twice in setup.xml. If the instrument has two locations then define the instrument once but with two locations. Exiting post processsing.", (const char*)m_setup.m_instrument[k].m_serial);
                throw std::invalid_argument(errorMessage.std_str());
            }
        }
    }

    // Check that, for each spectrometer, there's only one fit-window defined
    // at each instant
    for (int j = 0; j < m_setup.NumberOfInstruments(); ++j)
    {
        if (m_setup.m_instrument[j].m_eval.NumberOfFitWindows() == 1)
        {
            continue;
        }
        m_setup.m_instrument[j].m_eval.CheckSettings();
    }

    // Check that, for each spectrometer, there's only one location defined at each instant
    for (int j = 0; j < m_setup.NumberOfInstruments(); ++j)
    {
        if (m_setup.m_instrument[j].m_location.GetLocationNum() == 1)
        {
            continue;
        }
        m_setup.m_instrument[j].m_location.CheckSettings();
    }
}

void CPostProcessing::PrepareEvaluation()
{
    ShowMessage("--- Reading References --- ");

    CDateTime fromTime, toTime; //  these are not used but must be passed onto GetFitWindow...
    novac::CString errorMessage;

    // this is true if we failed to prepare the evaluation...
    bool failure = false;

    // Loop through each of the configured instruments
    for (int instrumentIndex = 0; instrumentIndex < m_setup.NumberOfInstruments(); ++instrumentIndex)
    {
        // For each instrument, loop through the fit-windows that are defined
        int fitWindowNum = m_setup.m_instrument[instrumentIndex].m_eval.NumberOfFitWindows();
        for (int fitWindowIndex = 0; fitWindowIndex < fitWindowNum; ++fitWindowIndex)
        {
            novac::CFitWindow window;

            // get the fit window
            m_setup.m_instrument[instrumentIndex].m_eval.GetFitWindow(fitWindowIndex, window, fromTime, toTime);

            // For each reference in the fit-window, read it in and make sure that it exists...
            for (int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
            {

                if (window.ref[referenceIndex].m_path.empty())
                {
                    // The reference file was not given in the configuration. Try to generate a configuration
                    //  from the cross section, slit-function and wavelength calibration. These three must then 
                    //  exist or the evaluation fails.
                    if (!ConvolveReference(window.ref[referenceIndex], m_setup.m_instrument[instrumentIndex].m_serial))
                    {
                        errorMessage.Format("Failed to create reference for fit window '%s' for spectrometer %s.",
                            window.ref[referenceIndex].m_specieName.c_str(), (const char*)m_setup.m_instrument[instrumentIndex].m_serial);
                        ShowMessage(errorMessage);
                        failure = true;
                        continue;
                    }
                }
                else
                {
                    if (!IsExistingFile(window.ref[referenceIndex].m_path))
                    {
                        // the file does not exist, try to change it to include the path of the configuration-directory...
                        std::string fileName = Filesystem::GetAbsolutePathFromRelative(window.ref[referenceIndex].m_path, this->m_exePath);

                        if (Filesystem::IsExistingFile(fileName))
                        {
                            window.ref[referenceIndex].m_path = fileName;
                        }
                        else
                        {
                            errorMessage.Format("Cannot find reference file %s", window.ref[referenceIndex].m_path.c_str());
                            ShowMessage(errorMessage);
                            failure = true;
                            continue;
                        }
                    }

                    // Read in the cross section
                    if (window.ref[referenceIndex].ReadCrossSectionDataFromFile())
                    {
                        errorMessage.Format("Failed to read cross section file: %s", window.ref[referenceIndex].m_path.c_str());
                        ShowMessage(errorMessage);
                        failure = true;
                        continue;
                    }

                    // If we are supposed to high-pass filter the spectra then
                    // we should also high-pass filter the cross-sections
                    if (window.fitType == novac::FIT_TYPE::FIT_HP_DIV || window.fitType == novac::FIT_TYPE::FIT_HP_SUB)
                    {
                        if (window.ref[referenceIndex].m_isFiltered == false)
                        {
                            if (novac::Equals(window.ref[referenceIndex].m_specieName, "ring"))
                            {
                                HighPassFilter_Ring(*window.ref[referenceIndex].m_data);
                            }
                            else
                            {
                                HighPassFilter(*window.ref[referenceIndex].m_data);
                            }
                        }
                        else
                        {
                            // The filtered cross sections are scaled to the unit of ppmm
                            // rescale them to molecules/cm2 as all other cross sections
                            Multiply(*window.ref[referenceIndex].m_data, (1.0 / 2.5e15));
                        }
                    }
                }// endif
            }

            // If the window also contains a fraunhofer-reference then read it too.
            if (window.fraunhoferRef.m_path.size() > 4)
            {
                if (!IsExistingFile(window.fraunhoferRef.m_path))
                {

                    // the file does not exist, try to change it to include the path of the configuration-directory...
                    std::string fileName = Filesystem::GetAbsolutePathFromRelative(window.fraunhoferRef.m_path, this->m_exePath);

                    if (Filesystem::IsExistingFile(fileName))
                    {
                        window.fraunhoferRef.m_path = fileName;
                    }
                    else
                    {
                        errorMessage.Format("Cannot find reference file %s", window.fraunhoferRef.m_path.c_str());
                        ShowMessage(errorMessage);
                        failure = true;
                        continue;
                    }
                }

                if (window.fraunhoferRef.ReadCrossSectionDataFromFile())
                {
                    errorMessage.Format("Failed to read fraunhofer-reference file: %s", window.fraunhoferRef.m_path.c_str());
                    ShowMessage(errorMessage);
                    failure = true;
                    continue;
                }
                if (window.fitType == novac::FIT_TYPE::FIT_HP_DIV || window.fitType == novac::FIT_TYPE::FIT_HP_SUB)
                {
                    HighPassFilter_Ring(*window.fraunhoferRef.m_data);
                }
                else
                {
                    Log(*window.fraunhoferRef.m_data);
                }
            }

            // If we've made it this far, then we've managed to read in all the references. Now
            // store the data in m_setup
            m_setup.m_instrument[instrumentIndex].m_eval.SetFitWindow(fitWindowIndex, window, fromTime, toTime);
        }
    }

    if (failure)
    {
        throw std::invalid_argument("failed to setup evaluation");
    }
}

void CPostProcessing::ReadWindField()
{
    ShowMessage("--- Reading Wind information --- ");

    novac::CString name1, name2, name3, path1, path2, path3, messageToUser;
    Common common;
    FileHandler::CXMLWindFileReader reader{ m_log, m_userSettings };

    if (m_userSettings.m_windFieldFileOption == 0)
    {
        // If the user has given a file-name, then try to use that one
        if (m_userSettings.m_windFieldFile.GetLength() > 3)
        {
            messageToUser.Format("Reading wind field from file: %s", (const char*)m_userSettings.m_windFieldFile);
            ShowMessage(messageToUser);

            reader.ReadWindFile(m_userSettings.m_windFieldFile, m_windDataBase);

            messageToUser.Format("Parsed %s containing %d wind data items", (const char*)m_userSettings.m_windFieldFile, m_windDataBase.GetDataBaseSize());
            ShowMessage(messageToUser);

            name1.Format("%sParsedWindField.wxml", (const char*)common.m_exePath);
            m_windDataBase.WriteToFile(name1);

            return;
        }

        // Get the name of the volcano that we are about to process...
        //  there are two options, either the full name or the simple name
        g_volcanoes.GetVolcanoName(m_userSettings.m_volcano, name1);
        g_volcanoes.GetSimpleVolcanoName(m_userSettings.m_volcano, name2);
        g_volcanoes.GetVolcanoCode(m_userSettings.m_volcano, name3);

        // Get the path to the executable, so that we know where to start looking
        path1.Format("%sconfiguration%c%s.wxml", (const char*)common.m_exePath, Poco::Path::separator(), (const char*)name1);
        path2.Format("%sconfiguration%c%s.wxml", (const char*)common.m_exePath, Poco::Path::separator(), (const char*)name2);
        path3.Format("%sconfiguration%c%s.wxml", (const char*)common.m_exePath, Poco::Path::separator(), (const char*)name3);

        // check which of the files exists
        if (Filesystem::IsExistingFile(path1))
        {
            messageToUser.Format("Reading wind field from file: %s", (const char*)path1);
            ShowMessage(messageToUser);

            reader.ReadWindFile(path1, m_windDataBase);

            return;
        }

        if (Filesystem::IsExistingFile(path2))
        {
            messageToUser.Format("Reading wind field from file: %s", (const char*)path2);
            ShowMessage(messageToUser);

            reader.ReadWindFile(path2, m_windDataBase);

            return;
        }

        if (Filesystem::IsExistingFile(path3))
        {
            messageToUser.Format("Reading wind field from file: %s", (const char*)path3);
            ShowMessage(messageToUser);

            reader.ReadWindFile(path3, m_windDataBase);

            return;
        }

        messageToUser.Format("Cannot find wind field file. Attempted to read: %s, %s and %s", (const char*)path1, (const char*)path2, (const char*)path3);
        throw std::invalid_argument(messageToUser.c_str());
    }

    // If the user has specified a directory of files...
    if (m_userSettings.m_windFieldFileOption == 1)
    {
        reader.ReadWindDirectory(m_userSettings.m_windFieldFile, m_windDataBase, &m_userSettings.m_fromDate, &m_userSettings.m_toDate);
        return;
    }

    // should never get to this point!
    throw std::logic_error("Invalid wind field file option, failed to read winf field.");
}

int CPostProcessing::PreparePlumeHeights()
{
    // we need to construct a default plume height to use, if there's nothing else...
    Geometry::CPlumeHeight plumeHeight;
    plumeHeight.m_plumeAltitude = g_volcanoes.GetPeakAltitude(m_userSettings.m_volcano);
    plumeHeight.m_plumeAltitudeSource = Meteorology::MET_DEFAULT;
    plumeHeight.m_validFrom = CDateTime(0, 0, 0, 0, 0, 0);
    plumeHeight.m_validTo = CDateTime(9999, 12, 31, 23, 59, 59);

    // the estimated plume height is half of the altitude difference between the highest
    // instrument for this volcano and the volcano altitude
    double maxInstrumentAltitude = -1e6;
    Configuration::CInstrumentLocation location;
    novac::CString volcanoName;
    g_volcanoes.GetVolcanoName(m_userSettings.m_volcano, volcanoName);
    for (int k = 0; k < m_setup.NumberOfInstruments(); ++k)
    {
        unsigned long N = m_setup.m_instrument[k].m_location.GetLocationNum();
        for (unsigned int j = 0; j < N; ++j)
        {
            m_setup.m_instrument[k].m_location.GetLocation(j, location);
            if (Equals(volcanoName, location.m_volcano))
            {
                maxInstrumentAltitude = std::max(maxInstrumentAltitude, double(location.m_altitude));
            }
        }
    }
    if (maxInstrumentAltitude > 0)
    {
        plumeHeight.m_plumeAltitudeError = fabs(g_volcanoes.GetPeakAltitude(m_userSettings.m_volcano) - maxInstrumentAltitude) / 2.0;
    }
    else
    {
        plumeHeight.m_plumeAltitudeError = g_volcanoes.GetPeakAltitude(m_userSettings.m_volcano) / 2.0;
    }

    m_plumeDataBase.InsertPlumeHeight(plumeHeight);

    return 0;
}

void CPostProcessing::CalculateGeometries(novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&>& evalLogFiles, novac::CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*>& geometryResults)
{
    novac::CString serial1, serial2, messageToUser;
    CDateTime startTime1, startTime2;
    MEASUREMENT_MODE measMode1, measMode2;
    int channel;
    unsigned long nFilesChecked1 = 0; // this is for debugging purposes...
    unsigned long nFilesChecked2 = 0; // this is for debugging purposes...
    unsigned long nCalculationsMade = 0; // this is for debugging purposes...
    unsigned long nTooLongdistance = 0; // this is for debugging purposes...
    unsigned long nTooLargeAbsoluteError = 0; // this is for debugging purposes...
    unsigned long nTooLargeRelativeError = 0; // this is for debugging purposes...
    Configuration::CInstrumentLocation location[2];

    // Tell the user what's happening
    ShowMessage("Begin to calculate plume heights from scans");

    // Loop through list with output text files from evaluation and apply geometrical corrections
    auto pos1 = evalLogFiles.GetHeadPosition();
    while (pos1 != nullptr)
    {
        const novac::CString& evalLog1 = evalLogFiles.GetAt(pos1).m_evalLogFile[m_userSettings.m_mainFitWindow];
        const CPlumeInScanProperty& plume1 = evalLogFiles.GetNext(pos1).m_scanProperties;

        ++nFilesChecked1; // for debugging...

        // if this is the last file in the list, then
        // quit. There's nothing more to compare to...
        if (pos1 == nullptr)
        {
            break;
        }

        // if this scan does not see a large enough portion of the plume, then ignore it...
        if (plume1.completeness < m_userSettings.m_calcGeometry_CompletenessLimit)
        {
            continue;
        }

        //  Get the information about evaluation log file #1
        novac::CFileUtils::GetInfoFromFileName(evalLog1, startTime1, serial1, channel, measMode1);

        // If this is not a flux-measurement, then there's no use in trying to use it...
        if (measMode1 != MODE_FLUX)
        {
            continue;
        }

        // try to combine this evaluation-log file with every other eval-log
        //  use the fact that the list of eval-logs is sorted by increasing start-time
        //  thus we start at the eval-log next after this one and compare with all
        //  eval-logs until the difference in start-time is too big.
        auto pos2 = pos1;
        evalLogFiles.GetNext(pos2);
        bool successfullyCombined = false; // this is true if evalLog1 was combined with (at least one) other eval-log to make a geomery calculation.
        while (pos2 != nullptr)
        {
            const novac::CString& evalLog2 = evalLogFiles.GetAt(pos2).m_evalLogFile[m_userSettings.m_mainFitWindow];
            const CPlumeInScanProperty& plume2 = evalLogFiles.GetNext(pos2).m_scanProperties;

            ++nFilesChecked2; // for debugging...

            // if this scan does not see a large enough portion of the plume, then ignore it...
            if (plume2.completeness < m_userSettings.m_calcGeometry_CompletenessLimit)
            {
                continue;
            }

            //  Get the information about evaluation log file # 2
            novac::CFileUtils::GetInfoFromFileName(evalLog2, startTime2, serial2, channel, measMode2);

            // The time elapsed between the two measurements must not be more than 
            // the user defined time-limit (in seconds)
            double timeDifference = fabs(CDateTime::Difference(startTime1, startTime2));
            if (timeDifference > m_userSettings.m_calcGeometry_MaxTimeDifference)
            {
                pos2 = nullptr;
                continue;
            }

            // If this is not a flux-measurement, then there's no use in trying to use it...
            if (measMode2 != MODE_FLUX)
            {
                continue;
            }

            // the serials must be different (i.e. the two measurements must be
            //  from two different instruments)
            if (Equals(serial1, serial2))
            {
                continue;
            }

            // Get the locations of the two instruments
            try
            {
                location[0] = m_setup.GetInstrumentLocation(serial1.std_str(), startTime1);
                location[1] = m_setup.GetInstrumentLocation(serial2.std_str(), startTime1);
            }
            catch (PPPLib::NotFoundException& ex)
            {
                ShowMessage(ex.message);
                continue;
            }

            // make sure that the distance between the instruments is not too long....
            double instrumentDistance = Common::GPSDistance(location[0].m_latitude, location[0].m_longitude, location[1].m_latitude, location[1].m_longitude);
            if (instrumentDistance < m_userSettings.m_calcGeometry_MinDistance || instrumentDistance > m_userSettings.m_calcGeometry_MaxDistance)
            {
                ++nTooLongdistance;
                continue;
            }

            // count the number of times we calculate a result, for improving the software...
            ++nCalculationsMade;

            // If the files have passed these tests then make a geometry-calculation
            Geometry::CGeometryResult* result = new Geometry::CGeometryResult();
            if (Geometry::CGeometryCalculator::CalculateGeometry(plume1, startTime1, plume2, startTime2, location, *result))
            {
                // Check the quality of the measurement before we insert it...
                if (result->m_plumeAltitudeError > m_userSettings.m_calcGeometry_MaxPlumeAltError)
                {
                    ++nTooLargeAbsoluteError;
                    delete result; // too bad, continue.
                }
                else if ((result->m_plumeAltitudeError > 0.5 * result->m_plumeAltitude) || (result->m_windDirectionError > m_userSettings.m_calcGeometry_MaxWindDirectionError))
                {
                    ++nTooLargeRelativeError;
                    delete result; // too bad, continue.
                }
                else if (result->m_windDirectionError > m_userSettings.m_calcGeometry_MaxWindDirectionError)
                {
                    delete result; // too bad, continue.
                }
                else
                {
                    // remember which instruments were used
                    result->m_instr1.Format(serial1);
                    result->m_instr2.Format(serial2);

                    geometryResults.AddTail(result);

                    messageToUser.Format(" + Calculated a plume altitude of %.0lf +- %.0lf meters and wind direction of %.0lf +- %.0lf degrees by combining measurements from %s and %s",
                        result->m_plumeAltitude, result->m_plumeAltitudeError, result->m_windDirection, result->m_windDirectionError, (const char*)serial1, (const char*)serial2);
                    ShowMessage(messageToUser);

                    successfullyCombined = true;
                }
            }
            else
            {
                // something went wrong... delete the 'info'
                delete result;
            }
        } // end while(pos2 != nullptr)

        // if it was not possible to combine this scan with any other to generate an
        // estimated plume height and wind direction we might still be able to use it to calculate
        // a wind direction given the plume height at the time of the measurement.
        if (!successfullyCombined)
        {
            Meteorology::CWindField windField;
            Geometry::CPlumeHeight plumeHeight;

            // Get the location of the instrument
            try
            {
                location[0] = m_setup.GetInstrumentLocation(serial1.std_str(), startTime1);
            }
            catch (PPPLib::NotFoundException& ex)
            {
                ShowMessage(ex.message);
                continue;
            }

            Geometry::CGeometryResult* result = new Geometry::CGeometryResult();

            // Get the altitude of the plume at this moment. First look into the
            // general database. Then have a look in the list of geometry-results
            // that we just generated to see if there's anything better there...
            m_plumeDataBase.GetPlumeHeight(startTime1, plumeHeight);
            auto gp = geometryResults.GetTailPosition();
            while (gp != nullptr)
            {
                const Geometry::CGeometryResult* oldResult = geometryResults.GetPrev(gp);
                if (fabs(CDateTime::Difference(oldResult->m_averageStartTime, startTime1)) < m_userSettings.m_calcGeometryValidTime)
                {
                    if ((oldResult->m_plumeAltitudeError < plumeHeight.m_plumeAltitudeError) && (oldResult->m_plumeAltitude > NOT_A_NUMBER))
                    {
                        plumeHeight.m_plumeAltitude = oldResult->m_plumeAltitude;
                        plumeHeight.m_plumeAltitudeError = oldResult->m_plumeAltitudeError;
                        plumeHeight.m_plumeAltitudeSource = oldResult->m_calculationType;
                    }
                }
            }

            // Try to calculate the wind-direction
            if (Geometry::CGeometryCalculator::CalculateWindDirection(evalLog1, 0, plumeHeight, location[0], *result))
            {
                // Success!!
                result->m_instr1.Format(serial1);
                geometryResults.AddTail(result);

                // tell the user   
                messageToUser.Format(" + Calculated a wind direction of %.0lf +- %.0lf degrees from a scan by instrument %s",
                    result->m_windDirection, result->m_windDirectionError, (const char*)serial1);
                ShowMessage(messageToUser);
            }
            else
            {
                delete result;
                continue;
            }
        }
    } // end while(pos1 != nullptr)

    // Tell the user what we have done
    if (geometryResults.GetCount() == 0)
    {
        ShowMessage("No plume heights could be calculated");
    }
    else
    {
        messageToUser.Format("Done calculating geometries. Plume height calculated on %d occasions", geometryResults.GetCount());
        ShowMessage(messageToUser);
    }
    messageToUser.Format("nFilesChecked1 = %ld, nFilesChecked2 = %ld, nCalculationsMade = %ld", nFilesChecked1, nFilesChecked2, nCalculationsMade);
    ShowMessage(messageToUser);
}

void CPostProcessing::CalculateFluxes(novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&>& evalLogFiles)
{
    CDateTime scanStartTime;
    novac::CString serial, messageToUser;
    Geometry::CPlumeHeight plumeHeight; // the altitude of the plume, in meters above sea level
    MEASUREMENT_MODE measMode;
    int channel;
    Flux::CFluxStatistics stat;

    // we keep the calculated fluxes in a list
    novac::CList <Flux::CFluxResult, Flux::CFluxResult&> calculatedFluxes;

    // Initiate the flux-calculator
    Flux::CFluxCalculator fluxCalc;

    // Loop through the list of evaluation log files. For each of them, find
    // the best available wind-speed, wind-direction and plume height and
    // calculate the flux.
    auto pos = evalLogFiles.GetHeadPosition();
    while (pos != nullptr)
    {
        // Get the name of this eval-log
        const novac::CString& evalLog = evalLogFiles.GetAt(pos).m_evalLogFile[m_userSettings.m_mainFitWindow];
        const CPlumeInScanProperty& plume = evalLogFiles.GetNext(pos).m_scanProperties;

        // if the completeness is too low then ignore this scan.
        if (plume.completeness < (m_userSettings.m_completenessLimitFlux + 0.01))
        {
            messageToUser.Format(" - Scan %s has completeness = %.2lf which is less than limit of %.2lf. Rejected!", (const char*)evalLog, plume.completeness, m_userSettings.m_completenessLimitFlux);
            ShowMessage(messageToUser);
            continue;
        }

        // Extract the date and time of day when the measurement was made
        novac::CFileUtils::GetInfoFromFileName(evalLog, scanStartTime, serial, channel, measMode);

        // If this is not a flux-measurement, then there's no point in calculating any flux for it
        if (measMode != MODE_FLUX)
            continue;

        // Extract a plume height at this time of day
        m_plumeDataBase.GetPlumeHeight(scanStartTime, plumeHeight);

        // tell the user
        messageToUser.Format("Calculating flux for measurement %s", (const char*)evalLog);
        ShowMessage(messageToUser);

        // Calculate the flux. This also takes care of writing
        // the results to file
        Flux::CFluxResult fluxResult;
        if (0 == fluxCalc.CalculateFlux(evalLog, m_windDataBase, plumeHeight, fluxResult))
        {
            calculatedFluxes.AddTail(fluxResult);
        }
    }

    // Now we can write the final fluxes to file
    ShowMessage("Writing flux log");
    WriteFluxResult_XML(calculatedFluxes);
    WriteFluxResult_Txt(calculatedFluxes);

    // Also write the statistics for the flux
    ShowMessage("Writing flux statistics");
    stat.AttachFluxList(calculatedFluxes);

    novac::CString fluxStatFileName;
    fluxStatFileName.Format("%s%c%s", m_userSettings.m_outputDirectory.c_str(), Poco::Path::separator(), "FluxStatistics.txt");
    stat.WriteFluxStat(fluxStatFileName);
}

void CPostProcessing::WriteFluxResult_XML(novac::CList <Flux::CFluxResult, Flux::CFluxResult&>& calculatedFluxes)
{
    novac::CString fluxLogFile, styleFile, wsSrc, wdSrc, phSrc, typeStr;
    CDateTime now;

    // get the current time
    now.SetToNow();

    // the name and path of the final flux-log file
    fluxLogFile.Format("%s%cFluxLog.xml", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());

    // Check if there's already a file like this, then archive it...
    Common::ArchiveFile(fluxLogFile);

    // Try to open the file
    FILE* f = fopen(fluxLogFile, "w");
    if (f == nullptr)
    {
        ShowMessage("Could not open flux log file for writing. Writing of results failed. ");
        return;
    }

    // Write the header and the starting comments
    fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(f, "<?xml-stylesheet type=\"text/xsl\" href=\"fluxresult.xsl\"?>\n");
    fprintf(f, "<!-- This is result of the flux calculations the NOVAC Post Processing Program -->\n");
    fprintf(f, "<!-- File generated on %04d.%02d.%02d at %02d:%02d:%02d -->\n\n", now.year, now.month, now.day, now.hour, now.minute, now.second);

    fprintf(f, "<NovacPPPFluxResults>\n");
    auto pos = calculatedFluxes.GetHeadPosition();
    while (pos != nullptr)
    {
        // Get the next flux result in the list
        const Flux::CFluxResult& fluxResult = calculatedFluxes.GetNext(pos);

        // extract the sources of information about wind-speed, wind-direction and plume-height
        fluxResult.m_windField.GetWindSpeedSource(wsSrc);
        fluxResult.m_windField.GetWindDirectionSource(wdSrc);
        Meteorology::MetSourceToString(fluxResult.m_plumeHeight.m_plumeAltitudeSource, phSrc);

        // write a <flux> section
        fprintf(f, "\t<flux>\n");

        fprintf(f, "\t\t<startTime>%04d.%02d.%02dT%02d:%02d:%02d</startTime>\n",
            fluxResult.m_startTime.year, fluxResult.m_startTime.month, fluxResult.m_startTime.day,
            fluxResult.m_startTime.hour, fluxResult.m_startTime.minute, fluxResult.m_startTime.second);
        fprintf(f, "\t\t<stopTime>%04d.%02d.%02dT%02d:%02d:%02d</stopTime>\n",
            fluxResult.m_stopTime.year, fluxResult.m_stopTime.month, fluxResult.m_stopTime.day,
            fluxResult.m_stopTime.hour, fluxResult.m_stopTime.minute, fluxResult.m_stopTime.second);

        fprintf(f, "\t\t<serial>%s</serial>\n", (const char*)fluxResult.m_instrument);

        // extract the instrument type
        if (fluxResult.m_instrumentType == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
        {
            typeStr.Format("heidelberg");
        }
        else
        {
            typeStr.Format("gothenburg");
        }
        fprintf(f, "\t\t<instrumentType>%s</instrumentType>\n", (const char*)typeStr);

        fprintf(f, "\t\t<value>%.2lf</value>\n", fluxResult.m_flux);

        // The judged quality of the calculated flux
        if (fluxResult.m_fluxQualityFlag == FLUX_QUALITY_GREEN)
        {
            fprintf(f, "\t\t<Quality>g</Quality>\n");
        }
        else if (fluxResult.m_fluxQualityFlag == FLUX_QUALITY_YELLOW)
        {
            fprintf(f, "\t\t<Quality>y</Quality>\n");
        }
        else
        {
            fprintf(f, "\t\t<Quality>r</Quality>\n");
        }

        // the errors
        fprintf(f, "\t\t<FluxError_Wind_kgs>%.2lf</FluxError_Wind_kgs>\n", fluxResult.m_fluxError_Wind);
        fprintf(f, "\t\t<FluxError_PlumeHeight_kgs>%.2lf</FluxError_PlumeHeight_kgs>\n", fluxResult.m_fluxError_PlumeHeight);

        // the wind speed
        fprintf(f, "\t\t<windspeed>%.2lf</windspeed>\n", fluxResult.m_windField.GetWindSpeed());
        fprintf(f, "\t\t<windspeedError>%.2lf</windspeedError>\n", fluxResult.m_windField.GetWindSpeedError());
        fprintf(f, "\t\t<windspeedSource>%s</windspeedSource>\n", (const char*)wsSrc);

        // the wind direction
        fprintf(f, "\t\t<winddirection>%.2lf</winddirection>\n", fluxResult.m_windField.GetWindDirection());
        fprintf(f, "\t\t<winddirectionError>%.2lf</winddirectionError>\n", fluxResult.m_windField.GetWindDirectionError());
        fprintf(f, "\t\t<winddirectionSource>%s</winddirectionSource>\n", (const char*)wdSrc);

        // the plume height
        fprintf(f, "\t\t<plumeheight>%.2lf</plumeheight>\n", fluxResult.m_plumeHeight.m_plumeAltitude);
        fprintf(f, "\t\t<plumeheightError>%.2lf</plumeheightError>\n", fluxResult.m_plumeHeight.m_plumeAltitudeError);
        fprintf(f, "\t\t<plumeheightSource>%s</plumeheightSource>\n", (const char*)phSrc);

        // some additional information about the scan
        fprintf(f, "\t\t<Compass>%.1lf<Compass>\n", fluxResult.m_compass);
        fprintf(f, "\t\t<ConeAngle>%.1lf<ConeAngle>\n", fluxResult.m_coneAngle);
        fprintf(f, "\t\t<Tilt>%.1lf<Tilt>\n", fluxResult.m_tilt);
        fprintf(f, "\t\t<nSpectra>%d<nSpectra>\n", fluxResult.m_numGoodSpectra);
        fprintf(f, "\t\t<PlumeCentre_1>%.1lf<PlumeCentre_1>\n", fluxResult.m_plumeCentre[0]);
        fprintf(f, "\t\t<PlumeCentre_2>%.1lf<PlumeCentre_2>\n", fluxResult.m_plumeCentre[1]);
        fprintf(f, "\t\t<PlumeCompleteness>%.2lf<PlumeCompleteness>\n", fluxResult.m_completeness);
        fprintf(f, "\t\t<ScanOffset>%.1e<ScanOffset>\n", fluxResult.m_scanOffset);

        fprintf(f, "\t</flux>\n");
    }

    fprintf(f, "</NovacPPPFluxResults>\n");

    // remember to close the file
    fclose(f);

    // ------------- we also need an xslt - file to display the output -----------------
    styleFile.Format("%s%cfluxresult.xsl", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());

    // Try to open the file
    f = fopen(styleFile, "w");
    if (f == nullptr)
    {
        return;
    }
    fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");

    fprintf(f, "<html xsl:version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\" xmlns=\"http://www.w3.org/1999/xhtml\">\n");
    fprintf(f, "<body style=\"font-family:Arial;font-size:12pt;background-color:#EEEEEE\">\n");
    fprintf(f, "\t<div style=\"background-color:white;color:black;padding:4px\">\n");
    fprintf(f, "\t\t<span style=\"font-weight:bold\">Result of flux calculation</span>\n");
    fprintf(f, "\t</div>\n");

    fprintf(f, "\t<xsl:for-each select=\"NovacPPPFluxResults/flux\">\n");
    fprintf(f, "\t<div style=\"background-color:white;color:teal;padding:4px\">\n");
    fprintf(f, "\t\t<span style=\"font-weight:bold\">Measurement from <xsl:value-of select=\"startTime\"/> to <xsl:value-of select=\"stopTime\"/> </span>\n");
    fprintf(f, "\t- <xsl:value-of select=\"value\"/> kg/s\n");
    fprintf(f, "\t</div>\n");
    fprintf(f, "\t<div style=\"margin-left:20px;margin-bottom:1em;font-size:10pt\">\n");

    fprintf(f, "\t\t<xsl:value-of select=\"description\"/>\n");
    fprintf(f, "\t\t<span style=\"font-style:italic\">\n");
    fprintf(f, "\t\t\tMade by <xsl:value-of select=\"serial\"/>\n");
    fprintf(f, "\t\t</span>\n");
    fprintf(f, "\t</div>\n");
    fprintf(f, "\t</xsl:for-each>\n");
    fprintf(f, "</body>\n");
    fprintf(f, "</html>\n");

    fclose(f);
}

void CPostProcessing::WriteFluxResult_Txt(novac::CList <Flux::CFluxResult, Flux::CFluxResult&>& calculatedFluxes)
{
    novac::CString fluxLogFile, wsSrc, wdSrc, phSrc, typeStr;
    CDateTime now;

    // get the current time
    now.SetToNow();

    // the name and path of the final flux-log file
    fluxLogFile.Format("%s%cFluxLog.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());

    // Try to open the file
    if (Filesystem::IsExistingFile(fluxLogFile))
    {
        Common::ArchiveFile(fluxLogFile);
    }

    FILE* f = fopen(fluxLogFile, "w");
    if (f == nullptr)
    {
        ShowMessage("Could not open flux log file for writing. Writing of results failed. ");
        return;
    }

    // Write the header and the starting comments
    fprintf(f, "# This is result of the flux calculations the NOVAC Post Processing Program \n");
    fprintf(f, "#   File generated on %04d.%02d.%02d at %02d:%02d:%02d \n\n", now.year, now.month, now.day, now.hour, now.minute, now.second);

    fprintf(f, "#StartTime\tStopTime\tSerial\tInstrumentType\tFlux_kgs\tFluxQuality\tFluxError_Wind_kgs\tFluxError_PlumeHeight_kgs\tWindSpeed_ms\tWindSpeedErr_ms\tWindSpeedSrc\tWindDir_deg\tWindDirErr_deg\tWindDirSrc\tPlumeHeight_m\tPlumeHeightErr_m\tPlumeHeightSrc\t");
    fprintf(f, "Compass\tConeAngle\tTilt\tnSpectra\tPlumeCentre_1\tPlumeCentre_2\tPlumeCompleteness\tScanOffset\n");

    auto pos = calculatedFluxes.GetHeadPosition();
    while (pos != nullptr)
    {
        // Get the next flux result in the list
        const Flux::CFluxResult& fluxResult = calculatedFluxes.GetNext(pos);

        // extract the instrument type
        if (fluxResult.m_instrumentType == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
        {
            typeStr.Format("heidelberg");
        }
        else
        {
            typeStr.Format("gothenburg");
        }

        // extract the sources of information about wind-speed, wind-direction and plume-height
        fluxResult.m_windField.GetWindSpeedSource(wsSrc);
        fluxResult.m_windField.GetWindDirectionSource(wdSrc);
        Meteorology::MetSourceToString(fluxResult.m_plumeHeight.m_plumeAltitudeSource, phSrc);

        // write the date and time when the measurement started and ended
        fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t",
            fluxResult.m_startTime.year, fluxResult.m_startTime.month, fluxResult.m_startTime.day,
            fluxResult.m_startTime.hour, fluxResult.m_startTime.minute, fluxResult.m_startTime.second);
        fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t",
            fluxResult.m_stopTime.year, fluxResult.m_stopTime.month, fluxResult.m_stopTime.day,
            fluxResult.m_stopTime.hour, fluxResult.m_stopTime.minute, fluxResult.m_stopTime.second);

        // the type of instrument and the serial-number
        fprintf(f, "%s\t", (const char*)fluxResult.m_instrument);
        fprintf(f, "%s\t", (const char*)typeStr);

        // The actual flux!!!
        fprintf(f, "%.2lf\t", fluxResult.m_flux);

        // The judged quality of the calculated flux
        if (fluxResult.m_fluxQualityFlag == FLUX_QUALITY_GREEN)
        {
            fprintf(f, "g\t");
        }
        else if (fluxResult.m_fluxQualityFlag == FLUX_QUALITY_YELLOW)
        {
            fprintf(f, "y\t");
        }
        else
        {
            fprintf(f, "r\t");
        }

        // the errors
        fprintf(f, "%.2lf\t", fluxResult.m_fluxError_Wind);
        fprintf(f, "%.2lf\t", fluxResult.m_fluxError_PlumeHeight);

        // the wind speed
        fprintf(f, "%.2lf\t", fluxResult.m_windField.GetWindSpeed());
        fprintf(f, "%.2lf\t", fluxResult.m_windField.GetWindSpeedError());
        fprintf(f, "%s\t", (const char*)wsSrc);

        // the wind direction
        fprintf(f, "%.2lf\t", fluxResult.m_windField.GetWindDirection());
        fprintf(f, "%.2lf\t", fluxResult.m_windField.GetWindDirectionError());
        fprintf(f, "%s\t", (const char*)wdSrc);

        // the plume height
        fprintf(f, "%.2lf\t", fluxResult.m_plumeHeight.m_plumeAltitude);
        fprintf(f, "%.2lf\t", fluxResult.m_plumeHeight.m_plumeAltitudeError);
        fprintf(f, "%s\t", (const char*)phSrc);

        // write additional information about the scan
        fprintf(f, "%.1lf\t", fluxResult.m_compass);
        fprintf(f, "%.1lf\t", fluxResult.m_coneAngle);
        fprintf(f, "%.1lf\t", fluxResult.m_tilt);
        fprintf(f, "%d\t", fluxResult.m_numGoodSpectra);
        fprintf(f, "%.1lf\t", fluxResult.m_plumeCentre[0]);
        fprintf(f, "%.1lf\t", fluxResult.m_plumeCentre[1]);
        fprintf(f, "%.2lf\t", fluxResult.m_completeness);
        fprintf(f, "%.1e\n", fluxResult.m_scanOffset);

    }

    // remember to close the file
    fclose(f);
}

void CPostProcessing::WriteCalculatedGeometriesToFile(novac::CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*>& geometryResults)
{
    if (geometryResults.GetCount() == 0)
        return; // nothing to write...

    FILE* f = nullptr;
    novac::CString geomLogFile;
    geomLogFile.Format("%s%cGeometryLog.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());

    if (Filesystem::IsExistingFile(geomLogFile))
    {
        f = fopen(geomLogFile, "a");
        if (f == nullptr)
        {
            ShowMessage("Could not open geometry log file for writing. Writing of results failed. ");
            return;
        }
    }
    else
    {
        f = fopen(geomLogFile, "w");
        if (f == nullptr)
        {
            ShowMessage("Could not open geometry log file for writing. Writing of results failed. ");
            return;
        }
        fprintf(f, "Date\tTime\tDifferenceInStartTime_minutes\tInstrument1\tInstrument2\tPlumeAltitude_masl\tPlumeHeightError_m\tWindDirection_deg\tWindDirectionError_deg\tPlumeCentre1_deg\tPlumeCentreError1_deg\tPlumeCentre2_deg\tPlumeCentreError2_deg\n");
    }

    auto pos = geometryResults.GetHeadPosition();
    while (pos != nullptr)
    {
        Geometry::CGeometryResult* result = geometryResults.GetNext(pos);
        // write the file
        if (result->m_calculationType == Meteorology::MET_GEOMETRY_CALCULATION)
        {
            fprintf(f, "%04d.%02d.%02d\t", result->m_averageStartTime.year, result->m_averageStartTime.month, result->m_averageStartTime.day);
            fprintf(f, "%02d:%02d:%02d\t", result->m_averageStartTime.hour, result->m_averageStartTime.minute, result->m_averageStartTime.second);
            fprintf(f, "%.1lf\t", result->m_startTimeDifference / 60.0);
            fprintf(f, "%s\t%s\t", (const char*)result->m_instr1, (const char*)result->m_instr2);
            fprintf(f, "%.0lf\t%.0lf\t", result->m_plumeAltitude, result->m_plumeAltitudeError);
            fprintf(f, "%.0lf\t%.0lf\t", result->m_windDirection, result->m_windDirectionError);

            fprintf(f, "%.1f\t%.1f\t", result->m_plumeCentre1, result->m_plumeCentreError1);
            fprintf(f, "%.1f\t%.1f\n", result->m_plumeCentre2, result->m_plumeCentreError2);
        }
        else
        {
            fprintf(f, "%04d.%02d.%02d\t", result->m_averageStartTime.year, result->m_averageStartTime.month, result->m_averageStartTime.day);
            fprintf(f, "%02d:%02d:%02d\t", result->m_averageStartTime.hour, result->m_averageStartTime.minute, result->m_averageStartTime.second);
            fprintf(f, "0\t");
            fprintf(f, "%s\t\t", (const char*)result->m_instr1);
            fprintf(f, "%.0lf\t%.0lf\t", result->m_plumeAltitude, result->m_plumeAltitudeError);
            fprintf(f, "%.0lf\t%.0lf\t", result->m_windDirection, result->m_windDirectionError);

            fprintf(f, "%.1f\t%.1f\t", result->m_plumeCentre1, result->m_plumeCentreError1);
            fprintf(f, "0\t0\n");
        }
    }
    fclose(f);
}

void CPostProcessing::InsertCalculatedGeometriesIntoDataBase(novac::CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*>& geometryResults)
{
    Meteorology::CWindField windField;
    CDateTime validFrom, validTo;
    Configuration::CInstrumentLocation location;

    auto pos = geometryResults.GetHeadPosition();
    while (pos != nullptr)
    {
        Geometry::CGeometryResult* result = geometryResults.GetNext(pos);

        if (result->m_plumeAltitude > 0.0)
        {
            // insert the plume height into the plume height database
            this->m_plumeDataBase.InsertPlumeHeight(*result);
        }

        if (result->m_windDirection > NOT_A_NUMBER)
        {
            try
            {
                // get the location of the instrument at the time of the measurement
                location = m_setup.GetInstrumentLocation(result->m_instr1.std_str(), result->m_averageStartTime);

                // get the time-interval that the measurement is valid for
                validFrom = CDateTime(result->m_averageStartTime);
                validFrom.Decrement(m_userSettings.m_calcGeometryValidTime);
                validTo = CDateTime(result->m_averageStartTime);
                validTo.Increment(m_userSettings.m_calcGeometryValidTime);

                // insert the wind-direction into the wind database
                m_windDataBase.InsertWindDirection(validFrom, validTo, result->m_windDirection, result->m_windDirectionError, result->m_calculationType, nullptr);
            }
            catch (PPPLib::NotFoundException& ex)
            {
                ShowMessage(ex.message);
            }
        }
    }
}

void CPostProcessing::CalculateDualBeamWindSpeeds(novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&>& evalLogs)
{
    novac::CList <novac::CString, novac::CString&> masterList; // list of wind-measurements from the master channel
    novac::CList <novac::CString, novac::CString&> slaveList;  // list of wind-measurements from the slave channel
    novac::CList <novac::CString, novac::CString&> heidelbergList;  // list of wind-measurements from the Heidelbergensis
    CDateTime validFrom, validTo;

    novac::CString serial, serial2, fileName, fileName2, nonsenseString;
    novac::CString userMessage, windLogFile;
    CDateTime startTime, startTime2;
    int channel, channel2, nWindMeasFound = 0;
    MEASUREMENT_MODE meas_mode, meas_mode2;
    Configuration::CInstrumentLocation location;
    WindSpeedMeasurement::CWindSpeedCalculator calculator;
    Geometry::CPlumeHeight plumeHeight;
    Meteorology::CWindField windField, oldWindField;

    // -------------------------------- step 1. -------------------------------------
    // search through 'evalLogs' for dual-beam measurements from master and from slave
    auto logPosition = evalLogs.GetHeadPosition();
    while (logPosition != nullptr)
    {
        try
        {
            const novac::CString& fileNameAndPath = evalLogs.GetNext(logPosition).m_evalLogFile[m_userSettings.m_mainFitWindow];

            // to know the start-time of the measurement, we need to 
            // extract just the file-name, i.e. remove the path
            fileName = novac::CString(fileNameAndPath);
            Common::GetFileName(fileName);

            novac::CFileUtils::GetInfoFromFileName(fileName, startTime, serial, channel, meas_mode);

            if (meas_mode == MODE_WINDSPEED)
            {
                ++nWindMeasFound;

                // first check if this is a heidelberg instrument
                location = m_setup.GetInstrumentLocation(serial.std_str(), startTime);

                if (location.m_instrumentType == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
                {
                    // this is a heidelberg instrument
                    heidelbergList.AddTail(novac::CString(fileNameAndPath));
                }
                else
                {
                    // this is a gothenburg instrument
                    if (channel == 0)
                    {
                        masterList.AddTail(novac::CString(fileNameAndPath));
                    }
                    else if (channel == 1)
                    {
                        slaveList.AddTail(novac::CString(fileNameAndPath));
                    }
                }
            }
        }
        catch (PPPLib::NotFoundException& ex)
        {
            ShowMessage(ex.message);
        }
    }
    if (nWindMeasFound == 0)
    {
        ShowMessage("No dual-beam wind speed measurements found.");
        return; // if nothing was found...
    }

    userMessage.Format("%d dual-beam wind speed measurements found. Calculating wind-speeds", nWindMeasFound);
    ShowMessage(userMessage);

    // Create the dual-beam log-file
    windLogFile.Format("%s%cDualBeamLog.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    calculator.WriteWindSpeedLogHeader(windLogFile);


    // -------------------------------- step 2. -------------------------------------
    // loop through each of the measurements from the heidelberg instruments
    // and calculate the wind speed for each measurement
    auto instrPos = heidelbergList.GetHeadPosition();
    while (instrPos != nullptr)
    {
        try
        {
            const novac::CString& fileNameAndPath = heidelbergList.GetNext(instrPos);

            // to know the start-time of the measurement, we need to 
            // extract just the file-name, i.e. remove the path
            fileName = novac::CString(fileNameAndPath);
            Common::GetFileName(fileName);
            novac::CFileUtils::GetInfoFromFileName(fileName, startTime, serial, channel, meas_mode);

            // Get the plume height at the time of the measurement
            m_plumeDataBase.GetPlumeHeight(startTime, plumeHeight);

            // Get the location of the instrument at the time of the measurement
            location = m_setup.GetInstrumentLocation(serial.std_str(), startTime);

            // calculate the speed of the wind at the time of the measurement
            if (0 == calculator.CalculateWindSpeed(fileNameAndPath, nonsenseString, location, plumeHeight, windField))
            {
                // append the results to file
                calculator.AppendResultToFile(windLogFile, startTime, location, plumeHeight, windField);

                // insert the newly calculated wind-speed into the database
                if (windField.GetWindSpeedError() > m_userSettings.m_dualBeam_MaxWindSpeedError)
                {
                    userMessage.Format("-Calculated a wind-speed of %.1lf +- %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Error too large, measurement discarded.", windField.GetWindSpeed(), windField.GetWindSpeedError(),
                        startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);
                }
                else
                {
                    // tell the user...
                    userMessage.Format("+Calculated a wind-speed of %.1lf +- %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Measurement accepted", windField.GetWindSpeed(), windField.GetWindSpeedError(),
                        startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);

                    // get the time-interval that the measurement is valid for
                    windField.GetValidTimeFrame(validFrom, validTo);

                    // insert the new wind speed into the database
                    m_windDataBase.InsertWindSpeed(validFrom, validTo, windField.GetWindSpeed(), windField.GetWindSpeedError(), Meteorology::MET_DUAL_BEAM_MEASUREMENT, nullptr);
                }
                ShowMessage(userMessage);
            }
            else
            {
                userMessage.Format("Failed to calculate wind speed from measurement: %s", (const char*)fileName);
                ShowMessage(userMessage);
            }
        }
        catch (PPPLib::NotFoundException& ex)
        {
            ShowMessage(ex.message);
        }
    }

    // -------------------------------- step 3. -------------------------------------
    // loop through each of the measurements from a master-channel and try to match them with a measurement
    // from a slave channel...
    auto masterPos = masterList.GetHeadPosition();
    while (masterPos != nullptr)
    {
        const novac::CString& fileNameAndPath = masterList.GetNext(masterPos);

        // extract just the file-name, i.e. remove the path
        fileName = novac::CString(fileNameAndPath);
        Common::GetFileName(fileName);

        novac::CFileUtils::GetInfoFromFileName(fileName, startTime, serial, channel, meas_mode);

        // now check if we can match this one with a file in the slave-channel
        auto pos2 = slaveList.GetHeadPosition();
        while (pos2 != nullptr)
        {
            const novac::CString& fileNameAndPath2 = slaveList.GetNext(pos2);

            // extract just the file-name, i.e. remove the path
            fileName2 = novac::CString(fileNameAndPath2);
            Common::GetFileName(fileName2);

            novac::CFileUtils::GetInfoFromFileName(fileName2, startTime2, serial2, channel2, meas_mode2);

            if (Equals(serial, serial2) && (startTime == startTime2))
            {
                // we have found a match!!!

                // Get the plume height at the time of the measurement
                m_plumeDataBase.GetPlumeHeight(startTime, plumeHeight);

                // Get the location of the instrument at the time of the measurement
                location = m_setup.GetInstrumentLocation(serial.std_str(), startTime);

                // calculate the speed of the wind at the time of the measurement
                if (0 == calculator.CalculateWindSpeed(fileNameAndPath, fileNameAndPath2, location, plumeHeight, windField))
                {
                    // append the results to file
                    calculator.AppendResultToFile(windLogFile, startTime, location, plumeHeight, windField);

                    // insert the newly calculated wind-speed into the database
                    if (windField.GetWindSpeedError() > m_userSettings.m_dualBeam_MaxWindSpeedError)
                    {
                        userMessage.Format("-Calculated a wind-speed of %.1lf +- %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Error too large, measurement discarded.", windField.GetWindSpeed(), windField.GetWindSpeedError(),
                            startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);
                    }
                    else
                    {
                        // tell the user...
                        userMessage.Format("+Calculated a wind-speed of %.1lf +- %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Measurement accepted", windField.GetWindSpeed(), windField.GetWindSpeedError(),
                            startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);

                        windField.GetValidTimeFrame(validFrom, validTo);

                        // insert the new wind speed into the database
                        m_windDataBase.InsertWindSpeed(validFrom, validTo, windField.GetWindSpeed(), windField.GetWindSpeedError(), Meteorology::MET_DUAL_BEAM_MEASUREMENT, nullptr);
                    }
                    ShowMessage(userMessage);

                }
                else
                {
                    userMessage.Format("Failed to calculate wind speed from measurement: %s", (const char*)fileName);
                    ShowMessage(userMessage);
                }
            }
        }
    }
}

void CPostProcessing::SortEvaluationLogs(novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&>& evalLogs)
{
    novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&> left;
    novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&> right;

    // If this list consists of only one element, then we're done
    if (evalLogs.GetCount() <= 1)
        return;

    // Divide the list into two, and sort each one of them
    auto pos1 = evalLogs.GetHeadPosition();
    int index = 0;
    while (pos1 != nullptr)
    {
        Evaluation::CExtendedScanResult& log = evalLogs.GetNext(pos1);
        if (index % 2 == 0)
            left.AddTail(log);
        else
            right.AddTail(log);
        ++index;
    }

    SortEvaluationLogs(left);
    SortEvaluationLogs(right);

    // Merge the two lists into one, sorted list
    evalLogs.RemoveAll();
    auto pos3 = left.GetHeadPosition();
    auto pos2 = right.GetHeadPosition();
    while (pos3 != nullptr && pos2 != nullptr)
    {
        Evaluation::CExtendedScanResult& log1 = left.GetAt(pos3);
        Evaluation::CExtendedScanResult& log2 = right.GetAt(pos2);

        if (log2.m_startTime < log1.m_startTime)
        {
            evalLogs.AddTail(log2);
            right.GetNext(pos2);
        }
        else
        {
            evalLogs.AddTail(log1);
            left.GetNext(pos3);
        }
    }

    while (pos3 != nullptr)
    {
        evalLogs.AddTail(left.GetNext(pos3));
    }

    while (pos2 != nullptr)
    {
        evalLogs.AddTail(right.GetNext(pos2));
    }

    return;
}

void CPostProcessing::UploadResultsToFTP()
{
    Communication::CFTPServerConnection connection;
    novac::CString fileName;

    // Generate a list with all the files we want to upload.
    novac::CList <novac::CString, novac::CString&> fileList;

    // 1. the geometry log file
    fileName.Format("%s%cGeometryLog.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    fileList.AddTail(fileName);

    // 2. the generated wind field
    fileName.Format("%s%cGeneratedWindField.wxml", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    fileList.AddTail(fileName);

    // 3. the generated flux log
    fileName.Format("%s%cFluxLog.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    fileList.AddTail(fileName);

    // upload the files
    connection.UploadResults(m_userSettings.m_FTPDirectory, m_userSettings.m_FTPUsername, m_userSettings.m_FTPPassword, fileList);

    return;
}

bool CPostProcessing::ConvolveReference(novac::CReferenceFile& ref, const novac::CString& instrumentSerial)
{
    // Make sure the high-res section do exist.
    if (!IsExistingFile(ref.m_crossSectionFile))
    {
        std::string fullPath = Filesystem::GetAbsolutePathFromRelative(ref.m_crossSectionFile, this->m_exePath);
        if (Filesystem::IsExistingFile(fullPath))
        {
            ref.m_crossSectionFile = fullPath;
        }
        else
        {
            novac::CString errorMessage;
            errorMessage.Format("Cannot find given cross section file '%s.", ref.m_crossSectionFile.c_str());
            ShowMessage(errorMessage);
            return false; // failed to find the file
        }
    }

    // Make sure the slit-function do exist.
    if (!IsExistingFile(ref.m_slitFunctionFile))
    {
        std::string fullPath = Filesystem::GetAbsolutePathFromRelative(ref.m_slitFunctionFile, this->m_exePath);
        if (Filesystem::IsExistingFile(fullPath))
        {
            ref.m_slitFunctionFile = fullPath;
        }
        else
        {
            novac::CString errorMessage;
            errorMessage.Format("Cannot find given slit-function file '%s.", ref.m_slitFunctionFile.c_str());
            ShowMessage(errorMessage);
            return false; // failed to find the file
        }
    }

    // Make sure the wavelength calibration do exist.
    if (!IsExistingFile(ref.m_wavelengthCalibrationFile))
    {
        std::string fullPath = Filesystem::GetAbsolutePathFromRelative(ref.m_wavelengthCalibrationFile, this->m_exePath);
        if (Filesystem::IsExistingFile(fullPath))
        {
            ref.m_wavelengthCalibrationFile = fullPath;
        }
        else
        {
            novac::CString errorMessage;
            errorMessage.Format("Cannot find given wavelength calibration file '%s.", ref.m_wavelengthCalibrationFile.c_str());
            ShowMessage(errorMessage);
            return false; // failed to find the file
        }
    }

    // Now do the convolution
    if (ref.ConvolveReference())
    {
        return false;
    }

    // Save the resulting reference, for reference...
    novac::CString tempFile;
    tempFile.Format("%s%s_%s.xs", (const char*)m_userSettings.m_tempDirectory, (const char*)instrumentSerial, ref.m_specieName.c_str());
    SaveCrossSectionFile(tempFile.ToStdString(), *ref.m_data);

    return true;
}

void CPostProcessing::LocateEvaluationLogFiles(const novac::CString& directory, novac::CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&>& evaluationLogFiles)
{
    std::vector<std::string> evalLogFiles;

    novac::CString messageToUser;
    messageToUser.Format("Searching for evaluation log - files in directory %s", (const char*)m_userSettings.m_outputDirectory);
    ShowMessage(messageToUser);

    const bool includeSubDirs = (m_userSettings.m_includeSubDirectories_Local > 0);
    Filesystem::FileSearchCriterion limits;
    limits.startTime = m_userSettings.m_fromDate;
    limits.endTime = m_userSettings.m_toDate;
    limits.fileExtension = "_flux.txt";
    Filesystem::SearchDirectoryForFiles(directory, includeSubDirs, evalLogFiles, &limits);


    messageToUser.Format("%d Evaluation log files found, starting reading", evalLogFiles.size());
    ShowMessage(messageToUser);

    size_t nofFailedLogReads = 0;

    for (std::string& f : evalLogFiles)
    {
        int channel;
        CDateTime startTime;
        MEASUREMENT_MODE mode;
        novac::CString serial;
        novac::CFileUtils::GetInfoFromFileName(f, startTime, serial, channel, mode);

        Evaluation::CExtendedScanResult result;
        result.m_evalLogFile[0] = f;
        result.m_startTime = startTime;

        FileHandler::CEvaluationLogFileHandler logReader;
        logReader.m_evaluationLog = novac::CString(f);
        if (RETURN_CODE::SUCCESS != logReader.ReadEvaluationLog() || logReader.m_scan.size() == 0)
        {
            ++nofFailedLogReads;
            continue;
        }
        Evaluation::CScanResult scanResult = logReader.m_scan[0];

        evaluationLogFiles.AddTail(result);
    }

    messageToUser.Format("%d Evaluation log files read successfully.", evalLogFiles.size() - nofFailedLogReads);
    ShowMessage(messageToUser);
}
