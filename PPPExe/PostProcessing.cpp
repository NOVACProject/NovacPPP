#include "stdafx.h"
#include "PostProcessing.h"
#include "Common/Common.h"
#include <SpectralEvaluation/File/File.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <list>
#include <sstream>
#include <thread>

#include <PPPLib/PostProcessingUtils.h>

// the PostEvaluationController takes care of the DOAS evaluations
#include <PPPLib/Evaluation/PostEvaluationController.h>

// The PostCalibration takes care of the instrument calibrations.
#include <PPPLib/Calibration/PostCalibration.h>
#include <PPPLib/Calibration/PostCalibrationStatistics.h>

// The FluxCalculator takes care of calculating the fluxes
#include <PPPLib/Flux/FluxCalculator.h>

// The Stratospherecalculator takes care of calculating Stratospheric VCD's
#include "Stratosphere/StratosphereCalculator.h"

// The flux CFluxStatistics takes care of the statistical part of the fluxes
#include <PPPLib/Flux/FluxStatistics.h>

// We also need to read the evaluation-log files
#include <PPPLib/File/EvaluationLogFileHandler.h>

#include <PPPLib/WindMeasurement/WindSpeedCalculator.h>
#include <PPPLib/Meteorology/XMLWindFileReader.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/VolcanoInfo.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <SpectralEvaluation/ThreadUtils.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Calibration/StandardCrossSectionSetup.h>
#include <SpectralEvaluation/VectorUtils.h>

// we need to be able to download data from the FTP-server
#include <PPPLib/Communication/FTPServerConnection.h>

#include <Poco/DirectoryIterator.h>
#include <Poco/Exception.h>

#undef min
#undef max

extern novac::CVolcanoInfo  g_volcanoes;   // <-- A list of all known volcanoes

using namespace novac;


// this is the working-thread that takes care of evaluating a portion of the scans
void EvaluateScansThread(
    ILogger& log,
    const Configuration::CNovacPPPConfiguration& setup,
    const Configuration::CUserConfiguration& userSettings,
    const CContinuationOfProcessing& continuation,
    CPostProcessingStatistics& processingStats);

CPostProcessing::CPostProcessing(ILogger& logger, Configuration::CNovacPPPConfiguration setup, Configuration::CUserConfiguration userSettings, const CContinuationOfProcessing& continuation)
    : m_log(logger), m_setup(setup), m_userSettings(userSettings), m_continuation(continuation)
{
    assert(m_setup.m_executableDirectory.size() > 3); // This should be set
}

void CPostProcessing::DoPostProcessing_Flux()
{
    std::vector<Evaluation::CExtendedScanResult> evaluatedScanResult;
    novac::CString messageToUser, windFileName;

    novac::LogContext context;

    m_log.Information("--- Prepairing to perform Flux Calculations --- ");
    // --------------- PREPARING FOR THE PROCESSING -----------

    // Checks that the evaluation is ok and that all settings makes sense.
    CheckProcessingSettings();

    // Prepare for the flux-calculation by reading in the wind-field
    ReadWindField(context);

    // Prepare for the flux-calculation by compiling a set of pausible plume heights
    PreparePlumeHeights(context);

    // --------------- DOING THE EVALUATIONS -----------

    if (m_userSettings.m_doEvaluations)
    {
        // Prepare for the evaluation by reading in the reference files
        PrepareEvaluation(m_log, m_userSettings.m_tempDirectory.std_str(), m_setup);

        // 1. Find all .pak files in the directory.
        m_log.Information(context, "--- Locating Pak Files --- ");
        std::vector<std::string> pakFileList;
        if (m_userSettings.m_LocalDirectory.GetLength() > 3)
        {
            novac::LogContext localContext = context.With(novac::LogContext::Directory, m_userSettings.m_LocalDirectory.std_str());
            m_log.Information(localContext, "Searching for .pak files");

            const bool includeSubDirs = (m_userSettings.m_includeSubDirectories_Local > 0);
            Filesystem::FileSearchCriterion limits;
            limits.startTime = m_userSettings.m_fromDate;
            limits.endTime = m_userSettings.m_toDate;
            limits.fileExtension = ".pak";
            Filesystem::SearchDirectoryForFiles(m_userSettings.m_LocalDirectory, includeSubDirs, pakFileList, &limits);

            std::stringstream msg;
            msg << pakFileList.size() << " .pak files found";
            m_log.Information(localContext, msg.str());
        }

        if (m_userSettings.m_FTPDirectory.GetLength() > 9)
        {
            novac::LogContext localContext = context.With("ftpDirectory", m_userSettings.m_FTPDirectory.std_str());
            m_log.Information(localContext, "Searching for .pak files on Ftp server");

            CheckForSpectraOnFTPServer(localContext, pakFileList);
        }
        if (pakFileList.size() == 0)
        {
            m_log.Information(context, "No spectrum files found. Exiting");
            return;
        }

        // Evaluate the scans. This at the same time generates a list of evaluation-log
        // files with the evaluated results
        m_log.Information(context, "--- Running Evaluations --- ");
        EvaluateScans(pakFileList, evaluatedScanResult);
        messageToUser.Format("%d evaluation log files accepted", evaluatedScanResult.size());
        m_log.Information(context, messageToUser.std_str());
    }
    else
    {
        // Don't evaluate the scans, just read the log-files and calculate fluxes from there.
        evaluatedScanResult = LocateEvaluationLogFiles(context, m_userSettings.m_outputDirectory.std_str());
    }

    // Sort the evaluation-logs in order of increasing start-time, this to make
    // the looking for matching files in 'CalculateGeometries' faster
    m_log.Information("Evaluation done. Sorting the evaluation results");
    SortEvaluationLogs(evaluatedScanResult);
    m_log.Information(context, "Sort done.");

    // 3. Loop through list with output text files from evaluation and calculate the geometries
    std::vector<Geometry::CGeometryResult> geometryResults;
    CalculateGeometries(context, evaluatedScanResult, geometryResults);

    // 4.1 write the calculations to file, for later checking or other uses...
    WriteCalculatedGeometriesToFile(context, geometryResults);

    // 4.2 Insert the calculated geometries into the plume height database
    InsertCalculatedGeometriesIntoDatabase(context, geometryResults);

    // 5. Calculate the wind-speeds from the wind-speed measurements
    //  the plume heights are taken from the database
    CalculateDualBeamWindSpeeds(context, evaluatedScanResult);

    // 6. Calculate flux from evaluation text files
    CalculateFluxes(context, evaluatedScanResult);

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
        UploadResultsToFTP(context);
    }
}

void CPostProcessing::DoPostProcessing_InstrumentCalibration()
{
    novac::LogContext context("mode", "instrumentCalibration");

    m_log.Information(context, "--- Prepairing to perform Instrument Calibrations --- ");

    m_log.Information("--- Validating settings --- ");
    if (CheckInstrumentCalibrationSettings())
    {
        m_log.Information(context, "Exiting post processing");
        return;
    }

    novac::StandardCrossSectionSetup standardCrossSections{ m_setup.m_executableDirectory };
    if (standardCrossSections.NumberOfReferences() == 0)
    {
        m_log.Information(context, "The StandardReferences folder is either missing or contains no references. No references can be created.");
        return;
    }
    else if (!IsExistingFile(standardCrossSections.FraunhoferReferenceFileName()))
    {
        m_log.Information(context, "Cannot locate the Fraunhofer reference in the StandardReferences folder. The file is either missing or path is invalid. No Fraunhofer reference can be created.");
    }


    // 1. Find all .pak files in the directory.
    m_log.Information(context, "--- Locating Pak Files --- ");
    std::vector<std::string> pakFileList;
    if (m_userSettings.m_LocalDirectory.GetLength() > 3)
    {
        novac::LogContext localContext = context.With(novac::LogContext::Directory, m_userSettings.m_LocalDirectory.std_str());
        m_log.Information(localContext, "Searching for .pak files");

        const bool includeSubDirs = (m_userSettings.m_includeSubDirectories_Local > 0);
        Filesystem::FileSearchCriterion limits;
        limits.startTime = m_userSettings.m_fromDate;
        limits.endTime = m_userSettings.m_toDate;
        limits.fileExtension = ".pak";
        Filesystem::SearchDirectoryForFiles(m_userSettings.m_LocalDirectory, includeSubDirs, pakFileList, &limits);

        std::stringstream msg;
        msg << pakFileList.size() << " .pak files found";
        m_log.Information(localContext, msg.str());
    }

    if (m_userSettings.m_FTPDirectory.GetLength() > 9)
    {
        novac::LogContext localContext = context.With("ftpDirectory", m_userSettings.m_FTPDirectory.std_str());
        m_log.Information(localContext, "Searching for .pak files on Ftp server");

        CheckForSpectraOnFTPServer(localContext, pakFileList);
    }
    if (pakFileList.size() == 0)
    {
        m_log.Information(context, "No spectrum files found. Exiting");
        return;
    }

    m_log.Information(context, "--- Running Calibrations --- ");

    novac::CPostCalibration calibrationController{ standardCrossSections, m_setup, m_userSettings, m_log };
    novac::CPostCalibrationStatistics calibrationStatistics;

    // Unlike other parts of the NovacPPP, this function is intended to be single threaded. The reason for this is that the
    //  classes which are called upon are strongly threaded themselves and further threading will not help the performance.
    int numberOfCalibrations = calibrationController.RunInstrumentCalibration(pakFileList, calibrationStatistics);

    {
        CString messageToUser;
        messageToUser.Format("%d instrument calibrations performed", numberOfCalibrations);
        m_log.Information(context, messageToUser.std_str());
    }
}

void CPostProcessing::DoPostProcessing_Strat()
{
    std::vector<Evaluation::CExtendedScanResult> evalLogFiles;
    novac::CString messageToUser, statFileName;
    Stratosphere::CStratosphereCalculator strat(m_log, m_setup, m_userSettings);

    // --------------- PREPARING FOR THE PROCESSING -----------

    novac::LogContext context("mode", "stratosphere");

    // Checks that the evaluation is ok and that all settings makes sense.
    CheckProcessingSettings();

    // Prepare for the evaluation by reading in the reference files
    PrepareEvaluation(m_log, m_userSettings.m_tempDirectory.std_str(), m_setup);

    // --------------- DOING THE PROCESSING -----------

    // 1. Find all .pak files in the directory.
    std::vector<std::string> pakFileList;
    if (m_userSettings.m_LocalDirectory.GetLength() > 3)
    {
        novac::LogContext localContext = context.With(novac::LogContext::Directory, m_userSettings.m_LocalDirectory.std_str());
        m_log.Information(localContext, "Searching for .pak files");

        const bool includeSubDirs = (m_userSettings.m_includeSubDirectories_Local > 0);
        Filesystem::FileSearchCriterion limits;
        limits.startTime = m_userSettings.m_fromDate;
        limits.endTime = m_userSettings.m_toDate;
        limits.fileExtension = ".pak";
        Filesystem::SearchDirectoryForFiles(m_userSettings.m_LocalDirectory, includeSubDirs, pakFileList, &limits);

        std::stringstream msg;
        msg << pakFileList.size() << " .pak files found";
        m_log.Information(localContext, msg.str());
    }
    if (m_userSettings.m_FTPDirectory.GetLength() > 9)
    {
        novac::LogContext localContext = context.With("ftpDirectory", m_userSettings.m_FTPDirectory.std_str());
        m_log.Information(localContext, "Searching for .pak files on Ftp server");

        CheckForSpectraOnFTPServer(localContext, pakFileList);
    }
    if (pakFileList.size() == 0)
    {
        m_log.Information(context, "No spectrum files found. Exiting");
        return;
    }

    // Evaluate the scans. This at the same time generates a list of evaluation-log
    // files with the evaluated results
    EvaluateScans(pakFileList, evalLogFiles);
    messageToUser.Format("%d evaluation log files accepted", evalLogFiles.size());
    m_log.Information(messageToUser.std_str());

    // Sort the evaluation-logs in order of increasing start-time, this to make
    // the looking for matching files in 'CalculateGeometries' faster
    m_log.Information("Evaluation done. Sorting the evaluation log files");
    SortEvaluationLogs(evalLogFiles);
    m_log.Information("Sort done.");

    // m_log.Information("Calculating stratospheric columns");
    // strat.CalculateVCDs(evalLogFiles);
    // m_log.Information("Calculation Done");

    // 7. Write the statistics
    statFileName.Format("%s%cProcessingStatistics.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    Common::ArchiveFile(statFileName);
    m_processingStats.WriteStatToFile(statFileName);
}

void CPostProcessing::CheckForSpectraOnFTPServer(novac::LogContext context, std::vector<std::string>& fileList)
{
    Communication::CFTPServerConnection serverDownload(m_log, m_userSettings);

    int ret = serverDownload.DownloadDataFromFTP(
        context,
        m_userSettings.m_FTPDirectory,
        m_userSettings.m_FTPUsername,
        m_userSettings.m_FTPPassword,
        fileList);

    if (ret == 0)
    {
        m_log.Information(context, "Successfully downloaded all data files.");
    }
    else
    {
        m_log.Information(context, "Error happened when downloading data from FTP.");
    }
}

novac::GuardedList<std::string> s_pakFilesRemaining;
novac::GuardedList<Evaluation::CExtendedScanResult> s_evalLogs;

volatile size_t s_nFilesToProcess;

void CPostProcessing::EvaluateScans(
    const std::vector<std::string>& pakFileList,
    std::vector<Evaluation::CExtendedScanResult>& evalLogFiles)
{
    s_nFilesToProcess = pakFileList.size();
    novac::CString messageToUser;

    // share the list of pak-files with the other functions around here
    for (const std::string& file : pakFileList)
    {
        s_pakFilesRemaining.AddItem(file);
    }

    // Keep the user informed about what we're doing
    messageToUser.Format("%ld spectrum files found. Begin evaluation using %d threads.", s_nFilesToProcess, m_userSettings.m_maxThreadNum);
    m_log.Information(messageToUser.std_str());

    // start the threads
    std::vector<std::thread> evalThreads(m_userSettings.m_maxThreadNum);
    for (unsigned int threadIdx = 0; threadIdx < m_userSettings.m_maxThreadNum; ++threadIdx)
    {
        std::thread t(EvaluateScansThread, std::ref(m_log), std::cref(m_setup), std::cref(m_userSettings), std::cref(m_continuation), std::ref(m_processingStats));
        evalThreads[threadIdx] = std::move(t);
    }

    // make sure that all threads have time to finish before we say that we're ready
    for (unsigned int threadIdx = 0; threadIdx < m_userSettings.m_maxThreadNum; ++threadIdx)
    {
        evalThreads[threadIdx].join();
    }

    // copy out the result
    s_evalLogs.CopyTo(evalLogFiles);

    messageToUser.Format("All %ld scans evaluated. Final number of results: %ld", s_nFilesToProcess, evalLogFiles.size());
    m_log.Information(messageToUser.std_str());
}

void EvaluateScansThread(
    ILogger& log,
    const Configuration::CNovacPPPConfiguration& setup,
    const Configuration::CUserConfiguration& userSettings,
    const CContinuationOfProcessing& continuation,
    CPostProcessingStatistics& processingStats)
{
    std::string fileName;

    // create a new CPostEvaluationController
    Evaluation::CPostEvaluationController eval{ log, setup, userSettings, continuation, processingStats };

    // while there are more .pak-files
    while (s_pakFilesRemaining.PopFront(fileName))
    {
        novac::LogContext context(novac::LogContext::FileName, novac::GetFileName(fileName));

        // evaluate the .pak-file in all the specified fit-windows and retrieve the name of the 
        // eval-logs. If any of the fit-windows fails then the scan is not inserted.
        bool evaluationSucceeded = true;
        Evaluation::CExtendedScanResult combinedResult;
        try
        {
            for (size_t fitWindowIndex = 0; fitWindowIndex < userSettings.m_nFitWindowsToUse; ++fitWindowIndex)
            {
                auto result = eval.EvaluateScan(fileName, userSettings.m_fitWindowsToUse[fitWindowIndex]);

                if (result == nullptr)
                {
                    evaluationSucceeded = false;
                    break;
                }

                combinedResult.m_pakFile = result->m_pakFile;
                combinedResult.m_startTime = result->m_startTime;
                combinedResult.m_instrumentSerial = result->m_instrumentSerial;
                combinedResult.m_measurementMode = result->m_measurementMode;
                combinedResult.m_evalLogFile.push_back(result->m_evalLogFile.front());
                combinedResult.m_fitWindowName.push_back(result->m_fitWindowName.front());
                if (fitWindowIndex == userSettings.m_mainFitWindow)
                {
                    combinedResult.m_scanProperties = result->m_scanProperties;
                }
            }
        }
        catch (std::exception& ex)
        {
            log.Error(context, ex.what());
            evaluationSucceeded = false;
        }

        if (evaluationSucceeded)
        {
            // If we made it this far then the measurement is ok, insert it into the list!
            s_evalLogs.AddItem(combinedResult);
            processingStats.InsertAcception(combinedResult.m_instrumentSerial);

            log.Information(context, "Inserted scan into list of evaluation logs");
        }
        else
        {
            log.Information(context, "No flux calculated for scan.");
        }
    }
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
        m_log.Error(message.str());
        returnCode = 1;
    }

    if (m_userSettings.m_calibrationIntervalTimeOfDayLow < 0 || m_userSettings.m_calibrationIntervalTimeOfDayLow > 86400 ||
        m_userSettings.m_calibrationIntervalTimeOfDayHigh < 0 || m_userSettings.m_calibrationIntervalTimeOfDayHigh > 86400)
    {
        std::stringstream message;
        message << "The calibration intervals lie outside of the allowed interval [0, 86400]" << std::endl;
        m_log.Error(message.str());
        returnCode = 1;
    }

    // Verify that each instrument has an initial instrument line shape file.
    for (size_t instrumentIdx = 0; instrumentIdx < m_setup.NumberOfInstruments(); ++instrumentIdx)
    {
        if (m_setup.m_instrument[instrumentIdx].m_instrumentCalibration.m_initialCalibrationFile.size() == 0)
        {
            std::stringstream message;
            message << "No initial calibration file defined for instrument " << (const char*)m_setup.m_instrument[instrumentIdx].m_serial << std::endl;
            m_log.Error(message.str());
            returnCode = 1;
        }
        if (!novac::IsExistingFile(m_setup.m_instrument[instrumentIdx].m_instrumentCalibration.m_initialCalibrationFile))
        {
            std::stringstream message;
            message << "Cannot locate the initial calibration file defined for instrument " << (const char*)m_setup.m_instrument[instrumentIdx].m_serial << std::endl;
            m_log.Error(message.str());
            returnCode = 1;
        }

        // The initial instrument line shape file is optional, but if it is defined then it must also exist.
        if (m_setup.m_instrument[instrumentIdx].m_instrumentCalibration.m_instrumentLineshapeFile.size() > 0 &&
            !novac::IsExistingFile(m_setup.m_instrument[instrumentIdx].m_instrumentCalibration.m_instrumentLineshapeFile))
        {
            std::stringstream message;
            message << "Cannot locate the initial instrument line shape file defined for instrument " << (const char*)m_setup.m_instrument[instrumentIdx].m_serial << std::endl;
            m_log.Error(message.str());
            returnCode = 1;
        }
    }
    return returnCode;
}

void CPostProcessing::CheckProcessingSettings() const
{
    novac::CString errorMessage;
    CDateTime now;

    m_log.Information("--- Validating settings --- ");

    // Get the volcano and make sure this is setup correctly.
    if (m_userSettings.m_volcano < 0)
    {
        throw std::invalid_argument("The volcano has not been setup correctly, or no volcano with the given name was found in the settings.");
    }
    {
        novac::CString volcanoName;
        g_volcanoes.GetVolcanoName(static_cast<unsigned int>(m_userSettings.m_volcano), volcanoName);
        const auto location = g_volcanoes.GetPeak(static_cast<unsigned int>(m_userSettings.m_volcano));

        std::stringstream msg;
        msg << "Monitoring volcano '" << volcanoName << " located at " << location;
        m_log.Information(msg.str());
    }

    // Check that no instrument is duplicated in the list of instruments...
    for (size_t j = 0; j < m_setup.NumberOfInstruments(); ++j)
    {
        for (size_t k = j + 1; k < m_setup.NumberOfInstruments(); ++k)
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
    for (size_t j = 0; j < m_setup.NumberOfInstruments(); ++j)
    {
        if (m_setup.m_instrument[j].m_eval.NumberOfFitWindows() == 1)
        {
            continue;
        }
        m_setup.m_instrument[j].m_eval.CheckSettings();
    }

    // Check that, for each spectrometer, there's only one location defined at each instant
    for (size_t j = 0; j < m_setup.NumberOfInstruments(); ++j)
    {
        if (m_setup.m_instrument[j].m_location.GetLocationNum() == 1)
        {
            continue;
        }
        m_setup.m_instrument[j].m_location.CheckSettings();
    }
}

void CPostProcessing::ReadWindField(novac::LogContext context)
{
    if (m_userSettings.m_volcano < 0)
    {
        throw std::invalid_argument("Volcano index has not been set when attempting to read wid field.");
    }

    m_log.Information("--- Reading Wind information --- ");

    novac::CString name1, name2, name3, path1, path2, path3, messageToUser;
    Common common;
    FileHandler::CXMLWindFileReader reader{ m_log, m_userSettings };

    if (m_userSettings.m_windFieldFileOption == 0)
    {
        // If the user has given a file-name, then try to use that one
        if (m_userSettings.m_windFieldFile.GetLength() > 3)
        {
            auto fileContext = context.With(novac::LogContext::FileName, m_userSettings.m_windFieldFile.std_str());

            m_log.Information(fileContext, "Reading wind field from file");

            reader.ReadWindFile(fileContext, m_userSettings.m_windFieldFile, m_windDataBase);

            messageToUser.Format("Parsed %d wind data items", m_windDataBase.GetDataBaseSize());
            m_log.Information(fileContext, messageToUser.std_str());

            name1.Format("%sParsedWindField.wxml", (const char*)common.m_exePath);
            m_windDataBase.WriteToFile(name1);

            return;
        }

        // Get the name of the volcano that we are about to process...
        //  there are two options, either the full name or the simple name
        const unsigned int volcanoIndex = static_cast<unsigned int>(m_userSettings.m_volcano);
        g_volcanoes.GetVolcanoName(volcanoIndex, name1);
        g_volcanoes.GetSimpleVolcanoName(volcanoIndex, name2);
        g_volcanoes.GetVolcanoCode(volcanoIndex, name3);

        // Get the path to the executable, so that we know where to start looking
        path1.Format("%sconfiguration%c%s.wxml", (const char*)common.m_exePath, Poco::Path::separator(), (const char*)name1);
        path2.Format("%sconfiguration%c%s.wxml", (const char*)common.m_exePath, Poco::Path::separator(), (const char*)name2);
        path3.Format("%sconfiguration%c%s.wxml", (const char*)common.m_exePath, Poco::Path::separator(), (const char*)name3);

        // check which of the files exists
        if (Filesystem::IsExistingFile(path1))
        {
            messageToUser.Format("Reading wind field from file: %s", (const char*)path1);
            m_log.Information(messageToUser.std_str());

            reader.ReadWindFile(context, path1, m_windDataBase);

            return;
        }

        if (Filesystem::IsExistingFile(path2))
        {
            messageToUser.Format("Reading wind field from file: %s", (const char*)path2);
            m_log.Information(messageToUser.std_str());

            reader.ReadWindFile(context, path2, m_windDataBase);

            return;
        }

        if (Filesystem::IsExistingFile(path3))
        {
            messageToUser.Format("Reading wind field from file: %s", (const char*)path3);
            m_log.Information(messageToUser.std_str());

            reader.ReadWindFile(context, path3, m_windDataBase);

            return;
        }

        messageToUser.Format("Cannot find wind field file. Attempted to read: %s, %s and %s", (const char*)path1, (const char*)path2, (const char*)path3);
        throw std::invalid_argument(messageToUser.c_str());
    }

    // If the user has specified a directory of files...
    if (m_userSettings.m_windFieldFileOption == 1)
    {
        reader.ReadWindDirectory(context, m_userSettings.m_windFieldFile, m_windDataBase, &m_userSettings.m_fromDate, &m_userSettings.m_toDate);
        return;
    }

    // should never get to this point!
    throw std::logic_error("Invalid wind field file option, failed to read winf field.");
}

void CPostProcessing::PreparePlumeHeights(novac::LogContext context)
{
    const unsigned int volcanoIndex = static_cast<unsigned int>(m_userSettings.m_volcano);

    // we need to construct a default plume height to use, if there's nothing else...
    Geometry::PlumeHeight plumeHeight;
    plumeHeight.m_plumeAltitude = g_volcanoes.GetPeakAltitude(volcanoIndex);
    plumeHeight.m_plumeAltitudeSource = Meteorology::MeteorologySource::Default;
    plumeHeight.m_validFrom = CDateTime(0, 0, 0, 0, 0, 0);
    plumeHeight.m_validTo = CDateTime(9999, 12, 31, 23, 59, 59);

    // the estimated plume height is half of the altitude difference between the highest
    // instrument for this volcano and the volcano altitude
    double maxInstrumentAltitude = -1e6;
    Configuration::CInstrumentLocation location;
    novac::CString volcanoName;
    g_volcanoes.GetVolcanoName(volcanoIndex, volcanoName);
    for (const auto& instrument : m_setup.m_instrument)
    {
        for (size_t j = 0; j < instrument.m_location.GetLocationNum(); ++j)
        {
            instrument.m_location.GetLocation(j, location);
            if (Equals(volcanoName, location.m_volcano))
            {
                maxInstrumentAltitude = std::max(maxInstrumentAltitude, double(location.m_altitude));
            }
        }
    }
    if (maxInstrumentAltitude > 0)
    {
        plumeHeight.m_plumeAltitudeError = std::abs(g_volcanoes.GetPeakAltitude(volcanoIndex) - maxInstrumentAltitude) / 2.0;
    }
    else
    {
        plumeHeight.m_plumeAltitudeError = g_volcanoes.GetPeakAltitude(volcanoIndex) / 2.0;
    }

    std::stringstream msg;
    msg << "Prepared default plume height of: " << plumeHeight.m_plumeAltitude << " +- " << plumeHeight.m_plumeAltitudeError << " masl ";
    m_log.Information(context, msg.str());

    m_plumeDataBase.InsertPlumeHeight(plumeHeight);
}

void CPostProcessing::CalculateGeometries(
    novac::LogContext context,
    std::vector<Evaluation::CExtendedScanResult>& scanResults,
    std::vector<Geometry::CGeometryResult>& geometryResults)
{
    novac::CString messageToUser;
    unsigned long nFilesChecked1 = 0; // this is for debugging purposes...
    unsigned long nFilesChecked2 = 0; // this is for debugging purposes...
    unsigned long nCalculationsMade = 0; // this is for debugging purposes...
    unsigned long nTooLongdistance = 0; // this is for debugging purposes...
    unsigned long nTooLargeAbsoluteError = 0; // this is for debugging purposes...
    unsigned long nTooLargeRelativeError = 0; // this is for debugging purposes...

    // Tell the user what's happening
    m_log.Information(context, "Begin to calculate plume heights from scans");

    // Loop through list with output text files from evaluation and apply geometrical corrections
    for (size_t pos1 = 0; pos1 < scanResults.size(); ++pos1)
    {
        const Evaluation::CExtendedScanResult& scanResult1 = scanResults[pos1];

        ++nFilesChecked1; // for debugging...

        // if this scan does not see a large enough portion of the plume, then ignore it...
        if (!scanResult1.m_scanProperties.completeness.HasValue() ||
            scanResult1.m_scanProperties.completeness.Value() < m_userSettings.m_calcGeometry_CompletenessLimit)
        {
            continue;
        }

        // If this is not a flux-measurement, then there's no use in trying to use it...
        if (scanResult1.m_measurementMode != MeasurementMode::Flux)
        {
            continue;
        }

        // try to combine this evaluation-log file with every other eval-log
        //  use the fact that the list of eval-logs is sorted by increasing start-time
        //  thus we start at the eval-log next after this one and compare with all
        //  eval-logs until the difference in start-time is too big.
        bool successfullyCombined = false; // this is true if evalLog1 was combined with (at least one) other eval-log to make a geomery calculation.
        for (size_t pos2 = pos1 + 1; pos2 < scanResults.size(); ++pos2)
        {
            const Evaluation::CExtendedScanResult& scanResult2 = scanResults[pos2];

            ++nFilesChecked2; // for debugging...

            // if this scan does not see a large enough portion of the plume, then ignore it...
            if (!scanResult2.m_scanProperties.completeness.HasValue() ||
                scanResult2.m_scanProperties.completeness.Value() < m_userSettings.m_calcGeometry_CompletenessLimit)
            {
                continue;
            }

            // The time elapsed between the two measurements must not be more than 
            // the user defined time-limit (in seconds)
            double timeDifference = std::abs(CDateTime::Difference(scanResult1.m_startTime, scanResult2.m_startTime));
            if (timeDifference > m_userSettings.m_calcGeometry_MaxTimeDifference)
            {
                break;
            }

            // If this is not a flux-measurement, then there's no use in trying to use it...
            if (scanResult2.m_measurementMode != MeasurementMode::Flux)
            {
                continue;
            }

            // the serials must be different (i.e. the two measurements must be
            //  from two different instruments)
            if (Equals(scanResult1.m_instrumentSerial, scanResult2.m_instrumentSerial))
            {
                continue;
            }

            // Get the locations of the two instruments
            Configuration::CInstrumentLocation location[2];
            try
            {
                location[0] = m_setup.GetInstrumentLocation(scanResult1.m_instrumentSerial, scanResult1.m_startTime);
                location[1] = m_setup.GetInstrumentLocation(scanResult2.m_instrumentSerial, scanResult2.m_startTime);
            }
            catch (novac::NotFoundException& ex)
            {
                m_log.Information(ex.message);
                continue;
            }

            // make sure that the distance between the instruments is not too long....
            const double instrumentDistance = novac::GpsMath::Distance(location[0].GpsData(), location[1].GpsData());
            if (instrumentDistance < m_userSettings.m_calcGeometry_MinDistance ||
                instrumentDistance > m_userSettings.m_calcGeometry_MaxDistance)
            {
                ++nTooLongdistance;
                continue;
            }

            // count the number of times we calculate a result, for improving the software...
            ++nCalculationsMade;

            // If the files have passed these tests then make a geometry-calculation
            Geometry::CGeometryResult result;
            Geometry::CGeometryCalculator geometryCalculator(m_log, m_userSettings);
            if (geometryCalculator.CalculateGeometry(scanResult1.m_scanProperties, scanResult1.m_startTime, scanResult2.m_scanProperties, scanResult2.m_startTime, location, result))
            {
                // Check the quality of the measurement before we insert it...
                if (result.m_plumeAltitudeError.Value() > m_userSettings.m_calcGeometry_MaxPlumeAltError)
                {
                    ++nTooLargeAbsoluteError; // too bad, continue.
                }
                else if ((result.m_plumeAltitudeError.Value() > 0.5 * result.m_plumeAltitude.Value()) || (result.m_windDirectionError.Value() > m_userSettings.m_calcGeometry_MaxWindDirectionError))
                {
                    ++nTooLargeRelativeError;  // too bad, continue.
                }
                else if (result.m_windDirectionError.Value() > m_userSettings.m_calcGeometry_MaxWindDirectionError)
                {
                    // ... too bad, continue.
                }
                else
                {
                    // remember which instruments were used
                    result.m_instrumentSerial1 = scanResult1.m_instrumentSerial;
                    result.m_instrumentSerial2 = scanResult2.m_instrumentSerial;

                    geometryResults.push_back(result);

                    messageToUser.Format("Calculated a plume altitude of %.0lf +- %.0lf masl and wind direction of %.0lf +- %.0lf degrees by combining measurements two instruments",
                        result.m_plumeAltitude.Value(),
                        result.m_plumeAltitudeError.Value(),
                        result.m_windDirection.Value(),
                        result.m_windDirectionError.Value());
                    m_log.Information(context.With("device1", scanResult1.m_instrumentSerial).With("device2", scanResult2.m_instrumentSerial).WithTimestamp(scanResult1.m_startTime), messageToUser.std_str());

                    successfullyCombined = true;
                }
            }
        } // end for (size_t pos2 = pos1 + 1;

        // if it was not possible to combine this scan with any other to generate an
        // estimated plume height and wind direction we might still be able to use it to calculate
        // a wind direction given the plume height at the time of the measurement.
        if (!successfullyCombined)
        {
            Meteorology::WindField windField;

            // Get the location of the instrument
            Configuration::CInstrumentLocation location;
            try
            {
                location = m_setup.GetInstrumentLocation(scanResult1.m_instrumentSerial, scanResult1.m_startTime);
            }
            catch (novac::NotFoundException& ex)
            {
                m_log.Information(context, ex.message);
                continue;
            }

            // Get the altitude of the plume at this moment. First look into the
            // general database. Then have a look in the list of geometry-results
            // that we just generated to see if there's anything better there...
            Geometry::PlumeHeight plumeHeight;
            m_plumeDataBase.GetPlumeHeight(scanResult1.m_startTime, plumeHeight);

            for (auto it = geometryResults.rbegin(); it != geometryResults.rend(); ++it)
            {
                const Geometry::CGeometryResult& oldResult = *it;
                if (std::abs(CDateTime::Difference(oldResult.m_averageStartTime, scanResult1.m_startTime)) < m_userSettings.m_calcGeometryValidTime)
                {
                    if (oldResult.m_plumeAltitude.HasValue() &&
                        (!plumeHeight.m_plumeAltitude || oldResult.m_plumeAltitudeError.Value() < plumeHeight.m_plumeAltitudeError))
                    {
                        plumeHeight.m_plumeAltitude = oldResult.m_plumeAltitude.Value();
                        plumeHeight.m_plumeAltitudeError = oldResult.m_plumeAltitudeError.Value();
                        plumeHeight.m_plumeAltitudeSource = oldResult.m_calculationType;
                    }
                }
            }

            // Try to calculate the wind-direction
            Geometry::CGeometryResult result;
            Geometry::CGeometryCalculator geometryCalculator(m_log, m_userSettings);
            if (geometryCalculator.CalculateWindDirection(scanResult1.m_scanProperties, scanResult1.m_startTime, plumeHeight, location, result))
            {
                if (result.m_windDirectionError.Value() > m_userSettings.m_calcGeometry_MaxWindDirectionError)
                {
                    continue;
                }

                // Success!!
                result.m_instrumentSerial1 = scanResult1.m_instrumentSerial;
                geometryResults.push_back(result);

                // tell the user   
                messageToUser.Format("Calculated a wind direction of %.0lf +- %.0lf degrees from a scan",
                    result.m_windDirection, result.m_windDirectionError);
                m_log.Information(context.With(novac::LogContext::Device, scanResult1.m_instrumentSerial).WithTimestamp(result.m_averageStartTime), messageToUser.std_str());
            }
            else
            {
                continue;
            }
        }
    } // end for (size_t pos1 = 0; pos1 < scanResults.size(); ++pos1)

    // Tell the user what we have done
    if (geometryResults.size() == 0)
    {
        m_log.Information(context, "No plume heights could be calculated");
    }
    else
    {
        std::vector<double> plumeHeights;
        std::vector<double> windDirections;
        plumeHeights.reserve(geometryResults.size());
        windDirections.reserve(geometryResults.size());
        for (const auto& g : geometryResults)
        {
            if (g.m_plumeAltitude.HasValue())
            {
                plumeHeights.push_back(g.m_plumeAltitude.Value());
            }
            if (g.m_windDirection.HasValue())
            {
                windDirections.push_back(g.m_windDirection.Value());
            }
        }

        std::stringstream msg;
        msg << "Done calculating geometries. Plume height calculated on " << plumeHeights.size() << " occasions. ";
        msg << "Median plume height: " << Median(plumeHeights) << ", min plume height: " << Min(plumeHeights) << ", max plume height: " << Max(plumeHeights) << " [masl]" << std::endl;
        msg << "Wind direction calculated on " << windDirections.size() << " occasions. ";
        msg << "Wind direction in range: [" << Min(windDirections) << "; " << Max(windDirections) << "] [deg]";
        m_log.Information(context, msg.str());
    }
    messageToUser.Format("nFilesChecked1 = %ld, nFilesChecked2 = %ld, nCalculationsMade = %ld", nFilesChecked1, nFilesChecked2, nCalculationsMade);
    m_log.Information(context, messageToUser.std_str());
}

void CPostProcessing::CalculateFluxes(novac::LogContext context, const std::vector<Evaluation::CExtendedScanResult>& scanResults)
{
    Flux::CFluxStatistics stat;

    // we keep the calculated fluxes in a list
    std::list<Flux::FluxResult> calculatedFluxes;

    // Initiate the flux-calculator
    Flux::CFluxCalculator fluxCalc(m_log, m_setup, m_userSettings);

    // Loop through the list of evaluation log files. For each of them, find
    // the best available wind-speed, wind-direction and plume height and
    // calculate the flux.
    for (const auto& scanResult : scanResults)
    {
        // Get the name of this eval-log
        const novac::CString& evalLog = scanResult.m_evalLogFile[m_userSettings.m_mainFitWindow];
        const CPlumeInScanProperty& plume = scanResult.m_scanProperties;

        const novac::LogContext fileContext = context.With(novac::LogContext::FileName, novac::GetFileName(evalLog.std_str()));

        // If this is not a flux-measurement, then there's no point in calculating any flux for it
        if (scanResult.m_measurementMode != MeasurementMode::Flux)
        {
            continue;
        }
        if (!plume.completeness.HasValue())
        {
            m_log.Information(fileContext, "Scan does not see the plume. Will not calculate any flux.");
            continue;
        }

        // if the completeness is too low then ignore this scan.
        if (plume.completeness.Value() < (m_userSettings.m_completenessLimitFlux + 0.01))
        {
            std::stringstream msg;
            msg << "Scan has completeness = " << plume.completeness.Value() << " which is less than limit of " << m_userSettings.m_completenessLimitFlux << ". Will not calculate any flux.";
            m_log.Information(fileContext, msg.str());
            continue;
        }

        // Extract a plume height at this time of day
        Geometry::PlumeHeight plumeHeight;
        if (!m_plumeDataBase.GetPlumeHeight(scanResult.m_startTime, plumeHeight))
        {
            std::stringstream msg;
            msg << "Failed to get plume height at the time of the measurement (" << scanResult.m_startTime << ") no flux calculated for scan.";
            m_log.Information(fileContext, msg.str());
            continue;
        }

        // Calculate the flux. This also takes care of writing
        // the results to file
        Flux::FluxResult fluxResult;
        if (fluxCalc.CalculateFlux(fileContext, scanResult, m_windDataBase, plumeHeight, fluxResult))
        {
            std::stringstream msg;
            msg << "Calculated flux of " << fluxResult.m_flux << " [kg/s] for scan.";
            m_log.Information(fileContext, msg.str());
            calculatedFluxes.push_back(fluxResult);
        }
        else
        {
            m_log.Information(fileContext, "No flux calculated for scan.");
        }
    }

    // Now we can write the final fluxes to file
    m_log.Information(context, "Writing flux log");
    WriteFluxResult_XML(calculatedFluxes);
    WriteFluxResult_Txt(calculatedFluxes);

    // Also write the statistics for the flux
    m_log.Information(context, "Writing flux statistics");
    stat.AttachFluxList(calculatedFluxes);

    novac::CString fluxStatFileName;
    fluxStatFileName.Format("%s%c%s", m_userSettings.m_outputDirectory.c_str(), Poco::Path::separator(), "FluxStatistics.txt");
    stat.WriteFluxStat(fluxStatFileName);
}

void CPostProcessing::WriteFluxResult_XML(const std::list<Flux::FluxResult>& calculatedFluxes)
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
        m_log.Information("Could not open flux log file for writing. Writing of results failed. ");
        return;
    }

    // Write the header and the starting comments
    fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(f, "<?xml-stylesheet type=\"text/xsl\" href=\"fluxresult.xsl\"?>\n");
    fprintf(f, "<!-- This is result of the flux calculations the NOVAC Post Processing Program -->\n");
    fprintf(f, "<!-- File generated on %04d.%02d.%02d at %02d:%02d:%02d -->\n\n", now.year, now.month, now.day, now.hour, now.minute, now.second);

    fprintf(f, "<NovacPPPFluxResults>\n");

    for (const Flux::FluxResult& fluxResult : calculatedFluxes)
    {
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

        fprintf(f, "\t\t<serial>%s</serial>\n", fluxResult.m_instrument.c_str());

        // extract the instrument type
        if (fluxResult.m_instrumentType == NovacInstrumentType::Heidelberg)
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
        if (fluxResult.m_fluxQualityFlag == FluxQuality::Green)
        {
            fprintf(f, "\t\t<Quality>g</Quality>\n");
        }
        else if (fluxResult.m_fluxQualityFlag == FluxQuality::Yellow)
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

void CPostProcessing::WriteFluxResult_Txt(const std::list<Flux::FluxResult>& calculatedFluxes)
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
        m_log.Information("Could not open flux log file for writing. Writing of results failed. ");
        return;
    }

    // Write the header and the starting comments
    fprintf(f, "# This is result of the flux calculations the NOVAC Post Processing Program \n");
    fprintf(f, "#   File generated on %04d.%02d.%02d at %02d:%02d:%02d \n\n", now.year, now.month, now.day, now.hour, now.minute, now.second);

    fprintf(f, "#StartTime\tStopTime\tSerial\tInstrumentType\tFlux_kgs\tFluxQuality\tFluxError_Wind_kgs\tFluxError_PlumeHeight_kgs\tWindSpeed_ms\tWindSpeedErr_ms\tWindSpeedSrc\tWindDir_deg\tWindDirErr_deg\tWindDirSrc\tPlumeHeight_m\tPlumeHeightErr_m\tPlumeHeightSrc\t");
    fprintf(f, "Compass\tConeAngle\tTilt\tnSpectra\tPlumeCentre_1\tPlumeCentre_2\tPlumeCompleteness\tScanOffset\n");

    for (const Flux::FluxResult& fluxResult : calculatedFluxes)
    {
        // extract the instrument type
        if (fluxResult.m_instrumentType == NovacInstrumentType::Heidelberg)
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
        fprintf(f, "%s\t", fluxResult.m_instrument.c_str());
        fprintf(f, "%s\t", (const char*)typeStr);

        // The actual flux!!!
        fprintf(f, "%.2lf\t", fluxResult.m_flux);

        // The judged quality of the calculated flux
        if (fluxResult.m_fluxQualityFlag == FluxQuality::Green)
        {
            fprintf(f, "g\t");
        }
        else if (fluxResult.m_fluxQualityFlag == FluxQuality::Yellow)
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

void CPostProcessing::WriteCalculatedGeometriesToFile(novac::LogContext context, const std::vector<Geometry::CGeometryResult>& geometryResults)
{
    if (geometryResults.size() == 0)
    {
        return; // nothing to write...
    }

    FILE* f = nullptr;
    const std::string geomLogFile = m_userSettings.m_outputDirectory.std_str() + "/" + "GeometryLog.csv";

    if (Filesystem::IsExistingFile(geomLogFile))
    {
        f = fopen(geomLogFile.c_str(), "a");
        if (f == nullptr)
        {
            m_log.Information(context, "Could not open geometry log file for writing. Writing of results failed. ");
            return;
        }
    }
    else
    {
        f = fopen(geomLogFile.c_str(), "w");
        if (f == nullptr)
        {
            m_log.Information(context, "Could not open geometry log file for writing. Writing of results failed. ");
            return;
        }
        fprintf(f, "Date;Time;DifferenceInStartTime_minutes;Instrument1;Instrument2;PlumeAltitude_masl;PlumeHeightError_m;WindDirection_deg;WindDirectionError_deg;PlumeCentre1_deg;PlumeCentreError1_deg;PlumeCentre2_deg;PlumeCentreError2_deg\n");
    }

    for (const auto& result : geometryResults)
    {
        // write the file
        if (result.m_calculationType == Meteorology::MeteorologySource::GeometryCalculationTwoInstruments)
        {
            fprintf(f, "%04d.%02d.%02d;", result.m_averageStartTime.year, result.m_averageStartTime.month, result.m_averageStartTime.day);
            fprintf(f, "%02d:%02d:%02d;", result.m_averageStartTime.hour, result.m_averageStartTime.minute, result.m_averageStartTime.second);
            fprintf(f, "%.1lf;", result.m_startTimeDifference / 60.0);
            fprintf(f, "%s;%s;", result.m_instrumentSerial1.c_str(), result.m_instrumentSerial2.c_str());
            fprintf(f, "%.0lf;%.0lf;", result.m_plumeAltitude.Value(), result.m_plumeAltitudeError.Value());
            fprintf(f, "%.0lf;%.0lf;", result.m_windDirection.Value(), result.m_windDirectionError.Value());

            fprintf(f, "%.1f;%.1f;", result.m_plumeCentre1.Value(), result.m_plumeCentreError1.Value());
            fprintf(f, "%.1f;%.1f\n", result.m_plumeCentre2.Value(), result.m_plumeCentreError2.Value());
        }
        else
        {
            fprintf(f, "%04d.%02d.%02d;", result.m_averageStartTime.year, result.m_averageStartTime.month, result.m_averageStartTime.day);
            fprintf(f, "%02d:%02d:%02d;", result.m_averageStartTime.hour, result.m_averageStartTime.minute, result.m_averageStartTime.second);
            fprintf(f, "0;");
            fprintf(f, "%s;;", result.m_instrumentSerial1.c_str());
            fprintf(f, "0;0;");
            fprintf(f, "%.0lf;%.0lf;", result.m_windDirection.Value(), result.m_windDirectionError.Value());

            fprintf(f, "%.1f;%.1f;", result.m_plumeCentre1.Value(), result.m_plumeCentreError1.Value());
            fprintf(f, ";\n");
        }
    }
    fclose(f);
}

void CPostProcessing::InsertCalculatedGeometriesIntoDatabase(novac::LogContext context, const std::vector<Geometry::CGeometryResult>& geometryResults)
{
    for (const auto& result : geometryResults)
    {
        if (result.m_plumeAltitude.HasValue())
        {
            this->m_plumeDataBase.InsertPlumeHeight(result);
        }

        if (result.m_windDirection.HasValue())
        {
            try
            {
                // get the location of the instrument at the time of the measurement
                // TODO: This was retrieved but not used further down. Should it be used?
                // const auto location = m_setup.GetInstrumentLocation(result.m_instrumentSerial1, result.m_averageStartTime);

                // get the time-interval that the measurement is valid for
                CDateTime validFrom = result.m_averageStartTime;
                validFrom.Decrement(m_userSettings.m_calcGeometryValidTime);
                CDateTime validTo = result.m_averageStartTime;
                validTo.Increment(m_userSettings.m_calcGeometryValidTime);

                // insert the wind-direction into the wind database
                m_windDataBase.InsertWindDirection(validFrom, validTo, result.m_windDirection.Value(), result.m_windDirectionError.Value(), result.m_calculationType, nullptr);
            }
            catch (novac::NotFoundException& ex)
            {
                m_log.Information(context, ex.message);
            }
        }
    }
}

void CPostProcessing::CalculateDualBeamWindSpeeds(novac::LogContext context, const std::vector<Evaluation::CExtendedScanResult>& evalLogs)
{
    std::vector<std::string> masterList; // list of wind-measurements from the master channel
    std::vector<std::string> slaveList;  // list of wind-measurements from the slave channel
    std::vector<std::string> heidelbergList;  // list of wind-measurements from the Heidelbergensis

    CDateTime validFrom, validTo;

    novac::CString serial, serial2, fileName, fileName2, nonsenseString;
    novac::CString userMessage, windLogFile;
    CDateTime startTime, startTime2;
    int channel, channel2, nWindMeasFound = 0;
    MeasurementMode meas_mode, meas_mode2;
    Configuration::CInstrumentLocation location;
    WindSpeedMeasurement::CWindSpeedCalculator calculator(m_log, m_userSettings);
    Geometry::PlumeHeight plumeHeight;
    Meteorology::WindField windField, oldWindField;

    // -------------------------------- step 1. -------------------------------------
    // search through 'evalLogs' for dual-beam measurements from master and from slave
    for (const auto& scanResult : evalLogs)
    {
        const novac::CString& fileNameAndPath = scanResult.m_evalLogFile[m_userSettings.m_mainFitWindow];

        try
        {
            // to know the start-time of the measurement, we need to 
            // extract just the file-name, i.e. remove the path
            fileName = novac::CString(fileNameAndPath);
            Common::GetFileName(fileName);

            novac::CFileUtils::GetInfoFromFileName(fileName, startTime, serial, channel, meas_mode);

            if (meas_mode == MeasurementMode::Windspeed)
            {
                ++nWindMeasFound;

                // first check if this is a heidelberg instrument
                location = m_setup.GetInstrumentLocation(serial.std_str(), startTime);

                if (location.m_instrumentType == NovacInstrumentType::Heidelberg)
                {
                    // this is a heidelberg instrument
                    heidelbergList.push_back(fileNameAndPath.std_str());
                }
                else
                {
                    // this is a gothenburg instrument
                    if (channel == 0)
                    {
                        masterList.push_back(fileNameAndPath.std_str());
                    }
                    else if (channel == 1)
                    {
                        slaveList.push_back(fileNameAndPath.std_str());
                    }
                }
            }
        }
        catch (novac::NotFoundException& ex)
        {
            novac::LogContext fileContext = context.With(novac::LogContext::FileName, fileNameAndPath.std_str());
            m_log.Information(fileContext, ex.message);
        }
    }
    if (nWindMeasFound == 0)
    {
        m_log.Information(context, "No dual-beam wind speed measurements found.");
        return; // if nothing was found...
    }

    userMessage.Format("%d dual-beam wind speed measurements found. Calculating wind-speeds", nWindMeasFound);
    m_log.Information(context, userMessage.std_str());

    // Create the dual-beam log-file
    windLogFile.Format("%s%cDualBeamLog.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    calculator.WriteWindSpeedLogHeader(windLogFile);


    // -------------------------------- step 2. -------------------------------------
    // loop through each of the measurements from the heidelberg instruments
    // and calculate the wind speed for each measurement
    for (const auto& fileNameAndPath : heidelbergList)
    {
        fileName = novac::CString(fileNameAndPath);
        Common::GetFileName(fileName);

        novac::LogContext fileContext = context.With(novac::LogContext::FileName, fileName.std_str());
        fileContext = fileContext.With("mode", "dualBeamWindSpeed");

        try
        {
            // to know the start-time of the measurement, we need to 
            // extract just the file-name, i.e. remove the path
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
                    userMessage.Format("Calculated a wind-speed of %.1lf +- %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Error too large, measurement discarded.", windField.GetWindSpeed(), windField.GetWindSpeedError(),
                        startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);
                }
                else
                {
                    // tell the user...
                    userMessage.Format("Calculated a wind-speed of %.1lf +- %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Measurement accepted", windField.GetWindSpeed(), windField.GetWindSpeedError(),
                        startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);

                    // get the time-interval that the measurement is valid for
                    windField.GetValidTimeFrame(validFrom, validTo);

                    // insert the new wind speed into the database
                    m_windDataBase.InsertWindSpeed(validFrom, validTo, windField.GetWindSpeed(), windField.GetWindSpeedError(), Meteorology::MeteorologySource::DualBeamMeasurement, nullptr);
                }
                m_log.Information(fileContext, userMessage.std_str());
            }
            else
            {
                m_log.Information(fileContext, "Failed to calculate wind speed from measurement.");
            }
        }
        catch (novac::NotFoundException& ex)
        {
            m_log.Information(fileContext, ex.message);
        }
    }

    // -------------------------------- step 3. -------------------------------------
    // loop through each of the measurements from a master-channel and try to match them with a measurement
    // from a slave channel...
    for (const auto& fileNameAndPath : masterList)
    {
        // extract just the file-name, i.e. remove the path
        fileName = novac::CString(fileNameAndPath);
        Common::GetFileName(fileName);

        novac::LogContext fileContext = context.With("file1", fileName.std_str());
        fileContext = fileContext.With("mode", "dualBeamWindSpeed");

        novac::CFileUtils::GetInfoFromFileName(fileName, startTime, serial, channel, meas_mode);

        // now check if we can match this one with a file in the slave-channel
        for (const auto& fileNameAndPath2 : slaveList)
        {
            // extract just the file-name, i.e. remove the path
            fileName2 = novac::CString(fileNameAndPath2);
            Common::GetFileName(fileName2);

            novac::CFileUtils::GetInfoFromFileName(fileName2, startTime2, serial2, channel2, meas_mode2);

            if (Equals(serial, serial2) && (startTime == startTime2))
            {
                // we have found a match!!!
                novac::LogContext calcContext = fileContext.With("file2", fileName2.std_str());

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
                        m_windDataBase.InsertWindSpeed(validFrom, validTo, windField.GetWindSpeed(), windField.GetWindSpeedError(), Meteorology::MeteorologySource::DualBeamMeasurement, nullptr);
                    }
                    m_log.Information(calcContext, userMessage.std_str());

                }
                else
                {
                    m_log.Information(calcContext, "Failed to calculate wind speed from measurement.");
                }
            }
        }
    }
}

void CPostProcessing::SortEvaluationLogs(std::vector<Evaluation::CExtendedScanResult>& evalLogs)
{
    std::vector<Evaluation::CExtendedScanResult> left;
    std::vector<Evaluation::CExtendedScanResult> right;

    // If this list consists of only one element, then we're done
    if (evalLogs.size() <= 1)
    {
        return;
    }

    // Divide the list into two, and sort each one of them
    int index = 0;
    for (auto& log : evalLogs)
    {
        if (index % 2 == 0)
            left.push_back(log);
        else
            right.push_back(log);
        ++index;
    }

    SortEvaluationLogs(left);
    SortEvaluationLogs(right);

    // Merge the two lists into one, sorted list
    evalLogs.clear();
    auto iteratorLeft = left.begin();
    auto iteratorRight = right.begin();
    while (iteratorLeft != left.end() && iteratorRight != right.end())
    {
        Evaluation::CExtendedScanResult& leftItem = *iteratorLeft;
        Evaluation::CExtendedScanResult& rightItem = *iteratorRight;

        // While we're looking through the data, do some sanity check of the data.
        assert(!leftItem.m_instrumentSerial.empty());
        assert(!rightItem.m_instrumentSerial.empty());
        assert(!leftItem.m_evalLogFile.empty());
        assert(!rightItem.m_evalLogFile.empty());

        if (rightItem.m_startTime < leftItem.m_startTime)
        {
            evalLogs.push_back(rightItem);
            ++iteratorRight;
        }
        else
        {
            evalLogs.push_back(leftItem);
            ++iteratorLeft;
        }
    }

    while (iteratorLeft != left.end())
    {
        evalLogs.push_back(*iteratorLeft);
        ++iteratorLeft;
    }

    while (iteratorRight != right.end())
    {
        evalLogs.push_back(*iteratorRight);
        ++iteratorRight;
    }

    return;
}

void CPostProcessing::UploadResultsToFTP(novac::LogContext context)
{
    Communication::CFTPServerConnection connection(m_log, m_userSettings);
    novac::CString fileName;

    // Generate a list with all the files we want to upload.
    std::vector<std::string> fileList;

    // 1. the geometry log file
    fileName.Format("%s%cGeometryLog.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    fileList.push_back(fileName.std_str());

    // 2. the generated wind field
    fileName.Format("%s%cGeneratedWindField.wxml", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    fileList.push_back(fileName.std_str());

    // 3. the generated flux log
    fileName.Format("%s%cFluxLog.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());
    fileList.push_back(fileName.std_str());

    // upload the files
    connection.UploadResults(context, m_userSettings.m_FTPDirectory, m_userSettings.m_FTPUsername, m_userSettings.m_FTPPassword, fileList);

    return;
}

std::vector<Evaluation::CExtendedScanResult> CPostProcessing::LocateEvaluationLogFiles(novac::LogContext context, const std::string& directory) const
{
    std::vector<std::string> filenames;
    std::vector<Evaluation::CExtendedScanResult> evaluationLogFiles;

    context = context.With(novac::LogContext::Directory, directory);
    m_log.Information(context, "Searching for evaluation log files in directory.");

    const bool includeSubDirs = (m_userSettings.m_includeSubDirectories_Local > 0);
    Filesystem::FileSearchCriterion limits;
    limits.startTime = m_userSettings.m_fromDate;
    limits.endTime = m_userSettings.m_toDate;
    limits.fileExtension = "_flux.txt";
    Filesystem::SearchDirectoryForFiles(directory, includeSubDirs, filenames, &limits);


    novac::CString messageToUser;
    messageToUser.Format("%d Evaluation log files found, starting reading", filenames.size());
    m_log.Information(context, messageToUser.std_str());

    size_t nofFailedLogReads = 0;

    for (const std::string& filename : filenames)
    {
        int channel;
        CDateTime startTime;
        MeasurementMode mode;
        novac::CString serial;
        novac::CFileUtils::GetInfoFromFileName(filename, startTime, serial, channel, mode);

        novac::LogContext filenameContext = context.With(LogContext::FileName, novac::GetFileName(filename));

        FileHandler::CEvaluationLogFileHandler logReader(m_log, filename, m_userSettings.m_molecule);
        if (RETURN_CODE::SUCCESS != logReader.ReadEvaluationLog() || logReader.m_scan.size() == 0)
        {
            ++nofFailedLogReads;
            continue;
        }
        if (logReader.m_scan.size() > 1)
        {
            m_log.Information(filenameContext, "File contained mored than one scan. Only the first will be used.");
        }

        Evaluation::CScanResult scanResult = logReader.m_scan[0];

        std::string message;
        auto plumeProperties = novac::CalculatePlumeProperties(scanResult, m_userSettings.m_molecule, message);
        if (plumeProperties == nullptr)
        {
            m_log.Information(filenameContext, message + " Scan does not see the plume, no flux will be calculated.");
            continue;
        }
        else
        {
            scanResult.m_plumeProperties = *plumeProperties;

            std::stringstream msg;
            msg << "Scan sees the plume at: " << scanResult.m_plumeProperties.plumeCenter;
            msg << " +- " << scanResult.m_plumeProperties.plumeCenterError;
            msg << " [deg]. Completeness:" << scanResult.m_plumeProperties.completeness;
            m_log.Information(filenameContext, msg.str());
        }

        Evaluation::CExtendedScanResult result(scanResult.GetSerial(), startTime, mode);
        result.m_pakFile = ""; // unknown
        result.m_evalLogFile.push_back(filename);
        result.m_fitWindowName.push_back(""); // unknown
        result.m_startTime = startTime;
        result.m_scanProperties = scanResult.m_plumeProperties;

        evaluationLogFiles.push_back(result);
    }

    messageToUser.Format("%d Evaluation log files read successfully. %d scans see the plume.", filenames.size() - nofFailedLogReads, evaluationLogFiles.size());
    m_log.Information(context, messageToUser.std_str());

    return evaluationLogFiles;
}
