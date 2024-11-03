#include "stdafx.h"
#include <SpectralEvaluation/StringUtils.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/MFC/CStringTokenizer.h>
#include <PPPLib/VolcanoInfo.h>
#include <PPPLib/File/Filesystem.h>

#include <PPPLib/Configuration/CommandLineParser.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/File/SetupFileReader.h>
#include <PPPLib/File/EvaluationConfigurationParser.h>
#include <PPPLib/File/ProcessingFileReader.h>
#include "PostProcessing.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <thread>
#include <Poco/Path.h>
#include <Poco/Logger.h>
#include <Poco/FileChannel.h>
#include <Poco/SplitterChannel.h>
#include <Poco/ConsoleChannel.h>
#include <Poco/PatternFormatter.h>
#include <Poco/FormattingChannel.h>
#include <Poco/Util/Application.h>
#include "Common/Common.h"

extern Configuration::CUserConfiguration            g_userSettings;// <-- The settings of the user

novac::CVolcanoInfo g_volcanoes;   // <-- A list of all known volcanoes
PocoLogger g_logger; // <-- global logger

static Configuration::CNovacPPPConfiguration s_setup; // <-- The settings

std::string s_exePath;
std::string s_exeFileName;

#undef min
#undef max

void LoadConfigurations(
    const novac::CString& workDir,
    Configuration::CNovacPPPConfiguration& configuration,
    Configuration::CUserConfiguration& userSettings);

void ReadEvaluationXmlFile(
    const novac::CString& workDir,
    FileHandler::CEvaluationConfigurationParser& eval_reader,
    Configuration::CInstrumentConfiguration& instrument);

void ReadProcessingXml(const novac::CString& workDir, Configuration::CUserConfiguration& userSettings);

void ReadSetupXml(const novac::CString& workDir, Configuration::CNovacPPPConfiguration& configuration);

void StartProcessing();
void CalculateAllFluxes(CContinuationOfProcessing continuation);

using namespace novac;

class NovacPPPApplication : public Poco::Util::Application
{
public:
    NovacPPPApplication() : Poco::Util::Application(){
        this->setUnixOptions(false);
    }

protected:
    void initialize(Poco::Util::Application& application)
    {
        this->loadConfiguration();
        Poco::Util::Application::initialize(application);
    }

    void uninitialize()
    {
        Poco::Util::Application::uninitialize();
    }

    void defineOptions(Poco::Util::OptionSet& optionSet)
    {
        Poco::Util::Application::defineOptions(optionSet);
        // 
        // optionSet.addOption(
        //         Poco::Util::Option("WorkDir", "W", "The working directory")
        //                 .required(false)
        //                 .repeatable(true)
        //                 .argument("file")
        //                 .callback(Poco::Util::OptionCallback<NovacPPPApplication>(this, &NovacPPPApplication::handleMyOpt))
        // );
    }

    void handleMyOpt(const std::string& name, const std::string& value)
    {
        std::cout << "Setting option " << name << " to " << value << std::endl;
        this->config().setString(name, value);
        std::cout << "The option is now " << this->config().getString(name) << std::endl;
    }

    int main(const std::vector<std::string>& arguments)
    {
        std::cout << "Novac Post Processing Program" << std::endl;

        try
        {
            auto& application = Poco::Util::Application::instance();
            Poco::Path executable(application.commandPath());

            s_exePath = executable.parent().toString();
            s_exeFileName = executable.getFileName();

            // Setup the logging
            Poco::AutoPtr<Poco::PatternFormatter> patternFormatter(new Poco::PatternFormatter);
            patternFormatter->setProperty("pattern", "%Y-%m-%d %H:%M:%S: %t");

            Poco::AutoPtr<Poco::SplitterChannel> splitterChannel(new Poco::SplitterChannel());
            splitterChannel->addChannel(new Poco::ConsoleChannel());

            Poco::AutoPtr<Poco::FormattingChannel> formattingChannel(new Poco::FormattingChannel(patternFormatter, splitterChannel));

            Poco::Logger::root().setChannel(formattingChannel);
            Poco::Logger& log = Poco::Logger::get("NovacPPP");

            // Get the options from the command line
            ShowMessage("Getting command line arguments");
            Configuration::CommandLineParser::ParseCommandLineOptions(arguments, g_userSettings, g_volcanoes, s_exePath, g_logger);
            ShowMessage(novac::CString::FormatString(" Executing %s in '%s'", s_exeFileName.c_str(), s_exePath.c_str()));

            // Read the configuration files
            ShowMessage("Loading configuration");
            LoadConfigurations(s_exePath, s_setup, g_userSettings);

            // Read the command line options again, in order to make sure the command line arguments override the configuration files.
            Configuration::CommandLineParser::ParseCommandLineOptions(arguments, g_userSettings, g_volcanoes, s_exePath, g_logger);

            splitterChannel->addChannel(new Poco::FileChannel(g_userSettings.m_outputDirectory.std_str() + "StatusLog.txt"));
            log.setChannel(formattingChannel);

            // Start calculating the fluxes, this is the old button handler
            ShowMessage("Setup done: starting calculations");
            StartProcessing();
        }
        catch (Poco::FileNotFoundException& e)
        {
            ShowMessage("FileNotFoundException: " + e.displayText());
            return 1;
        }
        catch (std::invalid_argument& e)
        {
            std::stringstream msg;
            msg << "Invalid argument exception caught: " << e.what();
            ShowMessage(msg.str());
            return 1;
        }
        catch (std::exception& e)
        {
            std::stringstream msg;
            msg << "General exception caught: " << e.what();
            ShowMessage(msg.str());
            return 1;
        }

        return 0;
    }
};

POCO_APP_MAIN(NovacPPPApplication)


void LoadConfigurations(
    const novac::CString& workDir,
    Configuration::CNovacPPPConfiguration& configuration,
    Configuration::CUserConfiguration& userSettings)
{
    ReadSetupXml(workDir, configuration);

    ReadProcessingXml(workDir, userSettings);

    // Check if there is a configuration file for every spectrometer serial number
    FileHandler::CEvaluationConfigurationParser eval_reader{ g_logger };
    for (size_t k = 0; k < configuration.NumberOfInstruments(); ++k)
    {
        ReadEvaluationXmlFile(workDir, eval_reader, configuration.m_instrument[k]);
    }
}

void ReadEvaluationXmlFile(
    const novac::CString& workDir,
    FileHandler::CEvaluationConfigurationParser& eval_reader,
    Configuration::CInstrumentConfiguration& instrument)
{
    novac::CString evalConfPath;
    evalConfPath.Format("%sconfiguration%c%s.exml", (const char*)workDir, Poco::Path::separator(), (const char*)instrument.m_serial);

    if (Filesystem::IsExistingFile(evalConfPath))
    {
        eval_reader.ReadConfigurationFile(
            evalConfPath,
            instrument.m_eval,
            instrument.m_darkCurrentCorrection,
            instrument.m_instrumentCalibration);
    }
    else
    {
        throw std::logic_error("Could not find configuration file: " + evalConfPath);
    }
}

void ReadProcessingXml(const novac::CString& workDir, Configuration::CUserConfiguration& userSettings)
{
    novac::CString processingPath;
    processingPath.Format("%sconfiguration%cprocessing.xml", (const char*)workDir, Poco::Path::separator());
    FileHandler::CProcessingFileReader processing_reader{ g_logger };
    processing_reader.ReadProcessingFile(processingPath, userSettings);
}

void ReadSetupXml(const novac::CString& workDir, Configuration::CNovacPPPConfiguration& configuration)
{
    novac::CString setupPath;
    setupPath.Format("%sconfiguration%csetup.xml", (const char*)workDir, Poco::Path::separator());

    FileHandler::CSetupFileReader reader{ g_logger };
    reader.ReadSetupFile(setupPath, configuration);

    ShowMessage(novac::CString::FormatString(" Parsed %s, %d instruments found.", setupPath.c_str(), configuration.NumberOfInstruments()));
}

void StartProcessing()
{
    // Make sure that the ftp-path ends with a '/'
    if (g_userSettings.m_FTPDirectory.GetLength() > 1)
    {
        if (!Equals(g_userSettings.m_FTPDirectory.Right(1), "/"))
        {
            g_userSettings.m_FTPDirectory.Append("/");
        }
    }

    CContinuationOfProcessing continuation(g_userSettings);

    // Run
#ifdef _MFC_VER 
    CWinThread* postProcessingthread = AfxBeginThread(CalculateAllFluxes, NULL, THREAD_PRIORITY_NORMAL, 0, 0, NULL);
    Common::SetThreadName(postProcessingthread->m_nThreadID, "PostProcessing");
#else
    std::thread postProcessingThread(CalculateAllFluxes, continuation);
    postProcessingThread.join();
#endif  // _MFC_VER 
}

void ArchiveSettingsFiles(const Configuration::CUserConfiguration& userSettings)
{
    // set the directory to which we want to copy the settings
    novac::CString confCopyDir;
    confCopyDir.Format("%scopiedConfiguration", (const char*)userSettings.m_outputDirectory);
    confCopyDir = Filesystem::AppendPathSeparator(confCopyDir);

    // make sure that the output directory exists
    if (Filesystem::CreateDirectoryStructure(userSettings.m_outputDirectory))
    {
        novac::CString userMessage;
        userMessage.Format("Could not create output directory: %s", (const char*)userSettings.m_outputDirectory);
        throw novac::FileIoException(userMessage.c_str());
    }

    if (Filesystem::CreateDirectoryStructure(confCopyDir))
    {
        novac::CString userMessage;
        userMessage.Format("Could not create directory for copied configuration: %s", (const char*)confCopyDir);
        throw novac::FileIoException(userMessage.c_str());
    }
    // we want to copy the setup and processing files to the confCopyDir
    novac::CString processingOutputFile, setupOutputFile;
    processingOutputFile.Format("%sprocessing.xml", (const char*)confCopyDir);
    setupOutputFile.Format("%ssetup.xml", (const char*)confCopyDir);

    Common::ArchiveFile(setupOutputFile);
    Common::ArchiveFile(processingOutputFile);

    FileHandler::CProcessingFileReader writer{ g_logger };
    writer.WriteProcessingFile(processingOutputFile, userSettings);

    Common common;
    Common::CopyFile(common.m_exePath + "configuration/setup.xml", setupOutputFile);
    for (size_t k = 0; k < s_setup.NumberOfInstruments(); ++k)
    {
        novac::CString serial(s_setup.m_instrument[k].m_serial);

        Common::CopyFile(common.m_exePath + "configuration/" + serial + ".exml", confCopyDir + serial + ".exml");
    }
}

// This is the starting point for all the processing modes.
void CalculateAllFluxes(CContinuationOfProcessing continuation)
{
    try
    {
        Common common;
        s_setup.m_executableDirectory = common.m_exePath.std_str();

        CPostProcessing post{ g_logger, s_setup, g_userSettings, continuation };

        // Copy the settings that we have read in from the 'configuration' directory
        //  to the output directory to make it easier for the user to remember 
        //  what has been done...
        ArchiveSettingsFiles(g_userSettings);

        // Do the post-processing
        if (g_userSettings.m_processingMode == PROCESSING_MODE::PROCESSING_MODE_COMPOSITION)
        {
            ShowMessage("Warning: Post processing of composition measurements is not yet fully implemented");
            post.DoPostProcessing_Flux(); // this uses the same code as the flux processing
        }
        else if (g_userSettings.m_processingMode == PROCESSING_MODE::PROCESSING_MODE_INSTRUMENT_CALIBRATION)
        {
            post.DoPostProcessing_InstrumentCalibration();
        }
        else if (g_userSettings.m_processingMode == PROCESSING_MODE::PROCESSING_MODE_STRATOSPHERE)
        {
            ShowMessage("Warning: Post processing of stratospheric measurements is not yet fully implemented");
            post.DoPostProcessing_Strat();
        }
        else
        {
            post.DoPostProcessing_Flux();
        }

        ShowMessage("-- Exit post processing --");
    }
    catch (Poco::FileNotFoundException& e)
    {
        ShowMessage("File not found exception: " + e.displayText());
        ShowMessage("-- Exit post processing --");
        return;
    }
    catch (std::invalid_argument& e)
    {
        std::stringstream msg;
        msg << "Invalid argument exception caught: " << e.what();
        ShowMessage(msg.str());
        ShowMessage("-- Exit post processing --");
        return;
    }
    catch (std::exception& e)
    {
        std::stringstream msg;
        msg << "General exception caught: " << e.what();
        ShowMessage(msg.str());
        ShowMessage("-- Exit post processing --");
        return;
    }
}
