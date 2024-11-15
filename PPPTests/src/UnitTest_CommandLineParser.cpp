#include <PPPLib/Configuration/CommandLineParser.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/VolcanoInfo.h>
#include "catch.hpp"

namespace Configuration
{
TEST_CASE("FromDate overrides set from date", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--FromDate=2024.05.31" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_fromDate == novac::CDateTime(2024, 5, 31, 0, 0, 0));
}

TEST_CASE("ToDate overrides set to date", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--ToDate=2024.05.31" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_fromDate == novac::CDateTime(2005, 10, 1, 0, 0, 0)); // Default value
    REQUIRE(userSettings.m_toDate == novac::CDateTime(2024, 5, 31, 0, 0, 0));
}

TEST_CASE("Volcano overrides set volcano - common name", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--Volcano=masaya" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_volcano == 8); // The index for Masaya

    novac::CString volcanoName;
    volcanoes.GetVolcanoName(userSettings.m_volcano, volcanoName);
    REQUIRE(volcanoName == "Masaya");
}

TEST_CASE("Volcano overrides set volcano - simple name", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--Volcano=nevado_del_ruiz" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_volcano == 10); // The index for Nevado del Ruiz

    novac::CString volcanoName;
    volcanoes.GetVolcanoName(userSettings.m_volcano, volcanoName);
    REQUIRE(volcanoName == "Nevado del Ruiz");
}

TEST_CASE("Volcano overrides set volcano - volcano number", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--Volcano=0101-04=" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_volcano == 23); // The index for Stromboli

    novac::CString volcanoName;
    volcanoes.GetVolcanoName(userSettings.m_volcano, volcanoName);
    REQUIRE(volcanoName == "Stromboli");
}

TEST_CASE("Workdir overrides set workdir", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--Workdir=/home/novacuser/novacppp" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(setExePath == "/home/novacuser/novacppp/");
}

TEST_CASE("outputdirectory overrides set default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--outputdirectory=/home/novacuser/novacppp/output/" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_outputDirectory == "/home/novacuser/novacppp/output/");
}

TEST_CASE("tempdirectory overrides set default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--tempdirectory=/home/novacuser/novacppp/temp/" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_tempDirectory == "/home/novacuser/novacppp/temp/");
}

TEST_CASE("WindFieldFile overrides set default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--WindFieldFile=/home/novacuser/novacppp/wind.wxml" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_windFieldFile == ""); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_windFieldFile == "/home/novacuser/novacppp/wind.wxml");
}

TEST_CASE("MaxThreadNum overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--MaxThreadNum=73" };
    Configuration::CUserConfiguration userSettings;
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_maxThreadNum == 73);
}

TEST_CASE("IncludeSubDirs_Local overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--IncludeSubDirs_Local=0" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_includeSubDirectories_Local); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_includeSubDirectories_Local == false);
}

TEST_CASE("IncludeSubDirs_FTP overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--IncludeSubDirs_FTP=0" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_includeSubDirectories_FTP); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_includeSubDirectories_FTP == false);
}

TEST_CASE("FTPDirectory overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--FTPDirectory=/data/novacppp" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_FTPDirectory == ""); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_FTPDirectory == "/data/novacppp");
}

TEST_CASE("FTPUsername overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--FTPUsername=the-novac-user" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_FTPUsername == ""); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_FTPUsername == "the-novac-user");
}

TEST_CASE("FTPPassword overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--FTPPassword=the-novac-user-password" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_FTPPassword == ""); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_FTPPassword == "the-novac-user-password");
}

TEST_CASE("UploadResults overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--UploadResults=1" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_uploadResults == false); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_uploadResults == true);
}

TEST_CASE("mode overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--mode=2" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_processingMode == ProcessingMode::Flux); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_processingMode == ProcessingMode::Stratosphere);
}

TEST_CASE("molecule O3 overrides default", "[CommandLineParser][Configuration]")
{
    // Arrange
    std::string setExePath;
    std::vector<std::string>arguments = { "--molecule=O3" };
    Configuration::CUserConfiguration userSettings;
    REQUIRE(userSettings.m_molecule == novac::StandardMolecule::SO2); // check assumption here
    novac::CVolcanoInfo volcanoes;
    novac::ConsoleLog logger;

    // Act
    CommandLineParser::ParseCommandLineOptions(arguments, userSettings, volcanoes, setExePath, logger);

    // Assert
    REQUIRE(userSettings.m_molecule == novac::StandardMolecule::O3);
}
}