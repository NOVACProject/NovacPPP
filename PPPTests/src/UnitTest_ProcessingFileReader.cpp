#include <PPPLib/File/ProcessingFileReader.h>
#include "catch.hpp"

namespace novac
{
static std::string GetTestDataDirectory()
{
#ifdef _MSC_VER
    return std::string("../testData/");
#else
    return std::string("testData/");
#endif // _MSC_VER 
}

static std::string GetProcessingConfigurationFile()
{
    return GetTestDataDirectory() + std::string("processing.xml");
}

TEST_CASE("ReadProcessingFile gives expected configuration", "[ProcessingFileReader][File]")
{
    novac::ConsoleLog logger;
    Configuration::CUserConfiguration resultingConfiguration;
    FileHandler::CProcessingFileReader sut{ logger };

    sut.ReadProcessingFile(
        GetProcessingConfigurationFile(),
        resultingConfiguration);

    // Expected settings
    {
        REQUIRE(4 == resultingConfiguration.m_maxThreadNum);

#ifdef _MSC_VER
        REQUIRE("~/Novac/Piton de la Fournaise/OutputFeb2017UTC/" == resultingConfiguration.m_outputDirectory.std_str());
        REQUIRE("~/Novac/Piton de la Fournaise/Temp/" == resultingConfiguration.m_tempDirectory.std_str());
#endif // _MSC_VER

        REQUIRE(ProcessingMode::Flux == resultingConfiguration.m_processingMode);
        REQUIRE(StandardMolecule::SO2 == resultingConfiguration.m_molecule);

        REQUIRE(novac::CDateTime(2017, 1, 29, 12, 50, 51) == resultingConfiguration.m_fromDate);
        REQUIRE(novac::CDateTime(2017, 3, 01, 23, 50, 51) == resultingConfiguration.m_toDate);

        REQUIRE("C:/Temp/" == resultingConfiguration.m_LocalDirectory);
        REQUIRE(resultingConfiguration.m_includeSubDirectories_Local);

        REQUIRE("ftp://129.16.35.206/piton_de_la_fournaise/" == resultingConfiguration.m_FTPDirectory);
        REQUIRE(resultingConfiguration.m_includeSubDirectories_FTP);

        REQUIRE(0.5 == resultingConfiguration.m_calcGeometry_CompletenessLimit);
        REQUIRE(3600 == resultingConfiguration.m_calcGeometryValidTime);
        REQUIRE(100 == resultingConfiguration.m_calcGeometry_MinDistance);

        REQUIRE(7200 == resultingConfiguration.m_dualBeam_ValidTime);
        REQUIRE(10.0 == resultingConfiguration.m_dualBeam_MaxWindSpeedError);

        REQUIRE(1 == resultingConfiguration.m_nFitWindowsToUse);
        REQUIRE(0 == resultingConfiguration.m_mainFitWindow);

        REQUIRE(false == resultingConfiguration.m_uploadResults);

        REQUIRE(Configuration::SKY_OPTION::USER_SUPPLIED == resultingConfiguration.sky.skyOption);
        REQUIRE("C:/Temp/Some_sky_spectrum.std" == resultingConfiguration.sky.skySpectrumFile);
    }
}

TEST_CASE("ReadProcessingFile ReadProcessingFile with invalid file path throws exception", "[ProcessingFileReader][File]")
{
    novac::ConsoleLog logger;
    Configuration::CUserConfiguration resultingConfiguration;
    FileHandler::CProcessingFileReader sut{ logger };

    REQUIRE_THROWS(sut.ReadProcessingFile("./some-not-existing-file.txtk", resultingConfiguration));
}
}
