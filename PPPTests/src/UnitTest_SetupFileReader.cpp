#include <PPPLib/File/SetupFileReader.h>
#include "catch.hpp"
#include "StdOutLogger.h"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4505)
#endif

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
    return GetTestDataDirectory() + std::string("setup.xml");
}

TEST_CASE("SetupFileReader ReadSetupFile gives expected configuration", "[SetupFileReader][File]")
{
    StdOutLogger logger;
    Configuration::CNovacPPPConfiguration resultingConfiguration;
    FileHandler::CSetupFileReader sut{ logger };

    sut.ReadSetupFile(
        GetProcessingConfigurationFile(),
        resultingConfiguration);

    REQUIRE(4 == resultingConfiguration.NumberOfInstruments());

    {
        const auto* instrumentConfiguration = resultingConfiguration.GetInstrument(std::string("D2J2134"));
        REQUIRE(nullptr != instrumentConfiguration);
        REQUIRE(1 == instrumentConfiguration->m_location.GetLocationNum());
        REQUIRE_NOTHROW(instrumentConfiguration->m_location.CheckSettings());
    }

    {
        const auto* instrumentConfiguration = resultingConfiguration.GetInstrument(std::string("I2J8550"));
        REQUIRE(nullptr != instrumentConfiguration);
        REQUIRE(1 == instrumentConfiguration->m_location.GetLocationNum());
        REQUIRE_NOTHROW(instrumentConfiguration->m_location.CheckSettings());
    }

    {
        const auto* instrumentConfiguration = resultingConfiguration.GetInstrument(std::string("I2J8549"));
        REQUIRE(nullptr != instrumentConfiguration);
        REQUIRE(1 == instrumentConfiguration->m_location.GetLocationNum());
        REQUIRE_NOTHROW(instrumentConfiguration->m_location.CheckSettings());
    }

    {
        const auto* instrumentConfiguration = resultingConfiguration.GetInstrument(std::string("I2J8552"));
        REQUIRE(nullptr != instrumentConfiguration);
        REQUIRE(1 == instrumentConfiguration->m_location.GetLocationNum());
        REQUIRE_NOTHROW(instrumentConfiguration->m_location.CheckSettings());
    }
}

TEST_CASE("SetupFileReader ReadSetupFile with invalid file path throws exception", "[SetupFileReader][File]")
{
    StdOutLogger logger;
    Configuration::CNovacPPPConfiguration resultingConfiguration;
    FileHandler::CSetupFileReader sut{ logger };

    REQUIRE_THROWS(sut.ReadSetupFile("./some-non-exisisting-file.txt", resultingConfiguration));
}

}