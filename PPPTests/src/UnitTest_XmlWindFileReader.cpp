#include <PPPLib/Meteorology/XMLWindFileReader.h>
#include <PPPLib/Configuration/UserConfiguration.h>
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

static std::string GetWindFieldFile()
{
    return GetTestDataDirectory() + std::string("wind_ruapehu.wxml");
}

TEST_CASE("ReadWindFile gives expected wind profile", "[XMLWindFileReader][File]")
{
    novac::ConsoleLog logger;
    novac::LogContext context;
    Configuration::CUserConfiguration userConfiguration;
    Meteorology::CWindDataBase resultingDatabase;
    FileHandler::CXMLWindFileReader sut{ logger, userConfiguration };

    sut.ReadWindFile(
    context,
        GetWindFieldFile(),
        resultingDatabase);

    // Expected properties of the database
    {
        REQUIRE(resultingDatabase.m_dataBaseName == "Ruapehu");

        REQUIRE(2 == resultingDatabase.GetDataBaseSize());
    }

    // Expected Wind fields
    {
        auto time = novac::CDateTime{ 2023, 1, 19, 22, 0, 0 };
        auto location = novac::CGPSData{ -39.281302 , 175.564254 , 2700.0 };
        auto windField = Meteorology::CWindField{};

        // Act
        bool exists = resultingDatabase.GetWindField(time, location, Meteorology::INTERP_EXACT, windField);

        // Assert
        REQUIRE(exists);
        REQUIRE(Approx(windField.GetWindSpeed()) == 5.13);
        REQUIRE(Approx(windField.GetWindDirection()) == 277.3);
    }

    {
        auto time = novac::CDateTime{ 2023, 1, 20, 1, 0, 0 };
        auto location = novac::CGPSData{ -39.281302 , 175.564254 , 2700.0 };
        auto windField = Meteorology::CWindField{};

        // Act
        bool exists = resultingDatabase.GetWindField(time, location, Meteorology::INTERP_EXACT, windField);

        // Assert
        REQUIRE(exists);
        REQUIRE(Approx(windField.GetWindSpeed()) == 7.56);
        REQUIRE(Approx(windField.GetWindDirection()) == 278.4);
    }

    // Requests of data outside of the time range retrieves 'not found'
    {
        auto time = novac::CDateTime{ 2023, 1, 21, 1, 0, 0 }; // Time out of range
        auto location = novac::CGPSData{ -39.281302 , 175.564254 , 2700.0 };
        auto windField = Meteorology::CWindField{};

        // Act
        bool exists = resultingDatabase.GetWindField(time, location, Meteorology::INTERP_EXACT, windField);

        // Assert
        REQUIRE(exists == false);
        REQUIRE(Approx(windField.GetWindSpeed()) == 10.0); // Default wind speed
        REQUIRE(Approx(windField.GetWindDirection()) == 0.0); // Default wind direction
    }

    {
        auto time = novac::CDateTime{ 2023, 1, 18, 1, 0, 0 }; // Time out of range
        auto location = novac::CGPSData{ -39.281302 , 175.564254 , 2700.0 };
        auto windField = Meteorology::CWindField{};

        // Act
        bool exists = resultingDatabase.GetWindField(time, location, Meteorology::INTERP_EXACT, windField);

        // Assert
        REQUIRE(exists == false);
        REQUIRE(Approx(windField.GetWindSpeed()) == 10.0); // Default wind speed
        REQUIRE(Approx(windField.GetWindDirection()) == 0.0); // Default wind direction
    }
}
}
