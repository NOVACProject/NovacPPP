#include <PPPLib/Flux/FluxCalculator.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Meteorology/WindDataBase.h>
#include <PPPLib/Geometry/PlumeHeight.h>
#include <SpectralEvaluation/File/File.h>
#include "catch.hpp"
#include <iostream>


static std::string GetTestDataDirectory()
{
#ifdef _MSC_VER
    return std::string("../testData/");
#else
    return std::string("testData/");
#endif // _MSC_VER 
}

// Region Helper methods


// Endregion Helper methods

TEST_CASE("CalculateFlux, valid scan calculates flux (Ruahepu, Avantes)", "[CFluxCalculator][CalculateFlux][IntegrationTest][Avantes]")
{
    // Arrange
    const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0.txt";
    novac::ConsoleLog logger;

    Configuration::CUserConfiguration userSettings;

    Configuration::CInstrumentLocation instrumentLocation;
    instrumentLocation.m_spectrometerModel = "AVASPEC";
    instrumentLocation.m_coneangle = 60.0;

    Configuration::CInstrumentConfiguration instrumentConfiguration;
    instrumentConfiguration.m_serial = "2002128M1";
    instrumentConfiguration.m_location.InsertLocation(instrumentLocation);

    Configuration::CNovacPPPConfiguration configuration;

    Meteorology::CWindDataBase windDataBase;
    Geometry::CPlumeHeight plumeAltitude;

    novac::LogContext context;

    Flux::CFluxCalculator sut(logger, configuration, userSettings);

    SECTION("Instrument does not have a configured location. Returns error")
    {
        Flux::CFluxResult fluxResult;

        // Act
        bool success = sut.CalculateFlux(context, filename, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(!success);
    }

    SECTION("Instrument does not have a configured wind field. Returns error")
    {
        Flux::CFluxResult fluxResult;
        configuration.m_instrument.push_back(instrumentConfiguration);

        // Act
        bool success = sut.CalculateFlux(context, filename, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(!success);
    }

    SECTION("Returns success")
    {
        Flux::CFluxResult fluxResult;
        configuration.m_instrument.push_back(instrumentConfiguration);

        double windDirection = 10.0;
        double windSpeed = 277.3;
        auto defaultSource = Meteorology::MET_SOURCE::MET_DEFAULT;
        novac::CDateTime validFrom(2020, 01, 01, 00, 00, 00);
        novac::CDateTime validTo(9999, 12, 31, 23, 59, 59);
        double instrumentLatitude = -39.281302;
        double instrumentLongitude = 175.564254;
        double altitude = 2700;
        Meteorology::CWindField windField(windSpeed, defaultSource, windDirection, defaultSource, validFrom, validTo, instrumentLatitude, instrumentLongitude, altitude);
        windDataBase.InsertWindField(windField);

        // Act
        bool success = sut.CalculateFlux(context, filename, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(success == true);
    }
}
