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

TEST_CASE("CalculateFlux, very low intensity in evaluation log, cannot calculate flux", "[CFluxCalculator][CalculateFlux][IntegrationTest][Avantes]")
{
    // Arrange
    Evaluation::CExtendedScanResult evaluationResult;
    evaluationResult.m_evalLogFile.push_back(GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0.txt");
    evaluationResult.m_startTime = novac::CDateTime(2023, 1, 20, 19, 07, 00);
    novac::ConsoleLog logger;
    novac::LogContext context;

    Configuration::CUserConfiguration userSettings;

    Configuration::CInstrumentLocation instrumentLocation;
    instrumentLocation.m_spectrometerModel = "AVASPEC";
    instrumentLocation.m_coneangle = 60.0;

    Configuration::CInstrumentConfiguration instrumentConfiguration;
    instrumentConfiguration.m_serial = "2002128M1";
    instrumentConfiguration.m_location.InsertLocation(instrumentLocation);

    Configuration::CNovacPPPConfiguration configuration;

    Meteorology::CWindDataBase windDataBase;
    Geometry::PlumeHeight plumeAltitude;

    Flux::CFluxCalculator sut(logger, configuration, userSettings);

    SECTION("Instrument does not have a configured location. Returns error")
    {
        Flux::FluxResult fluxResult;

        // Act
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(!success);
    }

    SECTION("Instrument does not have a configured wind field. Returns error")
    {
        Flux::FluxResult fluxResult;
        configuration.m_instrument.push_back(instrumentConfiguration);

        // Act
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(!success);
    }

    SECTION("All data required, returns failure")
    {
        Flux::FluxResult fluxResult;
        configuration.m_instrument.push_back(instrumentConfiguration);

        double windDirection = 10.0;
        double windSpeed = 277.3;
        auto defaultSource = Meteorology::MeteorologySource::Default;
        novac::CDateTime validFrom(2020, 01, 01, 00, 00, 00);
        novac::CDateTime validTo(9999, 12, 31, 23, 59, 59);
        double instrumentLatitude = -39.281302;
        double instrumentLongitude = 175.564254;
        double altitude = 2700;
        Meteorology::WindField windField(windSpeed, defaultSource, windDirection, defaultSource, validFrom, validTo, instrumentLatitude, instrumentLongitude, altitude);
        windDataBase.InsertWindField(windField);

        // Act
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(!success);
    }
}

TEST_CASE("CalculateFlux, reads values from configuration and verifies input", "[CFluxCalculator][CalculateFlux][IntegrationTest][Avantes][2002128M1_230120_1907_0_ReEvaluation]")
{
    // Arrange
    Evaluation::CExtendedScanResult evaluationResult;
    evaluationResult.m_evalLogFile.push_back(GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0_ReEvaluation.txt");
    evaluationResult.m_instrumentSerial = "2002128M1";
    evaluationResult.m_startTime = novac::CDateTime(2023, 1, 20, 19, 07, 00);
    novac::ConsoleLog logger;
    novac::LogContext context;

    Configuration::CUserConfiguration userSettings;

    Configuration::CInstrumentLocation instrumentLocation;
    instrumentLocation.m_spectrometerModel = "AVASPEC";
    instrumentLocation.m_coneangle = 60.0;
    instrumentLocation.m_compass = 123.456;

    Configuration::CInstrumentConfiguration instrumentConfiguration;
    instrumentConfiguration.m_serial = "2002128M1";

    Configuration::CNovacPPPConfiguration configuration;

    Meteorology::CWindDataBase windDataBase;
    Geometry::PlumeHeight plumeAltitude;

    SECTION("Instrument does not have a configured location. Returns error")
    {
        Flux::FluxResult fluxResult;
        Flux::CFluxCalculator sut(logger, configuration, userSettings);

        // Act
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(!success);
    }

    SECTION("Instrument does not have a configured wind field. Returns error")
    {
        Flux::FluxResult fluxResult;
        instrumentConfiguration.m_location.InsertLocation(instrumentLocation);
        configuration.m_instrument.push_back(instrumentConfiguration);
        Flux::CFluxCalculator sut(logger, configuration, userSettings);

        // Act
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(!success);
    }

    SECTION("All data required, completeness below limit. Returns failure")
    {
        Flux::FluxResult fluxResult;

        instrumentConfiguration.m_location.InsertLocation(instrumentLocation);
        configuration.m_instrument.push_back(instrumentConfiguration);

        double windDirection = 10.0;
        double windSpeed = 277.3;
        auto defaultSource = Meteorology::MeteorologySource::Default;
        novac::CDateTime validFrom(2020, 01, 01, 00, 00, 00);
        novac::CDateTime validTo(9999, 12, 31, 23, 59, 59);
        double instrumentLatitude = -39.281302;
        double instrumentLongitude = 175.564254;
        double altitude = 2700;
        Meteorology::WindField windField(windSpeed, defaultSource, windDirection, defaultSource, validFrom, validTo, instrumentLatitude, instrumentLongitude, altitude);
        windDataBase.InsertWindField(windField);

        userSettings.m_completenessLimitFlux = 0.90; // this is higher than the completess of the scan
        Flux::CFluxCalculator sut(logger, configuration, userSettings);

        // Act
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(!success);
    }

    SECTION("Takes instrument location from configuration, not from evaluation log file")
    {
        instrumentConfiguration.m_location.InsertLocation(instrumentLocation);
        configuration.m_instrument.push_back(instrumentConfiguration);

        double windDirection = 262.3;
        double windSpeed = 10.54;
        auto defaultSource = Meteorology::MeteorologySource::Default;
        novac::CDateTime validFrom(2020, 01, 01, 00, 00, 00);
        novac::CDateTime validTo(9999, 12, 31, 23, 59, 59);
        double instrumentLatitude = -39.281302;
        double instrumentLongitude = 175.564254;
        double altitude = 2700;
        Meteorology::WindField windField(windSpeed, defaultSource, windDirection, defaultSource, validFrom, validTo, instrumentLatitude, instrumentLongitude, altitude);
        windDataBase.InsertWindField(windField);

        plumeAltitude.m_plumeAltitude = 3500;

        userSettings.m_completenessLimitFlux = 0.80;
        Flux::CFluxCalculator sut(logger, configuration, userSettings);

        // Act
        Flux::FluxResult fluxResult;
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(success);

        REQUIRE(fluxResult.m_compass == 123.456); // this value should be taken frm the configuration, not from the contents of the scan.
        REQUIRE(fluxResult.m_coneAngle == 60.0); // same value as in the evaluation log file
        REQUIRE(fluxResult.m_instrument == "2002128M1"); // same value as in the evaluation log file
        REQUIRE(fluxResult.m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG); // same value as in the evaluation log file
    }
}

TEST_CASE("CalculateFlux, valid scan with column values in molec/cm2 calculates flux", "[CFluxCalculator][CalculateFlux][IntegrationTest][Avantes][2002128M1_230120_1907_0_ReEvaluation]")
{
    // Arrange
    Evaluation::CExtendedScanResult evaluationResult;
    evaluationResult.m_evalLogFile.push_back(GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0_ReEvaluation.txt");
    evaluationResult.m_instrumentSerial = "2002128M1";
    evaluationResult.m_startTime = novac::CDateTime(2023, 1, 20, 19, 07, 00);
    novac::ConsoleLog logger;
    novac::LogContext context;

    Configuration::CUserConfiguration userSettings;
    userSettings.m_completenessLimitFlux = 0.80;

    Configuration::CInstrumentLocation instrumentLocation;
    instrumentLocation.m_spectrometerModel = "AVASPEC";
    instrumentLocation.m_coneangle = 60.0;
    instrumentLocation.m_compass = 266.0;
    instrumentLocation.m_coneangle = 60.0;
    instrumentLocation.m_altitude = 2700;

    Configuration::CInstrumentConfiguration instrumentConfiguration;
    instrumentConfiguration.m_serial = "2002128M1";

    Configuration::CNovacPPPConfiguration configuration;

    double windDirection = 262.3;
    double windSpeed = 10.54;
    auto defaultSource = Meteorology::MeteorologySource::EcmwfForecast;
    novac::CDateTime validFrom(2020, 01, 01, 00, 00, 00);
    novac::CDateTime validTo(9999, 12, 31, 23, 59, 59);
    double instrumentLatitude = -39.281302;
    double instrumentLongitude = 175.564254;
    double altitude = 2700;
    Meteorology::WindField windField(windSpeed, defaultSource, windDirection, defaultSource, validFrom, validTo, instrumentLatitude, instrumentLongitude, altitude);
    Meteorology::CWindDataBase windDataBase;

    Geometry::PlumeHeight plumeAltitude;
    plumeAltitude.m_plumeAltitudeSource = Meteorology::MeteorologySource::GeometryCalculationTwoInstruments;
    plumeAltitude.m_plumeAltitude = 3500;

    SECTION("Returns calculated flux")
    {
        plumeAltitude.m_plumeAltitudeError = 0.0;

        instrumentConfiguration.m_location.InsertLocation(instrumentLocation);
        configuration.m_instrument.push_back(instrumentConfiguration);

        windDataBase.InsertWindField(windField);

        Flux::CFluxCalculator sut(logger, configuration, userSettings);

        // Act
        Flux::FluxResult fluxResult;
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(success);

        REQUIRE(fluxResult.m_compass == 266.0);
        REQUIRE(fluxResult.m_coneAngle == 60.0);
        REQUIRE(fluxResult.m_instrument == "2002128M1");
        REQUIRE(fluxResult.m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG);

        REQUIRE(fluxResult.m_completeness == Approx(0.852).margin(0.001));
        REQUIRE(fluxResult.m_plumeCentre[0] == Approx(50.0).margin(0.1));
        REQUIRE(fluxResult.m_plumeCentre[1] == 0.0);
        REQUIRE(fluxResult.m_plumeHeight.m_plumeAltitude == 800); // altitude of the plume - altitude of the instrument
        REQUIRE(fluxResult.m_scanOffset == Approx(-1.887e17).margin(1e15));
        REQUIRE(fluxResult.m_numGoodSpectra == 51);
        REQUIRE(fluxResult.m_startTime == novac::CDateTime(2023, 01, 20, 19, 07, 48));
        REQUIRE(fluxResult.m_stopTime == novac::CDateTime(2023, 01, 20, 19, 15, 41));
        REQUIRE(fluxResult.m_tilt == 0.0);

        REQUIRE(Approx(windSpeed) == fluxResult.m_windField.GetWindSpeed());
        REQUIRE(Approx(windDirection) == fluxResult.m_windField.GetWindDirection());

        // The flux here has been verified to match with the result from the NovacProgram
        REQUIRE(fluxResult.m_flux == Approx(5.068).margin(0.001)); // kg/s
        REQUIRE(fluxResult.m_fluxError_Wind == Approx(0.0)); // no given error in the wind speed => no given error in the flux due to wind error.
        REQUIRE(fluxResult.m_fluxError_PlumeHeight == Approx(0.0)); // zero, since we set the error in plume height to zero above
        REQUIRE(fluxResult.m_fluxQualityFlag == FluxQuality::Yellow); // result is ok
    }

    SECTION("Calculates error in flux from error in wind speed")
    {
        plumeAltitude.m_plumeAltitudeError = 100.0;
        instrumentConfiguration.m_location.InsertLocation(instrumentLocation);
        configuration.m_instrument.push_back(instrumentConfiguration);

        windField.SetWindSpeedError(5.0);
        windDataBase.InsertWindField(windField);

        Flux::CFluxCalculator sut(logger, configuration, userSettings);

        // Act
        Flux::FluxResult fluxResult;
        bool success = sut.CalculateFlux(context, evaluationResult, windDataBase, plumeAltitude, fluxResult);

        // Assert
        REQUIRE(success);

        REQUIRE(Approx(windSpeed) == fluxResult.m_windField.GetWindSpeed());
        REQUIRE(Approx(windDirection) == fluxResult.m_windField.GetWindDirection());

        // The flux here has been verified to match with the result from the NovacProgram
        REQUIRE(fluxResult.m_flux == Approx(5.068).margin(0.001)); // kg/s
        REQUIRE(fluxResult.m_fluxError_Wind == Approx(2.40).margin(0.01));
        REQUIRE(fluxResult.m_fluxError_PlumeHeight == Approx(0.6).margin(0.05)); // as calculated from the error in plume altitude
        REQUIRE(fluxResult.m_fluxQualityFlag == FluxQuality::Yellow); // result is ok
    }
}
