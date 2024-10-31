#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <cstring>
#include "catch.hpp"

namespace novac
{

TEST_CASE("CNovacPPPConfiguration GetInstrumentLocation returns expected value", "[CNovacPPPConfiguration][Configuration]")
{
    const std::string instrumentSerial = "I2J5678";

    Configuration::CInstrumentLocation configuredLocation;
    configuredLocation.m_altitude = 1633;
    configuredLocation.m_compass = 172.0;
    configuredLocation.m_coneangle = 60.0;
    configuredLocation.m_instrumentType = NovacInstrumentType::Gothenburg;
    configuredLocation.m_latitude = -39.237137;
    configuredLocation.m_longitude = 175.556395;
    configuredLocation.m_locationName = "RUD02_Iwikau_Village";
    configuredLocation.m_spectrometerModel = "Avaspec";
    configuredLocation.m_tilt = 0.0;
    configuredLocation.m_volcano = "Ruapehu";
    configuredLocation.m_validFrom = CDateTime(2022, 05, 06, 0, 0, 0);
    configuredLocation.m_validTo = CDateTime(9999, 01, 01, 0, 0, 0);

    SECTION("No instruments configured, throws NotFoundException with expected message")
    {
        // Arrange
        const CDateTime searchTime{ 2024, 03, 14, 15, 16, 17 };
        Configuration::CNovacPPPConfiguration sut;

        // Act & Assert
        try
        {
            sut.GetInstrumentLocation(instrumentSerial, searchTime);
            REQUIRE(false); // failure
        }
        catch (novac::NotFoundException& ex)
        {
            REQUIRE(strstr(ex.message.c_str(), "Cannot find configuration for instrument with serial number ") != nullptr);
        }
    }

    SECTION("One instrument configured - Query done within instrument valid time - Returns instrument location")
    {
        Configuration::CLocationConfiguration configuredInstrumentLocation;
        configuredInstrumentLocation.InsertLocation(configuredLocation);
        Configuration::CInstrumentConfiguration configuredInstrument;
        configuredInstrument.m_serial = instrumentSerial;
        configuredInstrument.m_location = configuredInstrumentLocation;
        const CDateTime searchTime{ 2022, 05, 06, 15, 16, 17 };

        Configuration::CNovacPPPConfiguration sut;
        sut.m_instrument.push_back(configuredInstrument);

        // Act
        Configuration::CInstrumentLocation result = sut.GetInstrumentLocation(instrumentSerial, searchTime);

        // Assert
        REQUIRE(std::abs(result.m_latitude + 39.237137) < 1e-9);
        REQUIRE(std::abs(result.m_longitude - 175.556395) < 1e-9);
        REQUIRE(result.m_altitude == 1633);
    }

    SECTION("One instrument configured - Query done outside of instrument valid time - Throws NotFoundException")
    {
        Configuration::CLocationConfiguration configuredInstrumentLocation;
        configuredInstrumentLocation.InsertLocation(configuredLocation);
        Configuration::CInstrumentConfiguration configuredInstrument;
        configuredInstrument.m_serial = instrumentSerial;
        configuredInstrument.m_location = configuredInstrumentLocation;
        const CDateTime searchTime{ 2022, 05, 05, 15, 16, 17 }; // the day before the instrument was installed

        Configuration::CNovacPPPConfiguration sut;
        sut.m_instrument.push_back(configuredInstrument);

        // Act & Assert
        try
        {
            sut.GetInstrumentLocation(instrumentSerial, searchTime);
            REQUIRE(false); // failure
        }
        catch (novac::NotFoundException& ex)
        {
            REQUIRE(strstr(ex.message.c_str(), "does not have a configured location on 2022.05.05") != nullptr);
        }
    }
}

TEST_CASE("CNovacPPPConfiguration GetFitWindow returns expected value", "[CNovacPPPConfiguration][Configuration]")
{
    const std::string instrumentSerial = "I2J5678";

    CDateTime fitWindowValidFrom(2022, 05, 05, 0, 0, 0);
    CDateTime fitWindowValidTo(9999, 01, 01, 0, 0, 0);
    novac::CFitWindow configuredFitWindow;
    configuredFitWindow.channel = 0;
    configuredFitWindow.fitLow = 464;
    configuredFitWindow.fitHigh = 630;

    SECTION("No instruments configured, throws NotFoundException with expected message")
    {
        // Arrange
        const CDateTime searchTime{ 2024, 03, 14, 15, 16, 17 };
        Configuration::CNovacPPPConfiguration sut;

        // Act & Assert
        try
        {
            sut.GetFitWindow(instrumentSerial, 0, searchTime);
            REQUIRE(false); // failure
        }
        catch (novac::NotFoundException& ex)
        {
            REQUIRE(strstr(ex.message.c_str(), "Cannot find configuration for instrument with serial number ") != nullptr);
        }
    }

    SECTION("One instrument configured - Query done within instrument valid time - Returns instrument location")
    {
        Configuration::CEvaluationConfiguration configuredInstrumentEvaluation;
        configuredInstrumentEvaluation.m_serial = instrumentSerial;
        configuredInstrumentEvaluation.InsertFitWindow(configuredFitWindow, fitWindowValidFrom, fitWindowValidTo);

        Configuration::CInstrumentConfiguration configuredInstrument;
        configuredInstrument.m_serial = instrumentSerial;
        configuredInstrument.m_eval = configuredInstrumentEvaluation;

        const CDateTime searchTime{ 2022, 05, 06, 15, 16, 17 };

        Configuration::CNovacPPPConfiguration sut;
        sut.m_instrument.push_back(configuredInstrument);

        // Act
        novac::CFitWindow result = sut.GetFitWindow(instrumentSerial, 0, searchTime);

        // Assert
        REQUIRE(result.fitLow == 464);
        REQUIRE(result.fitHigh == 630);
    }

    SECTION("One instrument configured - Query done outside of instrument valid time - Throws NotFoundException")
    {
        Configuration::CInstrumentConfiguration configuredInstrument;
        configuredInstrument.m_serial = instrumentSerial;
        const CDateTime searchTime{ 2022, 05, 05, 15, 16, 17 }; // the day before the instrument was installed

        Configuration::CNovacPPPConfiguration sut;
        sut.m_instrument.push_back(configuredInstrument);

        // Act & Assert
        try
        {
            sut.GetFitWindow(instrumentSerial, 0, searchTime);
            REQUIRE(false); // failure
        }
        catch (novac::NotFoundException& ex)
        {
            REQUIRE(strstr(ex.message.c_str(), "does not have a configured fit-window on 2022.05.05") != nullptr);
        }
    }

    SECTION("One instrument configured - Query done on wrong channel of instrument - Throws NotFoundException")
    {
        Configuration::CInstrumentConfiguration configuredInstrument;
        configuredInstrument.m_serial = instrumentSerial;
        const CDateTime searchTime{ 2022, 05, 07, 15, 16, 17 };

        Configuration::CNovacPPPConfiguration sut;
        sut.m_instrument.push_back(configuredInstrument);

        // Act & Assert
        try
        {
            sut.GetFitWindow(instrumentSerial, 1, searchTime);
            REQUIRE(false); // failure
        }
        catch (novac::NotFoundException& ex)
        {
            REQUIRE(strstr(ex.message.c_str(), "does not have a configured fit-window on 2022.05.07") != nullptr);
        }
    }
}

TEST_CASE("CNovacPPPConfiguration GetDarkCorrection returns expected value", "[CNovacPPPConfiguration][Configuration]")
{
    const std::string instrumentSerial = "I2J5678";

    Configuration::CInstrumentLocation configuredLocation;
    configuredLocation.m_altitude = 1633;
    configuredLocation.m_compass = 172.0;
    configuredLocation.m_coneangle = 60.0;
    configuredLocation.m_instrumentType = NovacInstrumentType::Gothenburg;
    configuredLocation.m_latitude = -39.237137;
    configuredLocation.m_longitude = 175.556395;
    configuredLocation.m_locationName = "RUD02_Iwikau_Village";
    configuredLocation.m_spectrometerModel = "Avaspec";
    configuredLocation.m_tilt = 0.0;
    configuredLocation.m_volcano = "Ruapehu";
    configuredLocation.m_validFrom = CDateTime(2022, 05, 06, 0, 0, 0);
    configuredLocation.m_validTo = CDateTime(9999, 01, 01, 0, 0, 0);

    SECTION("No instruments configured, throws NotFoundException with expected message")
    {
        // Arrange
        const CDateTime searchTime{ 2024, 03, 14, 15, 16, 17 };
        Configuration::CNovacPPPConfiguration sut;

        // Act & Assert
        try
        {
            sut.GetDarkCorrection(instrumentSerial, searchTime);
            REQUIRE(false); // failure
        }
        catch (novac::NotFoundException& ex)
        {
            REQUIRE(strstr(ex.message.c_str(), "Cannot find configuration for instrument with serial number ") != nullptr);
        }
    }

    SECTION("One instrument configured - No Dark settings explitly set - Query done within instrument valid time - Returns default instrument location")
    {
        Configuration::CLocationConfiguration configuredInstrumentLocation;
        configuredInstrumentLocation.InsertLocation(configuredLocation);
        Configuration::CInstrumentConfiguration configuredInstrument;
        configuredInstrument.m_serial = instrumentSerial;
        configuredInstrument.m_location = configuredInstrumentLocation;
        const CDateTime searchTime{ 2022, 05, 06, 15, 16, 17 };

        Configuration::CNovacPPPConfiguration sut;
        sut.m_instrument.push_back(configuredInstrument);

        // Act
        Configuration::CDarkSettings result = sut.GetDarkCorrection(instrumentSerial, searchTime);

        // Assert
        REQUIRE(result.m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MEASURED_IN_SCAN);
    }
}
}