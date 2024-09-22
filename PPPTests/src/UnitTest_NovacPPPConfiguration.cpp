#include <PPPLib/Configuration/NovacPPPConfiguration.h>
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
    configuredLocation.m_instrumentType = INSTRUMENT_TYPE::INSTR_GOTHENBURG;
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
        catch (PPPLib::NotFoundException& ex)
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

    SECTION("One instrument configured - Query done outside of instrument valid time - Returns non zero result")
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
        catch (PPPLib::NotFoundException& ex)
        {
            REQUIRE(strstr(ex.message.c_str(), "does not have a configured location on 2022.05.05") != nullptr);
        }
    }
}

}