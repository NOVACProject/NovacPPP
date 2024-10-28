#include <PPPLib/Geometry/GeometryCalculator.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include "catch.hpp"

namespace novac
{

TEST_CASE("CalculateGeometry with PlumeInScanProperty gives expected wind direction and plume height", "[GeometryCalculator][Geometry]")
{
    novac::ConsoleLog log;
    Configuration::CUserConfiguration userSettings;
    Geometry::CGeometryCalculator sut(log, userSettings);

    novac::CPlumeInScanProperty plume1;
    novac::CPlumeInScanProperty plume2;
    novac::CDateTime startTime1(2023, 01, 20, 15, 16, 30);
    novac::CDateTime startTime2(2023, 01, 20, 15, 17, 30);
    Configuration::CInstrumentLocation locations[2];
    // location0 is almost directly east of the volcano
    locations[0].m_latitude = -39.277528;
    locations[0].m_longitude = 175.608731;
    locations[0].m_altitude = 1756;
    locations[0].m_compass = 266.0;
    locations[0].m_coneangle = 60.0;
    locations[0].m_volcano = "ruapehu";
    // location1 is almost directly north of the volcano
    locations[1].m_latitude = -39.237137;
    locations[1].m_longitude = 175.556395;
    locations[1].m_altitude = 1633;
    locations[1].m_compass = 172.0;
    locations[1].m_coneangle = 60.0;
    locations[1].m_volcano = "ruapehu";

    SECTION("Plume directly above instrument one gives wind direction from volcano to instrument one")
    {
        plume1.plumeCenter = 0.0; // directly above this instrument, meaning the plume goes from the volcano to the location of this instrument
        plume1.plumeCenterError = 2.0;
        plume2.plumeCenter = -75.0; // the calculated plume height and wind direction is also determined from this angle, so this needs to be reasonable...
        plume2.plumeCenterError = 2.0;

        // Act
        Geometry::CGeometryResult result;
        const bool returnValue = sut.CalculateGeometry(plume1, startTime1, plume2, startTime2, locations, result);

        // Assert
        REQUIRE(returnValue == true);
        REQUIRE(261.6 == Approx(result.m_windDirection).margin(1.0));
        REQUIRE(2 == Approx(result.m_windDirectionError).margin(1.0));
        REQUIRE(3826 == Approx(result.m_plumeAltitude).margin(10.0));
        REQUIRE(280 == Approx(result.m_plumeAltitudeError).margin(5.0));

        // The other properties should be filled in as well
        REQUIRE(0.0 == Approx(result.m_plumeCentre1.Value()));
        REQUIRE(-75.0 == Approx(result.m_plumeCentre2.Value()));
        REQUIRE(CDateTime(2023, 01, 20, 15, 17, 0) == result.m_averageStartTime);
        REQUIRE(60 == result.m_startTimeDifference);
    }

    SECTION("Plume directly above instrument two gives wind direction from volcano to instrument one")
    {
        plume1.plumeCenter = 80.0; // the calculated plume height and wind direction is also determined from this angle, so this needs to be reasonable...
        plume1.plumeCenterError = 2.0;
        plume2.plumeCenter = 0.0; // directly above this instrument, meaning the plume goes from the volcano to the location of this instrument
        plume2.plumeCenterError = 2.0;

        // Act
        Geometry::CGeometryResult result;
        const bool returnValue = sut.CalculateGeometry(plume1, startTime1, plume2, startTime2, locations, result);

        // Assert
        REQUIRE(returnValue == true);
        REQUIRE(162.1 == Approx(result.m_windDirection).margin(1.0));
        REQUIRE(2 == Approx(result.m_windDirectionError).margin(1.0));
        REQUIRE(3349 == Approx(result.m_plumeAltitude).margin(10.0));
        REQUIRE(365 == Approx(result.m_plumeAltitudeError).margin(5.0));

        // The other properties should be filled in as well
        REQUIRE(80.0 == Approx(result.m_plumeCentre1.Value()));
        REQUIRE(0.0 == Approx(result.m_plumeCentre2.Value()));
        REQUIRE(CDateTime(2023, 01, 20, 15, 17, 0) == result.m_averageStartTime);
        REQUIRE(60 == result.m_startTimeDifference);
    }
}

TEST_CASE("CalculateWindDirection with PlumeInScanProperty gives expected wind direction and plume height", "[GeometryCalculator][Geometry]")
{
    novac::ConsoleLog log;
    Configuration::CUserConfiguration userSettings;
    Geometry::CGeometryCalculator sut(log, userSettings);

    novac::CPlumeInScanProperty plume;
    novac::CDateTime startTime(2023, 01, 20, 15, 16, 30);
    Configuration::CInstrumentLocation location;
    // location0 is almost directly east of the volcano
    location.m_latitude = -39.277528;
    location.m_longitude = 175.608731;
    location.m_altitude = 1756;
    location.m_compass = 266.0;
    location.m_coneangle = 60.0;
    location.m_volcano = "ruapehu";
    Geometry::CPlumeHeight plumeHeight;
    plumeHeight.m_plumeAltitude = 3600;

    SECTION("Plume directly above instrument one gives wind direction from volcano to instrument one")
    {
        plume.plumeCenter = 0.0; // directly above this instrument, meaning the plume goes from the volcano to the location of this instrument
        plume.plumeCenterError = 2.0;

        // Act
        Geometry::CGeometryResult result;
        const bool returnValue = sut.CalculateWindDirection(plume, startTime, plumeHeight, location, result);

        // Assert
        REQUIRE(returnValue == true);
        REQUIRE(261.6 == Approx(result.m_windDirection).margin(1.0));
        REQUIRE(4 == Approx(result.m_windDirectionError).margin(1.0));
        REQUIRE(NOT_A_NUMBER == result.m_plumeAltitude);

        // The other properties should be filled in as well
        REQUIRE(0.0 == Approx(result.m_plumeCentre1.Value()));
        REQUIRE(false == result.m_plumeCentre2.HasValue());
        REQUIRE(CDateTime(2023, 01, 20, 15, 16, 30) == result.m_averageStartTime);
        REQUIRE(0 == result.m_startTimeDifference);
    }



}
}