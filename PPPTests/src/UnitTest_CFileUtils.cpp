#include "catch.hpp"
#include <PPPLib/MFC/CFileUtils.h>

namespace novac
{
TEST_CASE("GetInfoFromFileName behaves as expected", "[CFileUtils]")
{
    CDateTime start;
    CString serial;
    int channel;
    MEASUREMENT_MODE mode;

    SECTION("Finds correct serial")
    {
        CFileUtils::GetInfoFromFileName("D2J2134_170129_0317_1.pak", start, serial, channel, mode);
        REQUIRE(serial.std_str() == "D2J2134");

        CFileUtils::GetInfoFromFileName("I2J98765_170129_0317_1.pak", start, serial, channel, mode);
        REQUIRE(serial.std_str() == "I2J98765");
    }

    SECTION("Finds correct date")
    {
        CFileUtils::GetInfoFromFileName("D2J2134_170129_0317_1.pak", start, serial, channel, mode);

        REQUIRE(start.year == 2017);
        REQUIRE(start.month == 1);
        REQUIRE(start.day == 29);
        REQUIRE(start.hour == 3);
        REQUIRE(start.minute == 17);
        REQUIRE(start.second == 0);
    }

    SECTION("Finds correct channel")
    {
        CFileUtils::GetInfoFromFileName("D2J2134_170129_0317_1.pak", start, serial, channel, mode);
        REQUIRE(channel == 1);

        CFileUtils::GetInfoFromFileName("D2J2134_170129_0317_0.pak", start, serial, channel, mode);
        REQUIRE(channel == 0);
    }

    SECTION("Finds correct mode")
    {
        CFileUtils::GetInfoFromFileName("D2J2134_170129_0317_1.pak", start, serial, channel, mode);
        REQUIRE(mode == MODE_FLUX); // default

        CFileUtils::GetInfoFromFileName("D2J2134_170129_0317_1_wind.pak", start, serial, channel, mode);
        REQUIRE(mode == MODE_WINDSPEED);
    }
}
}