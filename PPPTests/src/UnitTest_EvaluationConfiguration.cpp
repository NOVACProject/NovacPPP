#include <PPPLib/Configuration/EvaluationConfiguration.h>
#include "catch.hpp"

namespace novac
{

TEST_CASE("EvaluationConfiguration has reasonable initial values", "[EvaluationConfiguration][Configuration]")
{
    SECTION("Initially no fit windows defined")
    {
        Configuration::CEvaluationConfiguration sut;

        REQUIRE(0 == sut.NumberOfFitWindows());
    }

    SECTION("Requesting a fit window if none defined - returns non zero error code")
    {
        Configuration::CEvaluationConfiguration sut;
        novac::CFitWindow window;
        novac::CDateTime validFrom, validTo;

        REQUIRE(0 != sut.GetFitWindow(0, window, validFrom, validTo));
    }
}

TEST_CASE("EvaluationConfiguration can insert and retrieve one window", "[EvaluationConfiguration][Configuration]")
{
    novac::CFitWindow insertedWindow;
    insertedWindow.name = "BrO";
    insertedWindow.fitLow = 543;
    insertedWindow.fitHigh = 654;
    novac::CDateTime insertedValidFrom{ 2021, 2, 28, 11, 23, 59 };
    novac::CDateTime insertedValidTo{ 2021, 3, 15, 25, 59, 59 };

    Configuration::CEvaluationConfiguration sut;

    SECTION("Insert will increment number of windows")
    {
        sut.InsertFitWindow(insertedWindow, insertedValidFrom, insertedValidTo);

        REQUIRE(1 == sut.NumberOfFitWindows());
    }

    SECTION("Requesting inserted window returns same window and time")
    {
        sut.InsertFitWindow(insertedWindow, insertedValidFrom, insertedValidTo);
        novac::CFitWindow returnedWindow;
        novac::CDateTime returnedValidFrom, returnedValidTo;

        REQUIRE(0 == sut.GetFitWindow(0, returnedWindow, returnedValidFrom, returnedValidTo));
        REQUIRE(insertedValidFrom == returnedValidFrom);
        REQUIRE(insertedValidTo == returnedValidTo);
        REQUIRE(insertedWindow.name == returnedWindow.name);
        REQUIRE(insertedWindow.fitLow == returnedWindow.fitLow);
        REQUIRE(insertedWindow.fitHigh == returnedWindow.fitHigh);
    }

    SECTION("Requesting not-inserted window returns error code")
    {
        sut.InsertFitWindow(insertedWindow, insertedValidFrom, insertedValidTo);
        novac::CFitWindow returnedWindow;
        novac::CDateTime returnedValidFrom, returnedValidTo;

        REQUIRE(0 != sut.GetFitWindow(1, returnedWindow, returnedValidFrom, returnedValidTo));
    }
}

TEST_CASE("EvaluationConfiguration can set window at given index and retrieve window", "[EvaluationConfiguration][Configuration]")
{
    novac::CFitWindow insertedWindow;
    insertedWindow.name = "BrO";
    insertedWindow.fitLow = 543;
    insertedWindow.fitHigh = 654;
    novac::CDateTime insertedValidFrom{ 2021, 2, 28, 11, 23, 59 };
    novac::CDateTime insertedValidTo{ 2021, 3, 15, 25, 59, 59 };

    Configuration::CEvaluationConfiguration sut;

    SECTION("Insert at index zero will set number of windows to one")
    {
        sut.SetFitWindow(0, insertedWindow, insertedValidFrom, insertedValidTo);

        REQUIRE(1 == sut.NumberOfFitWindows());
    }

    SECTION("Insert at non-zero position will increment number of windows")
    {
        sut.SetFitWindow(2, insertedWindow, insertedValidFrom, insertedValidTo);

        REQUIRE(3 == sut.NumberOfFitWindows());
    }

    SECTION("Requesting inserted window returns same window and time")
    {
        sut.SetFitWindow(2, insertedWindow, insertedValidFrom, insertedValidTo);
        novac::CFitWindow returnedWindow;
        novac::CDateTime returnedValidFrom, returnedValidTo;

        REQUIRE(0 == sut.GetFitWindow(2, returnedWindow, returnedValidFrom, returnedValidTo));
        REQUIRE(insertedValidFrom == returnedValidFrom);
        REQUIRE(insertedValidTo == returnedValidTo);
        REQUIRE(insertedWindow.name == returnedWindow.name);
        REQUIRE(insertedWindow.fitLow == returnedWindow.fitLow);
        REQUIRE(insertedWindow.fitHigh == returnedWindow.fitHigh);
    }

    SECTION("Requesting not-inserted window returns error code")
    {
        sut.InsertFitWindow(insertedWindow, insertedValidFrom, insertedValidTo);
        novac::CFitWindow returnedWindow;
        novac::CDateTime returnedValidFrom, returnedValidTo;

        REQUIRE(0 != sut.GetFitWindow(1, returnedWindow, returnedValidFrom, returnedValidTo));
    }
}
}
