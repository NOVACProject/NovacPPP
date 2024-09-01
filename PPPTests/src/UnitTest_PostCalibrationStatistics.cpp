#include <PPPLib/Calibration/PostCalibrationStatistics.h>
#include "catch.hpp"

namespace novac
{
static void CreateCalibrations(CPostCalibrationStatistics& sut, const SpectrometerId& calibratedInstrument, const std::vector<CDateTime>& calibrationsPerformedAt)
{
    for (const auto& timestamp : calibrationsPerformedAt)
    {
        std::vector<CReferenceFile> calibratedReferences;
        std::vector<std::string> specie{ "SO2", "O3", "Ring" };
        for (int idx = 0; idx < 3; ++idx)
        {
            CReferenceFile ref;
            ref.m_path = "C:\\Temp\\" + specie[idx] + ".xs";
            calibratedReferences.push_back(ref);
        }

        sut.RememberCalibrationPerformed(calibratedInstrument, timestamp, calibratedReferences);
    }
}

TEST_CASE("PostCalibrationStatistics has reasonable initial values", "[PostCalibrationStatistics][Configuration]")
{
    SECTION("Initially no calibrations")
    {
        CPostCalibrationStatistics sut;
        SpectrometerId instrument{ "ABC80", 0 };

        REQUIRE(0 == sut.GetNumberOfCalibrationsPerformedFor(instrument));
    }

    SECTION("GetCalibration throws invalid argument if no calibration performed at all")
    {
        CPostCalibrationStatistics sut;
        SpectrometerId instrument{ "ABC80", 0 };
        int requestedIndex = 0;
        CDateTime validFrom, validTo;
        std::vector<CReferenceFile> references;

        REQUIRE_THROWS(sut.GetCalibration(instrument, requestedIndex, validFrom, validTo, references));
    }

    SECTION("GetCalibration throws invalid argument if no calibration performed for given instrument")
    {
        CPostCalibrationStatistics sut;
        SpectrometerId requestedInstrument{ "ABC80", 0 };
        int requestedIndex = 0;
        CDateTime returnedValidFrom, returnedValidTo;
        std::vector<CReferenceFile> returnedReferences;

        SpectrometerId calibratedInstrument{ "ABC80", 1 };
        CDateTime calibratedTime{ 2021, 2, 21, 18, 19, 20 };
        std::vector<CReferenceFile> calibratedReferences;
        sut.RememberCalibrationPerformed(calibratedInstrument, calibratedTime, calibratedReferences);

        REQUIRE_THROWS(sut.GetCalibration(requestedInstrument, requestedIndex, returnedValidFrom, returnedValidTo, returnedReferences));
    }
}

TEST_CASE("PostCalibrationStatistics remembers calibrations performed", "[PostCalibrationStatistics][Configuration]")
{
    CPostCalibrationStatistics sut;
    SpectrometerId calibratedInstrument{ "ABC80", 0 };
    std::vector<CDateTime> calibrationsPerformedAt
    {
        CDateTime { 2021, 2, 21, 18, 0, 0 },
        CDateTime { 2021, 2, 21, 9, 0, 0 },
        CDateTime { 2021, 2, 21, 15, 0, 0 },
        CDateTime { 2021, 2, 21, 12, 0, 0 },
    };
    CreateCalibrations(sut, calibratedInstrument, calibrationsPerformedAt);

    SECTION("Correct number of calibrations stored")
    {
        REQUIRE(4 == sut.GetNumberOfCalibrationsPerformedFor(calibratedInstrument));
    }

    SECTION("Retrieves calibrations in order of increasing time")
    {
        CDateTime firstValidFrom, firstValidTo;
        std::vector<CReferenceFile> references;
        sut.GetCalibration(calibratedInstrument, 0, firstValidFrom, firstValidTo, references);
        REQUIRE(firstValidFrom < firstValidTo);

        CDateTime secondValidFrom, secondValidTo;
        sut.GetCalibration(calibratedInstrument, 1, secondValidFrom, secondValidTo, references);
        REQUIRE(secondValidFrom < secondValidTo);
        REQUIRE(secondValidFrom >= firstValidTo);

        CDateTime thirdValidFrom, thirdValidTo;
        sut.GetCalibration(calibratedInstrument, 2, thirdValidFrom, thirdValidTo, references);
        REQUIRE(thirdValidFrom < thirdValidTo);
        REQUIRE(thirdValidFrom >= secondValidTo);

        CDateTime fourthValidFrom, fourthValidTo;
        sut.GetCalibration(calibratedInstrument, 3, fourthValidFrom, fourthValidTo, references);
        REQUIRE(fourthValidFrom < fourthValidTo);
        REQUIRE(fourthValidFrom >= thirdValidTo);
    }

    SECTION("Retrieved calibration times are between actual calibrations performed")
    {
        std::sort(begin(calibrationsPerformedAt), end(calibrationsPerformedAt)); // sort the timestamps in order of increasing time.

        CDateTime validFrom, validTo;
        std::vector<CReferenceFile> references;

        // get the first calibration
        sut.GetCalibration(calibratedInstrument, 0, validFrom, validTo, references);
        REQUIRE(validFrom < calibrationsPerformedAt[0]);
        REQUIRE(validTo > calibrationsPerformedAt[0]);
        REQUIRE(validTo < calibrationsPerformedAt[1]);

        // get the second calibration
        sut.GetCalibration(calibratedInstrument, 1, validFrom, validTo, references);
        REQUIRE(validFrom < calibrationsPerformedAt[1]);
        REQUIRE(validTo > calibrationsPerformedAt[1]);
        REQUIRE(validTo < calibrationsPerformedAt[2]);

        // get the third calibration
        sut.GetCalibration(calibratedInstrument, 2, validFrom, validTo, references);
        REQUIRE(validFrom < calibrationsPerformedAt[2]);
        REQUIRE(validTo > calibrationsPerformedAt[2]);
        REQUIRE(validTo < calibrationsPerformedAt[3]);

        // get the fourth and final calibration
        sut.GetCalibration(calibratedInstrument, 3, validFrom, validTo, references);
        REQUIRE(validFrom < calibrationsPerformedAt[3]);
        REQUIRE(validTo > calibrationsPerformedAt[3]);
    }

    SECTION("Invalid index throws invalid_argument")
    {
        CDateTime validFrom, validTo;
        std::vector<CReferenceFile> references;

        REQUIRE_THROWS(sut.GetCalibration(calibratedInstrument, -1, validFrom, validTo, references));
        REQUIRE_THROWS(sut.GetCalibration(calibratedInstrument, 4, validFrom, validTo, references));
    }
}
}
