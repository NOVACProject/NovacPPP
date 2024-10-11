#include <PPPLib/File/EvaluationLogFileHandler.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Log.h>
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

TEST_CASE("EvaluationLogFileHandler, valid scan file read successfully - Case 1", "[ReadEvaluationLog][FileHandler][IntegrationTest]")
{
    const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0.txt";
    novac::ConsoleLog logger;

    FileHandler::CEvaluationLogFileHandler sut(logger, filename, StandardMolecule::SO2);

    // Act
    auto returnCode = sut.ReadEvaluationLog();

    // Assert
    REQUIRE(returnCode == RETURN_CODE::SUCCESS);

    REQUIRE(sut.m_evaluationLog == filename);
    REQUIRE(sut.m_molecule.name == "SO2");
    REQUIRE(sut.m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG);
    REQUIRE(sut.m_spectrometerModel.modelName == "AVASPEC"); // this should have been retrieved from the data

    REQUIRE(sut.m_specieName.size() == 3);
    REQUIRE(sut.m_specieName[0] == "SO2");
    REQUIRE(sut.m_specieName[1] == "O3");
    REQUIRE(sut.m_specieName[2] == "RING");

    // Verify the spectrum information which is written in the file
    // Notice here that the starttime of the m_specInfo doesn't match the contents of the <scanInformation> section in the file
    // REQUIRE(sut.m_specInfo.m_startTime == novac::CDateTime(2023, 01, 20, 19, 07, 48));
    REQUIRE(sut.m_specInfo.m_gps == novac::CGPSData(-39.277530, 175.608730, 1755.00));
    REQUIRE(sut.m_specInfo.m_compass == 266.0);
    REQUIRE(sut.m_specInfo.m_coneAngle == 60.0);
    REQUIRE(sut.m_specInfo.m_device == "2002128M1");
    REQUIRE(sut.m_specInfo.m_observatory == "geonet");
    REQUIRE(sut.m_specInfo.m_volcano == "ruapehub");
    REQUIRE(sut.m_specInfo.m_channel == 0);
    REQUIRE(sut.m_specInfo.m_interlaceStep == 1);
    REQUIRE(sut.m_specInfo.m_pitch == 0.0);
    REQUIRE(sut.m_specInfo.m_temperature == 27.45F);
    REQUIRE(sut.m_specInfo.m_batteryVoltage == 14.38F);

    // Verify the contents of the scan
    REQUIRE(sut.m_scan.size() == 1);
    REQUIRE(sut.m_scan[0].m_spec.size() == 51);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 7, 48));
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_compass == 266.0F);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_coneAngle == 60.0F);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_device == "2002128M1");
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_channel == 0);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_interlaceStep == 1);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_pitch == 0.0);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_temperature == 27.45F);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_batteryVoltage == 14.38F);

    REQUIRE(sut.m_scan[0].m_spec[0].m_chiSquare == Approx(2.35e-03));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult.size() == 3);
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[0].m_column == Approx(-98.44));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[1].m_column == Approx(503.45));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[2].m_column == Approx(7.47e-03));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[0].m_columnError == Approx(5.64));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[1].m_columnError == Approx(43.42));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[2].m_columnError == Approx(1.66e-03));

    REQUIRE(sut.m_scan[0].m_spec[25].m_chiSquare == Approx(2.34e-03));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult.size() == 3);
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[0].m_column == Approx(-24.03));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[1].m_column == Approx(-216.49));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[2].m_column == Approx(-1.66e-03));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[0].m_columnError == Approx(5.63));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[1].m_columnError == Approx(43.33));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[2].m_columnError == Approx(1.66e-03));

    REQUIRE(sut.m_scan[0].m_spec[50].m_chiSquare == Approx(4.93e-02));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult.size() == 3);
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[0].m_column == Approx(-50.42));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[1].m_column == Approx(-340.78));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[2].m_column == Approx(-5.40e-02));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[0].m_columnError == Approx(25.83));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[1].m_columnError == Approx(198.95));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[2].m_columnError == Approx(7.61e-03));
}

TEST_CASE("EvaluationLogFileHandler, valid scan file from ReEvaluation in NovacProgram read successfully - Case 2", "[ReadEvaluationLog][FileHandler][IntegrationTest]")
{
    const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0_ReEvaluation.txt";
    novac::ConsoleLog logger;

    FileHandler::CEvaluationLogFileHandler sut(logger, filename, StandardMolecule::SO2);

    // Act
    auto returnCode = sut.ReadEvaluationLog();

    // Assert
    REQUIRE(returnCode == RETURN_CODE::SUCCESS);

    REQUIRE(sut.m_evaluationLog == filename);
    REQUIRE(sut.m_molecule.name == "SO2");
    REQUIRE(sut.m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG);
    REQUIRE(sut.m_spectrometerModel.modelName == "AVASPEC"); // this should have been retrieved from the data

    REQUIRE(sut.m_specieName.size() == 4);
    REQUIRE(sut.m_specieName[0] == "SO2");
    REQUIRE(sut.m_specieName[1] == "O3");
    REQUIRE(sut.m_specieName[2] == "RING");
    REQUIRE(sut.m_specieName[3] == "FRAUNHOFERREF");

    // Verify the spectrum information which is written in the file
    // Notice here that the starttime of the m_specInfo doesn't match the contents of the <scanInformation> section in the file
    // REQUIRE(sut.m_specInfo.m_startTime == novac::CDateTime(2023, 01, 20, 19, 07, 48));
    REQUIRE(sut.m_specInfo.m_gps == novac::CGPSData(-39.277530, 175.608730, 1755.00));
    REQUIRE(sut.m_specInfo.m_compass == 266.0);
    REQUIRE(sut.m_specInfo.m_coneAngle == 60.0);
    REQUIRE(sut.m_specInfo.m_device == "2002128M1");
    REQUIRE(sut.m_specInfo.m_channel == 0);
    REQUIRE(sut.m_specInfo.m_interlaceStep == 1);
    REQUIRE(sut.m_specInfo.m_pitch == 0.0);
    REQUIRE(sut.m_specInfo.m_temperature == 27.45F);
    REQUIRE(sut.m_specInfo.m_batteryVoltage == 14.38F);

    // Verify the contents of the scan
    REQUIRE(sut.m_scan.size() == 1);
    REQUIRE(sut.m_scan[0].m_spec.size() == 51);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 7, 48));
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_compass == 266.0F);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_coneAngle == 60.0F);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_device == "2002128M1");
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_channel == 0);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_interlaceStep == 1);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_pitch == 0.0);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_temperature == 27.45F);
    REQUIRE(sut.m_scan[0].m_skySpecInfo.m_batteryVoltage == 14.38F);

    REQUIRE(sut.m_scan[0].m_specInfo[0].m_startTime == novac::CDateTime(2023, 1, 20, 19, 8, 54));
    REQUIRE(sut.m_scan[0].m_specInfo[50].m_startTime == novac::CDateTime(2023, 1, 20, 19, 15, 34));

    REQUIRE(sut.m_scan[0].m_spec[0].m_chiSquare == Approx(2.63e-03));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult.size() == 4);
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[0].m_column == Approx(-1.85e+17));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[1].m_column == Approx(9.34e+17));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[2].m_column == Approx(-2.81e+24));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[0].m_columnError == Approx(1.37e+16));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[1].m_columnError == Approx(1.06e+17));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[2].m_columnError == Approx(5.06e+23));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[0].m_shift == Approx(0.585));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[1].m_shift == Approx(0.585));
    REQUIRE(sut.m_scan[0].m_spec[0].m_referenceResult[2].m_shift == Approx(0.585));

    REQUIRE(sut.m_scan[0].m_spec[25].m_chiSquare == Approx(2.04e-03));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult.size() == 4);
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[0].m_column == Approx(-5.95e+16));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[1].m_column == Approx(-5.13e+17));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[2].m_column == Approx(2.25e+23));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[0].m_columnError == Approx(1.21e+16));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[1].m_columnError == Approx(9.29e+16));
    REQUIRE(sut.m_scan[0].m_spec[25].m_referenceResult[2].m_columnError == Approx(4.45e+23));

    REQUIRE(sut.m_scan[0].m_spec[50].m_chiSquare == Approx(4.07e-02));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult.size() == 4);
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[0].m_column == Approx(-6.69e+15));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[1].m_column == Approx(-1.14e+18));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[2].m_column == Approx(4.54e+24));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[0].m_columnError == Approx(5.41e+16));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[1].m_columnError == Approx(4.15e+17));
    REQUIRE(sut.m_scan[0].m_spec[50].m_referenceResult[2].m_columnError == Approx(1.99e+24));
}