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

TEST_CASE("EvaluationLogFileHandler, valid scan file read successfully", "[ReadEvaluationLog][FileHandler][IntegrationTest]")
{
    const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0.txt";
    novac::ConsoleLog logger;

    FileHandler::CEvaluationLogFileHandler sut(logger, filename, STANDARD_MOLECULE::MOLEC_SO2);

    // Act
    auto returnCode = sut.ReadEvaluationLog();

    // Assert
    REQUIRE(returnCode == RETURN_CODE::SUCCESS);

    REQUIRE(sut.m_evaluationLog == filename);
    REQUIRE(sut.m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG);
    REQUIRE(sut.m_molecule.m_name == "SO2");

    REQUIRE(sut.m_specieName.size() == 3);
    REQUIRE(sut.m_specieName[0] == "SO2");
    REQUIRE(sut.m_specieName[1] == "O3");
    REQUIRE(sut.m_specieName[2] == "RING");

    // Verify the contents of the scan
    REQUIRE(sut.m_scan.size() == 1);
    REQUIRE(sut.m_scan[0].m_spec.size() == 51);

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