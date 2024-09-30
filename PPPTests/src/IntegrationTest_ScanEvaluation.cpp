#include <PPPLib/Evaluation/ScanEvaluation.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include "catch.hpp"

static std::string GetTestDataDirectory()
{
#ifdef _MSC_VER
    return std::string("../testData/");
#else
    return std::string("testData/");
#endif // _MSC_VER 
}

// Region Helper methods

extern void VerifyScanCanBeRead(novac::CScanFileHandler& scan, const std::string filename);

static void SetupFitWindow(novac::CFitWindow& window)
{
    window.fitLow = 464;
    window.fitHigh = 630;

    novac::CReferenceFile so2;
    so2.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    so2.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    so2.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_SO2_Bogumil_293K_HP500.txt";
    int success = so2.ReadCrossSectionDataFromFile();
    REQUIRE(success == 0);

    novac::CReferenceFile o3;
    o3.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    o3.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    o3.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_O3_Voigt_223K_HP500.txt";
    success = o3.ReadCrossSectionDataFromFile();
    REQUIRE(success == 0);

    novac::CReferenceFile ring;
    ring.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    ring.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    ring.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_Ring_HP500_PPMM_0.txt";
    success = ring.ReadCrossSectionDataFromFile();
    REQUIRE(success == 0);

    window.ref[0] = so2;
    window.ref[1] = o3;
    window.ref[2] = ring;
    window.nRef = 3;
}

// Endregion Helper methods


TEST_CASE("EvaluateScan, Invalid fit window - throws Exception (Ruahepu, Avantes)", "[ScanEvaluation][EvaluateScan][IntegrationTest][Avantes]")
{
    // Arrange
    const std::string filename = GetTestDataDirectory() + "2002128M1_230120_0148_0.pak";
    novac::ConsoleLog logger;

    novac::CScanFileHandler scan;
    VerifyScanCanBeRead(scan, filename);
    Configuration::CUserConfiguration userSettings;
    const Configuration::CDarkSettings* darkSettings = nullptr;
    novac::SpectrometerModel spectrometerModel = novac::CSpectrometerDatabase::GetInstance().SpectrometerModel_AVASPEC();

    novac::CFitWindow fitWindow;

    novac::CReferenceFile so2;
    so2.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    so2.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    so2.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_SO2_Bogumil_293K_HP500.txt";
    int success = so2.ReadCrossSectionDataFromFile();
    REQUIRE(success == 0);

    Evaluation::CScanEvaluation sut(userSettings, logger);

    SECTION("Fit window has empty range")
    {
        fitWindow.fitLow = 320;
        fitWindow.fitHigh = 320;
        fitWindow.ref[0] = so2;
        fitWindow.nRef = 1;

        // Act & Assert
        REQUIRE_THROWS(sut.EvaluateScan(scan, fitWindow, spectrometerModel, darkSettings));
    }

    SECTION("Fit window has no reference setup")
    {
        fitWindow.fitLow = 320;
        fitWindow.fitHigh = 460;
        fitWindow.ref[0] = so2;
        fitWindow.nRef = 0;

        // Act & Assert
        REQUIRE_THROWS(sut.EvaluateScan(scan, fitWindow, spectrometerModel, darkSettings));
    }

    SECTION("Fit window has same reference multiple times")
    {
        fitWindow.fitLow = 320;
        fitWindow.fitHigh = 460;
        fitWindow.ref[0] = so2;
        fitWindow.ref[1] = so2;
        fitWindow.nRef = 2;

        // Act & Assert
        REQUIRE_THROWS(sut.EvaluateScan(scan, fitWindow, spectrometerModel, darkSettings));
    }
}

TEST_CASE("EvaluateScan, scan with saturated sky spectrum expected result - case 1 (Ruahepu, Avantes)", "[ScanEvaluation][EvaluateScan][IntegrationTest][Avantes]")
{
    // Arrange
    const std::string filename = GetTestDataDirectory() + "2002128M1_230120_0148_0.pak";
    novac::ConsoleLog logger;

    novac::CScanFileHandler scan;
    VerifyScanCanBeRead(scan, filename);
    Configuration::CUserConfiguration userSettings;
    const Configuration::CDarkSettings *darkSettings = nullptr;

    novac::CFitWindow fitWindow;
    SetupFitWindow(fitWindow);

    novac::SpectrometerModel spectrometerModel = novac::CSpectrometerDatabase::GetInstance().SpectrometerModel_AVASPEC();

    Evaluation::CScanEvaluation sut(userSettings, logger);

    // Act
    auto result = sut.EvaluateScan(scan, fitWindow, spectrometerModel, darkSettings);

    // Assert
    REQUIRE(result != nullptr);
}

// TEST_CASE("EvaluateScan, scan with clearly visible plume expected result - case 2 (Ruahepu, Avantes)", "[ScanEvaluation][EvaluateScan][IntegrationTest][Avantes]")
// {
//     const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0.pak";
//     novac::CScanFileHandler scan;
//     VerifyScanCanBeRead(scan, filename);
// 
// 
// 
// }