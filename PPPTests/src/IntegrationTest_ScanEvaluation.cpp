#include <PPPLib/Evaluation/ScanEvaluation.h>
#include <PPPLib/PostProcessingUtils.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/VectorUtils.h>
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
    so2.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_SO2_Bogumil_293K.txt";

    novac::CReferenceFile o3;
    o3.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    o3.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    o3.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_O3_Voigt_223K.txt";

    novac::CReferenceFile ring;
    ring.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    ring.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    ring.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_Ring_HR.txt";

    window.ref[0] = so2;
    window.ref[1] = o3;
    window.ref[2] = ring;
    window.nRef = 3;
}

static void SetupFitWindowWithCalibratedReferences(novac::CFitWindow& window)
{
    window.fitLow = 464;
    window.fitHigh = 630;

    novac::CReferenceFile so2;
    so2.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    so2.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    so2.m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_SO2_Bogumil_293K.txt";

    novac::CReferenceFile o3;
    o3.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    o3.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    o3.m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_O3_Voigt_223K.txt";

    novac::CReferenceFile ring;
    ring.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    ring.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    ring.m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_Ring.txt";

    novac::CReferenceFile fraunhofer;
    fraunhofer.m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_Fraunhofer.txt";

    window.ref[0] = so2;
    window.ref[1] = o3;
    window.ref[2] = ring;
    window.nRef = 3;
    window.fraunhoferRef = fraunhofer;
}

// Endregion Helper methods

TEST_CASE("EvaluateScan, Invalid fit window - throws Exception", "[ScanEvaluation][EvaluateScan][IntegrationTest][Avantes]")
{
    // Arrange
    const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_0148_0.pak";
    novac::ConsoleLog logger;
    novac::CScanFileHandler scan(logger);
    VerifyScanCanBeRead(scan, filename);
    Configuration::CUserConfiguration userSettings;
    const Configuration::CDarkSettings* darkSettings = nullptr;
    novac::SpectrometerModel spectrometerModel = novac::CSpectrometerDatabase::GetInstance().SpectrometerModel_AVASPEC();

    novac::LogContext context;
    context = context.With(novac::LogContext::FileName, novac::GetFileName(filename));
    context = context.With(novac::LogContext::DeviceModel, spectrometerModel.modelName);

    novac::directorySetup setup;
    setup.executableDirectory = ".";
    setup.tempDirectory = ".";

    novac::CFitWindow fitWindow;

    novac::CReferenceFile so2;
    so2.SetShift(novac::SHIFT_TYPE::SHIFT_FIX, 0.0);
    so2.SetSqueeze(novac::SHIFT_TYPE::SHIFT_FIX, 1.0);
    so2.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_SO2_Bogumil_293K_HP500.txt";

    Evaluation::CScanEvaluation sut(userSettings, logger);

    SECTION("Fit window has empty range")
    {
        fitWindow.fitLow = 320;
        fitWindow.fitHigh = 320;
        fitWindow.ref[0] = so2;
        fitWindow.nRef = 1;

        // Act & Assert
        REQUIRE_THROWS(sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings));
    }

    SECTION("Fit window has no reference setup")
    {
        fitWindow.fitLow = 320;
        fitWindow.fitHigh = 460;
        fitWindow.ref[0] = so2;
        fitWindow.nRef = 0;

        // Act & Assert
        REQUIRE_THROWS(sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings));
    }

    SECTION("Fit window has same reference multiple times")
    {
        fitWindow.fitLow = 320;
        fitWindow.fitHigh = 460;
        fitWindow.ref[0] = so2;
        fitWindow.ref[1] = so2;
        fitWindow.nRef = 2;

        // Act & Assert
        REQUIRE_THROWS(sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings));
    }
}

TEST_CASE("EvaluateScan, scan with saturated sky spectrum expected result ", "[ScanEvaluation][EvaluateScan][IntegrationTest][Avantes][2002128M1_230120_0148_0]")
{
    // Arrange
    const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_0148_0.pak";
    novac::ConsoleLog logger;
    novac::CScanFileHandler scan(logger);
    VerifyScanCanBeRead(scan, filename);
    Configuration::CUserConfiguration userSettings;
    const Configuration::CDarkSettings* darkSettings = nullptr;

    const novac::SpectrometerModel spectrometerModel = novac::CSpectrometerDatabase::GetInstance().SpectrometerModel_AVASPEC();

    novac::directorySetup setup;
    setup.executableDirectory = ".";
    setup.tempDirectory = ".";

    novac::LogContext context;
    context = context.With(novac::LogContext::FileName, novac::GetFileName(filename));
    context = context.With(novac::LogContext::DeviceModel, spectrometerModel.modelName);

    novac::CFitWindow fitWindow;
    fitWindow.fitType = novac::FIT_TYPE::FIT_HP_DIV; // the references are HP500
    SetupFitWindow(fitWindow);
    PrepareFitWindow(logger, context, "2002128M1", fitWindow, setup);

    Evaluation::CScanEvaluation sut(userSettings, logger);

    // Act
    auto result = sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings);

    // Assert, the sky spectrum is bad hence the entire scan is bad. Skip
    REQUIRE(result == nullptr);
}

TEST_CASE("EvaluateScan, scan with visible plume)", "[ScanEvaluation][EvaluateScan][IntegrationTest][Avantes][2002128M1_230120_1907_0]")
{
    // Arrange
    const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0.pak";

    novac::ConsoleLog logger;
    novac::CScanFileHandler scan(logger);
    VerifyScanCanBeRead(scan, filename);
    Configuration::CUserConfiguration userSettings;
    const Configuration::CDarkSettings* darkSettings = nullptr;

    const novac::SpectrometerModel spectrometerModel = novac::CSpectrometerDatabase::GetInstance().SpectrometerModel_AVASPEC();

    novac::directorySetup setup;
    setup.executableDirectory = ".";
    setup.tempDirectory = ".";

    novac::LogContext context;
    context = context.With(novac::LogContext::FileName, novac::GetFileName(filename));
    context = context.With(novac::LogContext::DeviceModel, spectrometerModel.modelName);

    SECTION("Default settings")
    {
        novac::CFitWindow fitWindow;
        fitWindow.fitType = novac::FIT_TYPE::FIT_HP_DIV;
        SetupFitWindow(fitWindow);
        PrepareFitWindow(logger, context, "2002128M1", fitWindow, setup);

        Evaluation::CScanEvaluation sut(userSettings, logger);

        // Act
        auto result = sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(44 == result->GetEvaluatedNum());
        REQUIRE(-90.0 == result->GetScanAngle(0));
        REQUIRE(82.0 == result->GetScanAngle(43));

        REQUIRE(Approx(-1.593e17).margin(1e16) == result->GetColumn(0, 0));
        REQUIRE(0.0 == result->GetShift(0, 0));

        // All the evaluations should be reasonably good
        for (size_t specIdx = 0; specIdx < result->GetEvaluatedNum(); ++specIdx)
        {
            REQUIRE(result->GetChiSquare(specIdx) < 8e-2);
        }

        const std::vector<double> columns = novac::GetColumns(*result, 0);
        REQUIRE(Approx(1.3e17).margin(2e16) == Max(columns));
        REQUIRE(Approx(-1.8e17).margin(2e16) == Min(columns));

        // the sky spectrum info should be set
        const auto& skySpecInfo = result->GetSkySpectrumInfo();
        REQUIRE(skySpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 7, 48, 870));

        // the dark spectrum info should be set
        const auto& darkSpecInfo = result->GetDarkSpectrumInfo();
        REQUIRE(darkSpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 8, 29, 240));
    }

    SECTION("Default setup but FIT_HP_SUB")
    {
        novac::CFitWindow fitWindow;
        fitWindow.fitType = novac::FIT_TYPE::FIT_HP_SUB;
        SetupFitWindow(fitWindow);
        PrepareFitWindow(logger, context, "2002128M1", fitWindow, setup);

        Evaluation::CScanEvaluation sut(userSettings, logger);

        // Act
        auto result = sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(44 == result->GetEvaluatedNum());
        REQUIRE(-90.0 == result->GetScanAngle(0));
        REQUIRE(82.0 == result->GetScanAngle(43));

        REQUIRE(Approx(-1.593e17).margin(1e16) == result->GetColumn(0, 0));
        REQUIRE(0.0 == result->GetShift(0, 0));

        // All the evaluations should be reasonably good
        for (size_t specIdx = 0; specIdx < result->GetEvaluatedNum(); ++specIdx)
        {
            REQUIRE(result->GetChiSquare(specIdx) < 8e-2);
        }

        const std::vector<double> columns = novac::GetColumns(*result, 0);
        REQUIRE(Approx(1.4e17).margin(2e16) == Max(columns));
        REQUIRE(Approx(-1.8e17).margin(2e16) == Min(columns));

        // the sky spectrum info should be set
        const auto& skySpecInfo = result->GetSkySpectrumInfo();
        REQUIRE(skySpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 7, 48, 870));

        // the dark spectrum info should be set
        const auto& darkSpecInfo = result->GetDarkSpectrumInfo();
        REQUIRE(darkSpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 8, 29, 240));
    }

    SECTION("Find optimum shift of references (there is a shift between the references and the spectra)")
    {
        novac::CFitWindow fitWindow;
        fitWindow.fitType = novac::FIT_TYPE::FIT_HP_DIV;
        fitWindow.findOptimalShift = 1;
        SetupFitWindow(fitWindow);
        PrepareFitWindow(logger, context, "2002128M1", fitWindow, setup);

        Evaluation::CScanEvaluation sut(userSettings, logger);

        // Act
        auto result = sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(44 == result->GetEvaluatedNum());
        REQUIRE(-90.0 == result->GetScanAngle(0));
        REQUIRE(82.0 == result->GetScanAngle(43));

        REQUIRE(Approx(-1.593e17).margin(1e16) == result->GetColumn(0, 0));
        REQUIRE(Approx(0.07769).margin(0.001) == result->GetShift(0, 0));
        REQUIRE(Approx(1.00) == result->GetSqueeze(0, 0));

        // All the evaluations should be reasonably good and have the same shift/squeeze
        for (size_t specIdx = 0; specIdx < result->GetEvaluatedNum(); ++specIdx)
        {
            REQUIRE(result->GetChiSquare(specIdx) < 8e-2);

            REQUIRE(3 == result->GetSpecieNum(specIdx));
            for (size_t specieIdx = 0; specieIdx < result->GetSpecieNum(specIdx); ++specieIdx)
            {
                REQUIRE(Approx(0.07769).margin(0.001) == result->GetShift(specIdx, specieIdx));
                REQUIRE(Approx(1.00) == result->GetSqueeze(specIdx, specieIdx));
            }
        }

        const std::vector<double> columns = novac::GetColumns(*result, 0);
        REQUIRE(Approx(1.4e17).margin(2e16) == Max(columns));
        REQUIRE(Approx(-1.8e17).margin(2e16) == Min(columns));
    }

    SECTION("Find optimum shift from Fraunhofer reference")
    {
        novac::CFitWindow fitWindow;
        fitWindow.fitType = novac::FIT_TYPE::FIT_HP_SUB; // This is actually the setup used in the NZ setup
        fitWindow.findOptimalShift = 1;
        SetupFitWindow(fitWindow);

        novac::CReferenceFile fraunhofer;
        fraunhofer.m_path = GetTestDataDirectory() + "2002128M1/2002128M1_Solar.txt";
        fitWindow.fraunhoferRef = fraunhofer;

        PrepareFitWindow(logger, context, "2002128M1", fitWindow, setup);

        const double expectedShift = 1.460;

        Evaluation::CScanEvaluation sut(userSettings, logger);

        // Act
        auto result = sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(44 == result->GetEvaluatedNum());
        REQUIRE(-90.0 == result->GetScanAngle(0));
        REQUIRE(82.0 == result->GetScanAngle(43));

        REQUIRE(Approx(expectedShift).margin(0.001) == result->GetShift(0, 0));
        REQUIRE(Approx(1.00) == result->GetSqueeze(0, 0));

        // TODO: This does not agree with the realtime results
        REQUIRE(Approx(-2.293e+17).margin(1e16) == result->GetColumn(0, 0));

        // All the evaluations should be reasonably good and have the same shift/squeeze
        for (size_t specIdx = 0; specIdx < result->GetEvaluatedNum(); ++specIdx)
        {
            REQUIRE(result->GetChiSquare(specIdx) < 8e-2);

            REQUIRE(4 == result->GetSpecieNum(specIdx));
            for (size_t specieIdx = 0; specieIdx < 3; ++specieIdx)
            {
                REQUIRE(Approx(expectedShift).margin(0.001) == result->GetShift(specIdx, specieIdx));
                REQUIRE(Approx(1.00) == result->GetSqueeze(specIdx, specieIdx));
            }
        }

        const std::vector<double> columns = novac::GetColumns(*result, 0);
        REQUIRE(Approx(1.65e17).margin(2e16) == Max(columns));
        REQUIRE(Approx(-2.36e17).margin(2e16) == Min(columns));
    }

    SECTION("Lower minimum saturation ratio evaluates all spectra")
    {
        userSettings.m_minimumSaturationInFitRegion = 0.01; // less than the default value of 5%

        novac::CFitWindow fitWindow;
        fitWindow.fitType = novac::FIT_TYPE::FIT_HP_DIV; // the references are HP500
        SetupFitWindow(fitWindow);
        PrepareFitWindow(logger, context, "2002128M1", fitWindow, setup);

        Evaluation::CScanEvaluation sut(userSettings, logger);

        // Act
        auto result = sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(51 == result->GetEvaluatedNum()); // all spectra evaluated
        REQUIRE(-90.0 == result->GetScanAngle(0));
        REQUIRE(90.0 == result->GetScanAngle(50));
    }
}

TEST_CASE("EvaluateScan, scan with visible plume and calibrated references", "[ScanEvaluation][EvaluateScan][IntegrationTest][Avantes][2002128M1_230120_1907_0]")
{
    // Arrange
    const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_1907_0.pak";

    novac::ConsoleLog logger;

    novac::CScanFileHandler scan(logger);
    VerifyScanCanBeRead(scan, filename);
    Configuration::CUserConfiguration userSettings;
    const Configuration::CDarkSettings* darkSettings = nullptr;

    const novac::SpectrometerModel spectrometerModel = novac::CSpectrometerDatabase::GetInstance().SpectrometerModel_AVASPEC();

    novac::directorySetup setup;
    setup.executableDirectory = ".";
    setup.tempDirectory = ".";

    novac::LogContext context;
    context = context.With(novac::LogContext::FileName, novac::GetFileName(filename));
    context = context.With(novac::LogContext::DeviceModel, spectrometerModel.modelName);

    SECTION("Calibrated references with Polynomial fit and Fraunhofer Reference")
    {
        novac::CFitWindow fitWindow;
        fitWindow.fitType = novac::FIT_TYPE::FIT_POLY;
        SetupFitWindowWithCalibratedReferences(fitWindow);
        PrepareFitWindow(logger, context, "2002128M1", fitWindow, setup);

        Evaluation::CScanEvaluation sut(userSettings, logger);

        // Act
        auto result = sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(44 == result->GetEvaluatedNum());
        REQUIRE(-90.0 == result->GetScanAngle(0));
        REQUIRE(82.0 == result->GetScanAngle(43));

        // A shift has been applied to the DOAS fit.
        REQUIRE(Approx(-0.446).margin(0.01) == result->GetShift(0, 0));

        REQUIRE(Approx(-2.238e17).margin(1e15) == result->GetColumn(0, 0));

        // the sky spectrum info should be set
        auto skySpecInfo = result->GetSkySpectrumInfo();
        REQUIRE(skySpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 7, 48, 870));

        // the dark spectrum info should be set
        auto darkSpecInfo = result->GetDarkSpectrumInfo();
        REQUIRE(darkSpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 8, 29, 240));
    }

    SECTION("Calibrated references with HP500 fit and Fraunhofer Reference")
    {
        novac::CFitWindow fitWindow;
        fitWindow.fitType = novac::FIT_TYPE::FIT_HP_DIV;
        SetupFitWindowWithCalibratedReferences(fitWindow);
        PrepareFitWindow(logger, context, "2002128M1", fitWindow, setup);

        Evaluation::CScanEvaluation sut(userSettings, logger);

        // Act
        auto result = sut.EvaluateScan(context, scan, fitWindow, spectrometerModel, darkSettings);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(44 == result->GetEvaluatedNum());
        REQUIRE(-90.0 == result->GetScanAngle(0));
        REQUIRE(82.0 == result->GetScanAngle(43));

        // A shift has been applied to the DOAS fit.
        REQUIRE(Approx(-0.43).margin(0.01) == result->GetShift(0, 0));

        REQUIRE(Approx(-89.43 * 2.5e15).margin(2.5e15) == result->GetColumn(0, 0));

        // the sky spectrum info should be set
        auto& skySpecInfo = result->GetSkySpectrumInfo();
        REQUIRE(skySpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 7, 48, 870));

        // the dark spectrum info should be set
        auto& darkSpecInfo = result->GetDarkSpectrumInfo();
        REQUIRE(darkSpecInfo.m_startTime == novac::CDateTime(2023, 1, 20, 19, 8, 29, 240));
    }
}

