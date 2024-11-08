#include <PPPLib/PostProcessingUtils.h>
#include <SpectralEvaluation/Log.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>

#include <PPPLib/File/Filesystem.h>

#include <stdexcept>
#include <iostream>

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

// Endregion Helper methods

TEST_CASE("PrepareEvaluation, reference file not found - throws exception", "[PrepareEvaluation]")
{
    novac::ConsoleLog logger;
    std::string tempDirectory = ".";

    Configuration::CNovacPPPConfiguration setup;
    Configuration::CInstrumentConfiguration instrument;
    instrument.m_serial = "ABC123";

    SECTION("Single instrument, single fit window, missing file")
    {
        // Arrange
        novac::CFitWindow window;
        window.nRef = 1;
        window.ref[0].m_path = "/some/none-existing/reference.txt";
        instrument.m_eval.InsertFitWindow(window, novac::CDateTime::MinValue(), novac::CDateTime::MaxValue());

        setup.m_instrument.push_back(instrument);

        // Act & assert
        REQUIRE_THROWS_AS(novac::PrepareEvaluation(logger, tempDirectory, setup), novac::InvalidReferenceException);
    }

    SECTION("Single instrument, single fit window, second reference missing")
    {
        // Arrange
        novac::CFitWindow window;
        window.nRef = 2;
        window.ref[0].m_path = GetTestDataDirectory() + "2002128M1/2002128M1_SO2_Bogumil_293K.txt"; // exists
        window.ref[1].m_path = "/some/none-existing/reference.txt";
        instrument.m_eval.InsertFitWindow(window, novac::CDateTime::MinValue(), novac::CDateTime::MaxValue());

        setup.m_instrument.push_back(instrument);

        // Act & assert
        REQUIRE_THROWS_AS(novac::PrepareEvaluation(logger, tempDirectory, setup), novac::InvalidReferenceException);
    }

    SECTION("Single instrument, Fraunhofer reference missing")
    {
        // Arrange
        novac::CFitWindow window;
        window.nRef = 2;
        window.ref[0].m_path = GetTestDataDirectory() + "2002128M1/2002128M1_SO2_Bogumil_293K.txt"; // exists
        window.ref[1].m_path = GetTestDataDirectory() + "2002128M1/2002128M1_O3_Voigt_223K.txt"; // exists
        window.fraunhoferRef.m_path = "/some/none-existing/reference.txt";
        instrument.m_eval.InsertFitWindow(window, novac::CDateTime::MinValue(), novac::CDateTime::MaxValue());

        setup.m_instrument.push_back(instrument);

        // Act & assert
        REQUIRE_THROWS_AS(novac::PrepareEvaluation(logger, tempDirectory, setup), novac::InvalidReferenceException);
    }

    SECTION("Single instrument, reference is already filtered.")
    {
        // Arrange
        novac::CFitWindow window;
        window.nRef = 2;
        window.ref[0].m_path = GetTestDataDirectory() + "2002128M1/2002128M1_SO2_Bogumil_293K.txt"; // exists
        window.ref[0].m_isFiltered = true;
        window.ref[1].m_path = GetTestDataDirectory() + "2002128M1/2002128M1_O3_Voigt_223K.txt"; // exists
        window.ref[1].m_isFiltered = true;
        instrument.m_eval.InsertFitWindow(window, novac::CDateTime::MinValue(), novac::CDateTime::MaxValue());

        setup.m_instrument.push_back(instrument);

        // Act & assert
        REQUIRE_THROWS_AS(novac::PrepareEvaluation(logger, tempDirectory, setup), novac::InvalidReferenceException);
    }
}

TEST_CASE("PrepareEvaluation, reads references", "[PrepareEvaluation][Test1]")
{
    novac::ConsoleLog logger;
    std::string tempDirectory = ".";

    Configuration::CNovacPPPConfiguration setup;
    Configuration::CInstrumentConfiguration instrument;
    instrument.m_serial = "ABC123";

    SECTION("Single instrument, Fit includes polynomial and references are not filtered, does not alter references.")
    {
        // Arrange
        novac::CFitWindow window;
        window.fitType = novac::FIT_TYPE::FIT_POLY;
        window.nRef = 2;
        window.ref[0].m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_SO2_Bogumil_293K.txt"; // exists
        window.ref[0].m_isFiltered = false;
        window.ref[1].m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_O3_Voigt_223K.txt"; // exists
        window.ref[1].m_isFiltered = false;
        window.fraunhoferRef.m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_Fraunhofer.txt"; // exists
        instrument.m_eval.InsertFitWindow(window, novac::CDateTime::MinValue(), novac::CDateTime::MaxValue());

        setup.m_instrument.push_back(instrument);

        // Act
        novac::PrepareEvaluation(logger, tempDirectory, setup);

        // Assert
        novac::CFitWindow configuredFitWindow;
        novac::CDateTime ignored1, ignored2;
        setup.m_instrument[0].m_eval.GetFitWindow(0, configuredFitWindow, ignored1, ignored2);

        REQUIRE(configuredFitWindow.nRef == 2);

        REQUIRE(configuredFitWindow.ref[0].m_data != nullptr);
        REQUIRE(configuredFitWindow.ref[0].m_data->m_crossSection.size() == 2048);
        REQUIRE(configuredFitWindow.ref[0].m_data->m_waveLength.size() == 2048);
        REQUIRE(configuredFitWindow.ref[0].m_data->m_crossSection[0] == Approx(4.342150823e-19));
        REQUIRE(configuredFitWindow.ref[0].m_data->m_waveLength[0] == Approx(267.231877000));

        REQUIRE(configuredFitWindow.ref[1].m_data != nullptr);
        REQUIRE(configuredFitWindow.ref[1].m_data->m_crossSection.size() == 2048);
        REQUIRE(configuredFitWindow.ref[1].m_data->m_waveLength.size() == 2048);
        REQUIRE(configuredFitWindow.ref[1].m_data->m_crossSection[0] == Approx(8.742648988e-18));
        REQUIRE(configuredFitWindow.ref[1].m_data->m_waveLength[0] == Approx(267.231877000));

        REQUIRE(configuredFitWindow.fraunhoferRef.m_data != nullptr);
        REQUIRE(configuredFitWindow.fraunhoferRef.m_data->m_crossSection.size() == 2048);
        REQUIRE(configuredFitWindow.fraunhoferRef.m_data->m_waveLength.size() == 2048);
        REQUIRE(configuredFitWindow.fraunhoferRef.m_data->m_crossSection[380] == Approx(std::log(4.192427697e+04)));
        REQUIRE(configuredFitWindow.fraunhoferRef.m_data->m_waveLength[380] == Approx(302.498797860));
    }

    SECTION("Single instrument, Fit is HP500 filtered, references are not filtered, performs filtering on references.")
    {
        // Arrange
        novac::CFitWindow window;
        window.fitType = novac::FIT_TYPE::FIT_HP_DIV;
        window.nRef = 2;
        window.ref[0].m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_SO2_Bogumil_293K.txt"; // exists
        window.ref[0].m_isFiltered = false;
        window.ref[1].m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_O3_Voigt_223K.txt"; // exists
        window.ref[1].m_isFiltered = false;
        window.fraunhoferRef.m_path = GetTestDataDirectory() + "2002128M1/Calibrated/2002128M1_Fraunhofer.txt"; // exists
        instrument.m_eval.InsertFitWindow(window, novac::CDateTime::MinValue(), novac::CDateTime::MaxValue());

        setup.m_instrument.push_back(instrument);

        // Act
        novac::PrepareEvaluation(logger, tempDirectory, setup);

        // Assert
        novac::CFitWindow configuredFitWindow;
        novac::CDateTime ignored1, ignored2;
        setup.m_instrument[0].m_eval.GetFitWindow(0, configuredFitWindow, ignored1, ignored2);

        REQUIRE(configuredFitWindow.nRef == 2);

        REQUIRE(configuredFitWindow.ref[0].m_data != nullptr);
        REQUIRE(configuredFitWindow.ref[0].m_data->m_crossSection.size() == 2048);
        REQUIRE(configuredFitWindow.ref[0].m_data->m_waveLength.size() == 2048);
        REQUIRE(configuredFitWindow.ref[0].m_data->m_crossSection[0] == Approx(2.93e-20).margin(1e-21));
        REQUIRE(configuredFitWindow.ref[0].m_data->m_waveLength[0] == Approx(267.231877000));

        REQUIRE(configuredFitWindow.ref[1].m_data != nullptr);
        REQUIRE(configuredFitWindow.ref[1].m_data->m_crossSection.size() == 2048);
        REQUIRE(configuredFitWindow.ref[1].m_data->m_waveLength.size() == 2048);
        REQUIRE(configuredFitWindow.ref[1].m_data->m_crossSection[0] == Approx(-2.79e-19).margin(1e-20));
        REQUIRE(configuredFitWindow.ref[1].m_data->m_waveLength[0] == Approx(267.231877000));

        REQUIRE(configuredFitWindow.fraunhoferRef.m_data != nullptr);
        REQUIRE(configuredFitWindow.fraunhoferRef.m_data->m_crossSection.size() == 2048);
        REQUIRE(configuredFitWindow.fraunhoferRef.m_data->m_waveLength.size() == 2048);
        REQUIRE(configuredFitWindow.fraunhoferRef.m_data->m_crossSection[380] == Approx(-0.1075877939));
        REQUIRE(configuredFitWindow.fraunhoferRef.m_data->m_waveLength[380] == Approx(302.498797860));
    }
}