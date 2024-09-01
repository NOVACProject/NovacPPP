#include <PPPLib/File/EvaluationConfigurationParser.h>
#include "catch.hpp"
#include "StdOutLogger.h"

namespace novac
{
static std::string GetTestDataDirectory()
{
#ifdef _MSC_VER
    return std::string("../testData/");
#else
    return std::string("testData/");
#endif // _MSC_VER 
}

static std::string GetEvaluationConfigurationFile()
{
    return GetTestDataDirectory() + std::string("I2J8552.exml");
}

static Configuration::CEvaluationConfiguration CreateEvaluationSettings()
{
    Configuration::CEvaluationConfiguration settings;
    settings.m_serial = "ABC80";

    novac::CFitWindow so2Evaluation;
    so2Evaluation.fitLow = 405;
    so2Evaluation.fitHigh = 503;
    so2Evaluation.fitType = novac::FIT_TYPE::FIT_POLY;
    so2Evaluation.name = "SO2";
    so2Evaluation.polyOrder = 3;
    so2Evaluation.nRef = 1;
    so2Evaluation.ref[0].m_path = "C:\\SO2_reference.txt";
    settings.InsertFitWindow(so2Evaluation, novac::CDateTime(1999, 01, 01, 12, 13, 14), novac::CDateTime(2029, 12, 31, 14, 15, 16));


    novac::CFitWindow brOEvaluation;
    brOEvaluation.fitLow = 405;
    brOEvaluation.fitHigh = 503;
    brOEvaluation.fitType = novac::FIT_TYPE::FIT_POLY;
    brOEvaluation.name = "BrO";
    brOEvaluation.polyOrder = 3;
    brOEvaluation.nRef = 3;
    brOEvaluation.ref[0].m_path = "C:\\SO2_reference.txt";
    brOEvaluation.ref[1].m_path = "C:\\Ring_reference.txt";
    brOEvaluation.ref[2].m_path = "C:\\BrO_reference.txt";
    settings.InsertFitWindow(so2Evaluation, novac::CDateTime(1999, 01, 01, 12, 13, 14), novac::CDateTime(2029, 12, 31, 14, 15, 16));

    return settings;
}

static Configuration::CDarkCorrectionConfiguration CreateDarkCorrectionSettings()
{
    Configuration::CDarkCorrectionConfiguration settings;
    Configuration::CDarkSettings darkSettings;
    darkSettings.m_darkSpecOption = Configuration::DARK_SPEC_OPTION::USER_SUPPLIED;
    darkSettings.m_offsetSpec = "C:\\some_user_supplied_dark_spectrum_file.txt";

    settings.InsertDarkCurrentCorrectionSettings(darkSettings, novac::CDateTime(1999, 01, 01, 12, 13, 14), novac::CDateTime(2029, 12, 31, 14, 15, 16));

    return settings;
}

static Configuration::CInstrumentCalibrationConfiguration CreateInstrumentCalibrationSettings()
{
    Configuration::CInstrumentCalibrationConfiguration settings;

    settings.m_initialCalibrationFile = "/mnt/this/is_the_initial/instrument_calibration_file.clb";
    settings.m_instrumentLineshapeFile = "/mnt/this/is_the_initial/instrument_lineshape_file.slf";

    return settings;
}

TEST_CASE("ReadConfigurationFile gives expected configuration", "[EvaluationConfigurationParser][File]")
{
    StdOutLogger logger;
    Configuration::CEvaluationConfiguration resultingEvaluationSettings;
    Configuration::CDarkCorrectionConfiguration resultingDarkSettings;
    Configuration::CInstrumentCalibrationConfiguration resultingCalibrationSettings;
    FileHandler::CEvaluationConfigurationParser sut{ logger };

    RETURN_CODE returnCode = sut.ReadConfigurationFile(
        GetEvaluationConfigurationFile(),
        resultingEvaluationSettings,
        resultingDarkSettings,
        resultingCalibrationSettings);

    REQUIRE(returnCode == RETURN_CODE::SUCCESS);

    // Expected evaluation settings
    {
        REQUIRE("I2J8552" == resultingEvaluationSettings.m_serial);
        REQUIRE(3 == resultingEvaluationSettings.NumberOfFitWindows());

    }

    REQUIRE_NOTHROW(resultingEvaluationSettings.CheckSettings());

    // Expected first evaluation fit window
    {
        novac::CFitWindow window;
        novac::CDateTime validFrom;
        novac::CDateTime validTo;

        int result = resultingEvaluationSettings.GetFitWindow(0, window, validFrom, validTo);

        REQUIRE(result == 0);
        REQUIRE(validFrom == novac::CDateTime(0, 0, 0, 0, 0, 0));
        REQUIRE(validTo == novac::CDateTime(2017, 2, 20, 5, 49, 1));

        REQUIRE(window.nRef == 3);
        REQUIRE(window.ref[0].m_path == "D:\\NovacPostProcessingProgram\\TestRun_2021_12\\OutputFeb2017UTC\\2017.02.20/I2J8552/I2J8552_SO2_Bogumil_293K_170220_0348.txt");
        REQUIRE(window.fraunhoferRef.m_path == "D:\\NovacPostProcessingProgram\\TestRun_2021_12\\OutputFeb2017UTC\\2017.02.20/I2J8552/I2J8552_Fraunhofer_170220_0348.txt");
    }

    // Expected second evaluation fit window
    {
        novac::CFitWindow window;
        novac::CDateTime validFrom;
        novac::CDateTime validTo;

        int result = resultingEvaluationSettings.GetFitWindow(1, window, validFrom, validTo);

        REQUIRE(result == 0);
        REQUIRE(validFrom == novac::CDateTime(2017, 2, 20, 5, 49, 1));
        REQUIRE(validTo == novac::CDateTime(2017, 2, 20, 9, 53, 24));

        REQUIRE(window.nRef == 3);
        REQUIRE(window.ref[0].m_path == "D:\\NovacPostProcessingProgram\\TestRun_2021_12\\OutputFeb2017UTC\\2017.02.20/I2J8552/I2J8552_SO2_Bogumil_293K_170220_0749.txt");
        REQUIRE(window.fraunhoferRef.m_path == "D:\\NovacPostProcessingProgram\\TestRun_2021_12\\OutputFeb2017UTC\\2017.02.20/I2J8552/I2J8552_Fraunhofer_170220_0749.txt");
    }

    // Expected third evaluation fit window
    {
        novac::CFitWindow window;
        novac::CDateTime validFrom;
        novac::CDateTime validTo;

        int result = resultingEvaluationSettings.GetFitWindow(2, window, validFrom, validTo);

        REQUIRE(result == 0);
        REQUIRE(validFrom == novac::CDateTime(2017, 2, 20, 9, 53, 24));
        REQUIRE(validTo == novac::CDateTime(9999, 12, 31, 23, 59, 59));

        REQUIRE(window.nRef == 3);
        REQUIRE(window.ref[0].m_path == "D:\\NovacPostProcessingProgram\\TestRun_2021_12\\OutputFeb2017UTC\\2017.02.20/I2J8552/I2J8552_SO2_Bogumil_293K_170220_1157.txt");
        REQUIRE(window.fraunhoferRef.m_path == "D:\\NovacPostProcessingProgram\\TestRun_2021_12\\OutputFeb2017UTC\\2017.02.20/I2J8552/I2J8552_Fraunhofer_170220_1157.txt");
    }

    // Expected dark settings
    {
        REQUIRE(1 == resultingDarkSettings.GetSettingsNum());
        novac::CDateTime requestTime{ 2017, 2, 20, 10, 0, 0 };
        Configuration::CDarkSettings darkSettings;

        int result = resultingDarkSettings.GetDarkSettings(darkSettings, requestTime);

        REQUIRE(result == 0);
        REQUIRE(darkSettings.m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MEASURED_IN_SCAN);
    }

    // Expected calibration settings
    {
        REQUIRE("C:\\NOVAC\\novacP3\\Cross sections\\I2J8552_SolarSpec.xs" == resultingCalibrationSettings.m_initialCalibrationFile);
    }
}

TEST_CASE("WriteConfigurationFile gives a file which can be read back in again", "[EvaluationConfigurationParser][File]")
{
    StdOutLogger logger;
    FileHandler::CEvaluationConfigurationParser sut{ logger };
    const std::string temporaryEvaluationFileName = GetTestDataDirectory() + "Temporary_evaluationSettings.exml";

    // Create some original settings
    Configuration::CEvaluationConfiguration originalEvaluationSettings = CreateEvaluationSettings();
    Configuration::CDarkCorrectionConfiguration originalDarkSettings = CreateDarkCorrectionSettings();
    Configuration::CInstrumentCalibrationConfiguration originalCalibrationSettings = CreateInstrumentCalibrationSettings();

    // Write these settings to file
    (void)sut.WriteConfigurationFile(
        temporaryEvaluationFileName,
        originalEvaluationSettings,
        originalDarkSettings,
        originalCalibrationSettings);

    // Read the file back again
    Configuration::CEvaluationConfiguration resultingEvaluationSettings;
    Configuration::CDarkCorrectionConfiguration resultingDarkSettings;
    Configuration::CInstrumentCalibrationConfiguration resultingCalibrationSettings;
    (void)sut.ReadConfigurationFile(
        temporaryEvaluationFileName,
        resultingEvaluationSettings,
        resultingDarkSettings,
        resultingCalibrationSettings);

    // Compare the evaluation settings
    {
        REQUIRE(resultingEvaluationSettings.m_serial == originalEvaluationSettings.m_serial);
        REQUIRE(resultingEvaluationSettings.NumberOfFitWindows() == originalEvaluationSettings.NumberOfFitWindows());

        for (int index = 0; index < originalEvaluationSettings.NumberOfFitWindows(); ++index)
        {
            novac::CFitWindow originalWindow;
            novac::CDateTime originalValidFrom, originalValidTo;
            originalEvaluationSettings.GetFitWindow(index, originalWindow, originalValidFrom, originalValidTo);

            novac::CFitWindow resultingWindow;
            novac::CDateTime resultingValidFrom, resultingValidTo;
            resultingEvaluationSettings.GetFitWindow(index, resultingWindow, resultingValidFrom, resultingValidTo);

            REQUIRE(originalValidFrom == resultingValidFrom);
            REQUIRE(originalValidTo == resultingValidTo);

            REQUIRE(originalWindow.fitLow == resultingWindow.fitLow);
            REQUIRE(originalWindow.fitHigh == resultingWindow.fitHigh);
            REQUIRE(originalWindow.fitType == resultingWindow.fitType);
            REQUIRE(originalWindow.polyOrder == resultingWindow.polyOrder);
            REQUIRE(originalWindow.includeIntensitySpacePolyominal == resultingWindow.includeIntensitySpacePolyominal);
            REQUIRE(originalWindow.nRef == resultingWindow.nRef);

            for (int refIndex = 0; refIndex < originalWindow.nRef; ++refIndex)
            {
                REQUIRE(originalWindow.ref[refIndex].m_path == resultingWindow.ref[refIndex].m_path);
            }
        }
    }

    // Compare the dark settings
    {
        REQUIRE(resultingDarkSettings.GetSettingsNum() == originalDarkSettings.GetSettingsNum());

        for (int index = 0; index < originalDarkSettings.GetSettingsNum(); ++index)
        {
            Configuration::CDarkSettings originalsettings;
            novac::CDateTime originalValidFrom, originalValidTo;
            originalDarkSettings.GetDarkSettings(index, originalsettings, originalValidFrom, originalValidTo);

            Configuration::CDarkSettings resultingsettings;
            novac::CDateTime resultingValidFrom, resultingValidTo;
            resultingDarkSettings.GetDarkSettings(index, resultingsettings, resultingValidFrom, resultingValidTo);

            REQUIRE(originalValidFrom == resultingValidFrom);
            REQUIRE(originalValidTo == resultingValidTo);

            REQUIRE(originalsettings.m_darkCurrentOption == resultingsettings.m_darkCurrentOption);
            REQUIRE(originalsettings.m_darkCurrentSpec == resultingsettings.m_darkCurrentSpec);
            REQUIRE(originalsettings.m_darkSpecOption == resultingsettings.m_darkSpecOption);
            REQUIRE(originalsettings.m_offsetOption == resultingsettings.m_offsetOption);
            REQUIRE(originalsettings.m_offsetSpec == resultingsettings.m_offsetSpec);
        }
    }

    // Compare the calibration settings
    REQUIRE(resultingCalibrationSettings.m_initialCalibrationFile == originalCalibrationSettings.m_initialCalibrationFile);
    REQUIRE(resultingCalibrationSettings.m_instrumentLineshapeFile == originalCalibrationSettings.m_instrumentLineshapeFile);
}
}