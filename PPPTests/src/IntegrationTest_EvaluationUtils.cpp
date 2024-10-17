#include <PPPLib/Evaluation/EvaluationUtils.h>
#include "catch.hpp"

static std::string GetTestDataDirectory()
{
#ifdef _MSC_VER
    return std::string("../testData/");
#else
    return std::string("testData/");
#endif // _MSC_VER 
}

void VerifyScanCanBeRead(novac::CScanFileHandler& scan, const std::string filename)
{
    novac::LogContext context(novac::LogContext::FileName, filename);
    bool readOk = scan.CheckScanFile(context, filename); // TODO: Add unit tests to this as well
    if (!readOk)
    {
        printf(scan.m_lastErrorMessage.c_str());
        REQUIRE(readOk == true);
    }
}

TEST_CASE("IsGoodEnoughToEvaluate gives expected configuration", "[EvaluationUtils][IsGoodEnoughToEvaluate]")
{
    novac::CFitWindow fitWindow;
    novac::SpectrometerModel model = novac::CSpectrometerDatabase::SpectrometerModel_AVASPEC();
    Configuration::CInstrumentLocation instrumentLocation;
    Configuration::CUserConfiguration userSettings;
    ReasonForScanRejection reason;
    std::string message;

    SECTION("Scan with saturated sky spectrum in fit region - Returns expected reason")
    {
        const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_0148_0.pak";
        novac::ConsoleLog log;
        novac::CScanFileHandler scan(log);
        VerifyScanCanBeRead(scan, filename);

        fitWindow.fitLow = 464;
        fitWindow.fitHigh = 630;

        bool result = Evaluation::IsGoodEnoughToEvaluate(scan, fitWindow, model, instrumentLocation, userSettings, reason, message);
        REQUIRE_FALSE(result);
        REQUIRE(reason == ReasonForScanRejection::SkySpectrumSaturated);
        REQUIRE(message.length() > 1);
    }

    SECTION("Scan with ok spectrum - Returns true and no failure reason")
    {
        const std::string filename = GetTestDataDirectory() + "2002128M1/2002128M1_230120_0148_0.pak";
        novac::ConsoleLog log;
        novac::CScanFileHandler scan(log);
        VerifyScanCanBeRead(scan, filename);

        // setup a fit region where the scan is not saturated
        fitWindow.fitLow = 320;
        fitWindow.fitHigh = 460;

        bool result = Evaluation::IsGoodEnoughToEvaluate(scan, fitWindow, model, instrumentLocation, userSettings, reason, message);
        REQUIRE(result == true);
        REQUIRE(message.length() == 0);
    }
}
