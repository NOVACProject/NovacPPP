#include <PPPLib/Meteorology/XMLWindFileReader.h>
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

    static std::string GetWindFieldFile()
    {
        return GetTestDataDirectory() + std::string("wind_pdf_2017_meteoFrance_UTC_novac.wxml");
    }

    TEST_CASE("ReadWindFile gives expected wind profile", "[XMLWindFileReader][File]")
    {
        StdOutLogger logger;
        Meteorology::CWindDataBase resultingDatabase;
        FileHandler::CXMLWindFileReader sut{ logger };

        int returnCode = sut.ReadWindFile(
            GetWindFieldFile(),
            resultingDatabase);

        REQUIRE(returnCode == 0);
        // Expected settings
        {
            REQUIRE(resultingDatabase.m_dataBaseName == "Ruapehu");

            /*
#ifdef _MSC_VER
            REQUIRE("~/Novac/Piton de la Fournaise/OutputFeb2017UTC/" == resultingConfiguration.m_outputDirectory.std_str());
            REQUIRE("~/Novac/Piton de la Fournaise/Temp/" == resultingConfiguration.m_tempDirectory.std_str());
#endif // _MSC_VER

            REQUIRE(PROCESSING_MODE::PROCESSING_MODE_FLUX == resultingConfiguration.m_processingMode);
            REQUIRE(STANDARD_MOLECULE::MOLEC_SO2 == resultingConfiguration.m_molecule);

            REQUIRE(novac::CDateTime(2017, 1, 29, 12, 50, 51) == resultingConfiguration.m_fromDate);
            REQUIRE(novac::CDateTime(2017, 3, 01, 23, 50, 51) == resultingConfiguration.m_toDate);

            REQUIRE("C:\\Temp\\" == resultingConfiguration.m_LocalDirectory.std_str());
            REQUIRE(1 == resultingConfiguration.m_includeSubDirectories_Local);

            REQUIRE("ftp://129.16.35.206/piton_de_la_fournaise/" == resultingConfiguration.m_FTPDirectory.std_str());
            REQUIRE(1 == resultingConfiguration.m_includeSubDirectories_FTP);

            REQUIRE(0.5 == resultingConfiguration.m_calcGeometry_CompletenessLimit);
            REQUIRE(3600 == resultingConfiguration.m_calcGeometryValidTime);
            REQUIRE(100 == resultingConfiguration.m_calcGeometry_MinDistance);

            REQUIRE(7200 == resultingConfiguration.m_dualBeam_ValidTime);
            REQUIRE(10.0 == resultingConfiguration.m_dualBeam_MaxWindSpeedError);

            REQUIRE(1 == resultingConfiguration.m_nFitWindowsToUse);
            REQUIRE(0 == resultingConfiguration.m_mainFitWindow);

            REQUIRE(Configuration::SKY_OPTION::USER_SUPPLIED == resultingConfiguration.sky.skyOption);
            REQUIRE("C:/Temp/Some_sky_spectrum.std" == resultingConfiguration.sky.skySpectrumFile);

    */
        }
    }
}
