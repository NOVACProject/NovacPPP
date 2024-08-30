#pragma once

#include <PPPLib/File/XMLFileReader.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/PPPLib.h>

namespace FileHandler {

/** The class <b>CSetupFileReader</b> is used to read and write the setup
    files for NovacPPP. The setup files contains the setup data for
    the instruments, such as site-name, lat&long of the position
    and the cone-angle of the instrument. */

class CProcessingFileReader : public CXMLFileReader
{
public:
    CProcessingFileReader(ILogger& logger);

    /** This reads in the contents of a processing.xml-file into the supplied data-structure.
        @param fileName - the name of the file to read from. This must be a .xml file
            in the correct format!
        @param settings - will on successful parsing of the file contain the processing-information
            found in the file.
        @return SUCCESS - if successful.
    */
    RETURN_CODE ReadProcessingFile(const novac::CString& fileName, Configuration::CUserConfiguration& settings);

    /** This takes care of writing the contents of a settings data-structure to file
        Only the part regarding the processing of the data will be written to the file */
    RETURN_CODE WriteProcessingFile(const novac::CString& fileName, const Configuration::CUserConfiguration& settings);

private:
    /** Parses an individual fit-window section */
    void Parse_FitWindow(Configuration::CUserConfiguration& settings);

    /** Parses an individual sky-spectrum section */
    void Parse_SkySpectrum(Configuration::CUserConfiguration& settings);

    /** Parses an individual configuration section */
    void Parse_CalibrationSetting(Configuration::CUserConfiguration& settings);

    /** Parses an individual geometry-calculation section */
    void Parse_GeometryCalc(Configuration::CUserConfiguration& settings);

    /** Parses an individual dual-beam section */
    void Parse_DualBeam(Configuration::CUserConfiguration& settings);

    /** Parses an individual quality-judgement section */
    void Parse_DiscardSettings(Configuration::CUserConfiguration& settings);

    void PrintParameter(FILE* f, int nTabs, const novac::CString& tag, const novac::CString& value);
    void PrintParameter(FILE* f, int nTabs, const novac::CString& tag, const int& value);
    void PrintParameter(FILE* f, int nTabs, const novac::CString& tag, const unsigned int& value);
    void PrintParameter(FILE* f, int nTabs, const novac::CString& tag, const unsigned long& value);
    void PrintParameter(FILE* f, int nTabs, const novac::CString& tag, const double& value);
    void PrintParameter(FILE* f, int nTabs, const novac::CString& tag, const novac::CDateTime& value);

};
}