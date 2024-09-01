#pragma once

#include <PPPLib/PPPLib.h>
#include <PPPLib/File/XMLFileReader.h>
#include <PPPLib/Configuration/EvaluationConfiguration.h>
#include <PPPLib/Configuration/DarkCorrectionConfiguration.h>
#include <PPPLib/Configuration/InstrumentCalibrationConfiguration.h>

namespace FileHandler
{
class CEvaluationConfigurationParser : public CXMLFileReader
{
public:
    CEvaluationConfigurationParser(ILogger& logger)
        : CXMLFileReader(logger)
    {
    }

    /** Reads in an evaluation-configuration file.
        In the format specified for the NovacPostProcessingProgram (NPPP)
        @return SUCCESS on sucess */
    RETURN_CODE ReadConfigurationFile(
        const novac::CString& fileName,
        Configuration::CEvaluationConfiguration& settings,
        Configuration::CDarkCorrectionConfiguration& darkSettings,
        Configuration::CInstrumentCalibrationConfiguration& calibrationSettings);

    /** Writes an evaluation configuration file in the NPPP-format
        @return SUCCESS on success */
    RETURN_CODE WriteConfigurationFile(
        const novac::CString& fileName,
        const Configuration::CEvaluationConfiguration& settings,
        const Configuration::CDarkCorrectionConfiguration& darkSettings,
        const Configuration::CInstrumentCalibrationConfiguration& calibrationSettings);

private:

    /** Reads a 'fitWindow' section */
    RETURN_CODE Parse_FitWindow(novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo);

    /** Reads a 'Reference' section */
    RETURN_CODE Parse_Reference(novac::CFitWindow& window);

    /** Reads a 'dark-correction' section */
    RETURN_CODE Parse_DarkCorrection(Configuration::CDarkSettings& dSettings, novac::CDateTime& validFrom, novac::CDateTime& validTo);

    RETURN_CODE Parse_CalibrationSettings(Configuration::CInstrumentCalibrationConfiguration& calibrationSettings);

};
}
