#pragma once

#include <PPPLib/PPPLib.h>
#include <PPPLib/File/XMLFileReader.h>
#include <PPPLib/Configuration/EvaluationConfiguration.h>
#include <PPPLib/Configuration/DarkCorrectionConfiguration.h>
#include <PPPLib/Configuration/InstrumentCalibrationConfiguration.h>

namespace FileHandler
{

// Exception signalling that there is something wrong with the settings in the EvaluationConfiguration file
class EvaluationConfigurationException : public std::invalid_argument
{
public:
    EvaluationConfigurationException(const std::string& filename, const std::string& msg) :
        std::invalid_argument(std::string("[file=") + filename + "] " + msg)
    {}

    EvaluationConfigurationException(const std::string& file, const std::string& msg, const std::string& token) :
        std::invalid_argument(std::string("[file=") + file + "] (string '" + token + "'): " + msg) {}
};

class CEvaluationConfigurationParser : public CXMLFileReader
{
public:
    CEvaluationConfigurationParser(novac::ILogger& logger)
        : CXMLFileReader(logger)
    {}

    /** Reads in an evaluation-configuration file.
        In the format specified for the NovacPostProcessingProgram (NPPP)
        @return SUCCESS on sucess.
        @throw EvaluationConfigurationException if there are invalid values in the file. */
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
