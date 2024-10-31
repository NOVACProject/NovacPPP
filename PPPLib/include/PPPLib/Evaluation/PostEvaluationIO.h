#pragma once

#include <memory>

#include <PPPLib/PPPLib.h>
#include <PPPLib/Evaluation/ScanResult.h>
#include <PPPLib/Configuration/InstrumentConfiguration.h>
#include <SpectralEvaluation/Log.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

namespace Evaluation
{

/// <summary>
/// Helper class for misc IO operations for the PostEvaluationController.
/// </summary>
class PostEvaluationIO
{
public:

    /** Gets the filename under which the scan-file should be stored.
        @return true if a filename is found. */
    static bool GetArchivingfileName(
        novac::ILogger& log,
        novac::CString& pakFile,
        novac::CString& txtFile,
        const novac::CString& fitWindowName,
        const novac::CString& temporaryScanFile,
        std::string outputDirectory,
        novac::MeasurementMode mode);

    /** Writes the evaluation result to the appropriate log file.
        @param result - a CScanResult holding information about the result
        @param scan - the scan itself, also containing information about the evaluation and the flux.
        @param scanningInstrument - information about the scanning instrument that generated the scan.
        @param txtFileName - if not null, this will on successful writing of the file be filled
            with the full path and filename of the txt - file generated
        @return SUCCESS if operation completed sucessfully. */
    static RETURN_CODE WriteEvaluationResult(
        novac::ILogger& log,
        novac::LogContext context,
        const std::string& outputDirectory,
        novac::SpectrometerModel spectrometerModel,
        const std::unique_ptr<CScanResult>& result,
        const novac::CScanFileHandler* scan,
        const Configuration::CInstrumentLocation* instrLocation,
        const novac::CFitWindow* window,
        Meteorology::WindField& windField,
        novac::CString* txtFileName = nullptr);

    /** Creates the 'pluem spectrum file' which is a text file containing a list of which spectra are judged to be _in_ the plume
        and which spectra are judged to be _out_ of the plume. Useful for determining plume composition at a later stage */
    static void CreatePlumespectrumFile(
        novac::ILogger& log,
        novac::LogContext context,
        const std::string& outputDirectory,
        const std::unique_ptr<CScanResult>& result,
        const novac::CString& fitWindowName,
        novac::CScanFileHandler& scan,
        const novac::SpectrometerModel& spectrometerModel,
        novac::CPlumeInScanProperty* plumeProperties,
        int specieIndex);

    /** Appends the evaluation result to the evaluation summary log file.
        @param result - a CScanResult holding information about the result
        @param scan - the scan itself
        @param scanningInstrument - information about the scanning instrument that generated the scan.
        @return SUCCESS if operation completed sucessfully. */
    static RETURN_CODE AppendToEvaluationSummaryFile(
        const std::string& outputDirectory,
        const std::unique_ptr<CScanResult>& result,
        const novac::CScanFileHandler* scan,
        const Configuration::CInstrumentLocation* instrLocation,
        const novac::CFitWindow* window,
        Meteorology::WindField& windField);

    /** Appends the evaluation result to the pak-file summary log file.
        @param result - a CScanResult holding information about the result
        @param scan - the scan itself
        @param scanningInstrument - information about the scanning instrument that generated the scan.
        @return SUCCESS if operation completed sucessfully. */
    static RETURN_CODE AppendToPakFileSummaryFile(
        const std::string& outputDirectory,
        const std::unique_ptr<CScanResult>& result,
        const novac::CScanFileHandler* scan,
        const Configuration::CInstrumentLocation* instrLocation,
        const novac::CFitWindow* window,
        Meteorology::WindField& windField);
};
}
