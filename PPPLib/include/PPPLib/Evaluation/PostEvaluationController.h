#pragma once

#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/Ratio.h>

#include <PPPLib/Logging.h>
#include <PPPLib/Evaluation/ScanResult.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/ContinuationOfProcessing.h>
#include <PPPLib/PostProcessingStatistics.h>

namespace Evaluation
{
/** <b>CPostEvaluationController</b> is used to to perform the
    evaluation of the scans.

    The main function called (from the outside) is
    <b>EvaluateScan</b> which takes care of checking the spectra
    and calling the help-classes <b>CScanEvaluation</b> to perform the
    actual evaluation. This class takes care of writing the results
    to file and performing other useful things...
*/

class CPostEvaluationController
{
public:
    CPostEvaluationController(
        ILogger& log,
        const Configuration::CNovacPPPConfiguration& setup,
        const Configuration::CUserConfiguration& userSettings,
        const CContinuationOfProcessing& continuation,
        CPostProcessingStatistics& processingStats)
        : m_log(log), m_setup(setup), m_userSettings(userSettings), m_continuation(continuation), m_processingStats(processingStats)
    {
    }

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Evaluates the spectra of one scan and writes the results to file.
        @param pakFileName - the name of the .pak-file that should be evaluated
        @param fitWindowName - the name of the fit-window that should be used in the evaluation
            there can be more than one valid fit-window for each spectrometer at each
            given time. The evaluation will only be performed for the fit-window with the
            given name. If this is empty then the first valid fit-window will be used
        @param txtfileName - if not NULL then this novac::CString will on return be filled
            with the full path and filename of the generated txt-file containing the evaluation
            results
        @param plumeProperties - if not NULL then this will on return be filled
            with the properties of the evaluated scan.
        @return true on success */
    bool EvaluateScan(const novac::CString& pakFileName, const novac::CString& fitWindowName, novac::CString* txtFileName = NULL, novac::CPlumeInScanProperty* plumeProperties = NULL);


private:
    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    ILogger& m_log;

    Configuration::CNovacPPPConfiguration m_setup;

    Configuration::CUserConfiguration m_userSettings;

    const CContinuationOfProcessing& m_continuation;

    CPostProcessingStatistics& m_processingStats;

    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------

    /** Writes the evaluation result to the appropriate log file.
        @param result - a CScanResult holding information about the result
        @param scan - the scan itself, also containing information about the evaluation and the flux.
        @param scanningInstrument - information about the scanning instrument that generated the scan.
        @param txtFileName - if not null, this will on successful writing of the file be filled
            with the full path and filename of the txt - file generated
        @return SUCCESS if operation completed sucessfully. */
    RETURN_CODE WriteEvaluationResult(const std::unique_ptr<CScanResult>& result, const novac::CScanFileHandler* scan, const Configuration::CInstrumentLocation* instrLocation, const novac::CFitWindow* window, Meteorology::CWindField& windField, novac::CString* txtFileName = nullptr);

    /** Writes the evaluation result of one ratio calculation to the appropriate log file.
        @param result - a vector of calculated ratios.
        @param scan - the scan itself, also containing information about the evaluation and the flux.
        @param scanningInstrument - information about the scanning instrument that generated the scan.
        @param txtFileName - if not null, this will on successful writing of the file be filled
            with the full path and filename of the txt - file generated
        @return SUCCESS if operation completed sucessfully. */
    RETURN_CODE WriteRatioResult(const std::vector<novac::Ratio>& result, const novac::CScanFileHandler& scan, const novac::CFitWindow& window);

    /** Appends the evaluation result to the evaluation summary log file.
        @param result - a CScanResult holding information about the result
        @param scan - the scan itself
        @param scanningInstrument - information about the scanning instrument that generated the scan.
        @return SUCCESS if operation completed sucessfully. */
    RETURN_CODE AppendToEvaluationSummaryFile(const std::unique_ptr<CScanResult>& result, const novac::CScanFileHandler* scan, const Configuration::CInstrumentLocation* instrLocation, const novac::CFitWindow* window, Meteorology::CWindField& windField);

    /** Appends the evaluation result to the pak-file summary log file.
        @param result - a CScanResult holding information about the result
        @param scan - the scan itself
        @param scanningInstrument - information about the scanning instrument that generated the scan.
        @return SUCCESS if operation completed sucessfully. */
    RETURN_CODE AppendToPakFileSummaryFile(const std::unique_ptr<CScanResult>& result, const novac::CScanFileHandler* scan, const Configuration::CInstrumentLocation* instrLocation, const novac::CFitWindow* window, Meteorology::CWindField& windField);

    /** This function takes as input parameter an eval-log containing the result of a flux - measurement
        and checks the quality of the measurement.
        @param evalLog - the full path and filename of the flux measurement
        @return 0 - if the measurement should be rejected.
        @return -1 - if the measurement is not a flux measurement.
        */
    int CheckQualityOfFluxMeasurement(std::unique_ptr<CScanResult>& result, const novac::CString& pakFileName) const;

    /** Creates the 'pluem spectrum file' which is a text file containing a list of which spectra are judged to be _in_ the plume
        and which spectra are judged to be _out_ of the plume. Useful for determining plume composition at a later stage */
    void CreatePlumespectrumFile(const std::unique_ptr<CScanResult>& result, const novac::CString& fitWindowName, novac::CScanFileHandler& scan, const novac::SpectrometerModel& spectrometerModel, novac::CPlumeInScanProperty* plumeProperties, int specieIndex);

};
}