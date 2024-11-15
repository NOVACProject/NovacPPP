#pragma once

#include <PPPLib/Logging.h>
#include <PPPLib/Evaluation/ScanResult.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/ContinuationOfProcessing.h>
#include <PPPLib/PostProcessingStatistics.h>
#include <PPPLib/Evaluation/ExtendedScanResult.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>

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
        novac::ILogger& log,
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
        @param plumeProperties - if not NULL then this will on return be filled
            with the properties of the evaluated scan.
        @return the collected scan result, if the evaluation succeeded and the scan is good enough to calculate a flux from.
        @return nullptr if the scan evaluation failed or the scan is not good enough to calculate a flux from. */
    std::unique_ptr<CExtendedScanResult> EvaluateScan(const std::string& pakFileName, const novac::CString& fitWindowName);


private:
    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    novac::ILogger& m_log;

    Configuration::CNovacPPPConfiguration m_setup;

    Configuration::CUserConfiguration m_userSettings;

    const CContinuationOfProcessing& m_continuation;

    CPostProcessingStatistics& m_processingStats;

    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------

    /** Checks if a scan is goood enough to use for flux calculation.
        @param result an evaluated scan.
        @return true if the measurement should be used to calculate a flux */
    bool IsGoodEnoughToCalculateFlux(novac::LogContext context, std::unique_ptr<CScanResult>& result) const;

};
}