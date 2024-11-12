#pragma once

#include <SpectralEvaluation/DateTime.h>
#include <PPPLib/ContinuationOfProcessing.h>
#include <PPPLib/Geometry/GeometryCalculator.h>
#include <PPPLib/Meteorology/WindDataBase.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/PostProcessingStatistics.h>
#include <PPPLib/Logging.h>

#include <PPPLib/Geometry/PlumeDataBase.h>
#include <PPPLib/Flux/FluxResult.h>
#include <PPPLib/Evaluation/ExtendedScanResult.h>
#include <PPPLib/MFC/CList.h>
#include <PPPLib/MFC/CString.h>

/** The class <b>CPostProcessing</b> is the main class in the NovacPPP
    This is where all the processing takes place (or at least the control
    of the processing).
    See the functions 'DoPostProcessing_Flux' and 'DoPostProcessing_Strat'
*/

namespace novac
{
class CReferenceFile;
}

class CPostProcessing
{
public:
    CPostProcessing(
        novac::ILogger& logger,
        Configuration::CNovacPPPConfiguration setup,
        Configuration::CUserConfiguration userSettings,
        const CContinuationOfProcessing& continuation);

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Performs an post processing of the data in order to extract
        good flux data */
    void DoPostProcessing_Flux();

    /** Performs a post processing of the data in order to create good
        calibrations for the instrument(s) with the purpose of improving the quality
        of the data for other processing modes */
    void DoPostProcessing_InstrumentCalibration();

    /** Performs an post processing of the data in order to extract
        good stratospheric data */
    void DoPostProcessing_Strat();

private:

    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    /** The database of wind-fields to use for the flux calculations */
    Meteorology::CWindDataBase m_windDataBase;

    /** The database of plume-heights to use for the flux calculations */
    Geometry::CPlumeDataBase m_plumeDataBase;

    novac::ILogger& m_log;

    Configuration::CNovacPPPConfiguration m_setup;

    Configuration::CUserConfiguration m_userSettings;

    CContinuationOfProcessing m_continuation;

    // The statistics of the processing itself (number of successfully processed scans, vs number of rejected etc)
    CPostProcessingStatistics m_processingStats;

    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------

    /** Prepares for the post-processing by first making sure that
        all settings found in the configuration are ok and that
        they make sense.
        @throw std::invalid_argument if the settings are not ok and the processing cannot continue */
    void CheckProcessingSettings() const;

    /** Prepares for post-processing by making sure that all settings
        relevant to instrument calibration are ok and make sense.
        @return 0 if all settings are ok otherwise non-zero */
    int CheckInstrumentCalibrationSettings() const;

    /** Prepares for the flux calculations by reading in the relevant
        wind-field file.
        @throw std::invalid_argument if the wind field could not be read properly and the processing cannot continue */
    void ReadWindField(novac::LogContext context);

    /** Prepares for the flux calculation by setting up a reasonable
        set of plume heights. This could also read in a set from file...? */
    void PreparePlumeHeights(novac::LogContext context);

    /** Runs through the supplied list of .pak-files and evaluates
        each one using the setups found in m_setup and m_userSettings.
        @param pakFileList - the list of pak-files to evaluate.
        @param evalLogFiles - will on successful return be filled
            with the path's and filenames of each evaluation log
            file generated and the properties of each scan. */
    void EvaluateScans(
        const std::vector<std::string>& pakFileList,
        std::vector<Evaluation::CExtendedScanResult>& evalLogFiles);

    /** Runs through the supplied list of evaluation - logs and performs
        geometry calculations on the ones which does match. The results
        are returned in the list geometryResults.
        The evaluations are assumed to be sorted in increasing start time of the scan.
        @param evalLogs - list of CExtendedScanResult, each holding the full path and filename
            of an evaluation-log file that should be considered for geometrical
            calculations and the properties of the scan (plume centre position etc)
        @param geometryResults - will on successful return be filled with the
            calculated plume heights and wind-directions. */
    void CalculateGeometries(
        novac::LogContext context,
        std::vector<Evaluation::CExtendedScanResult>& evalLogs,
        std::vector<Geometry::CGeometryResult>& geometryResults);

    /** Writes each of the calculated geometry results to the GeometryLog file */
    void WriteCalculatedGeometriesToFile(novac::LogContext context, const std::vector<Geometry::CGeometryResult>& geometryResults);

    /** Inserts the calculated geometry results into the databases.
        The wind directions will be inserted into m_windDataBase
        The plume altitudes will be inserted into m_plumeDataBase */
    void InsertCalculatedGeometriesIntoDatabase(novac::LogContext context, const std::vector<Geometry::CGeometryResult>& geometryResults);

    /** This calculates the wind speeds from the dual-beam measurements that has been made
        @param evalLogs - list of CExtendedScanResult, each holding the full path and filename
            of an evaluation-log file. Only the measurements containing a
            dual-beam measurement will be considered.
        The plume heights are taken from the database 'm_plumeDataBase' and the
            results are written to the database 'm_windDataBase' */
    void CalculateDualBeamWindSpeeds(novac::LogContext context, const std::vector<Evaluation::CExtendedScanResult>& evalLogs);

    /** Runs through the supplied list of evaluation-results and
        calculates the flux for each scan. The resulting fluxes are written
        to a flux-log file in the output directory.
        @param evalLogs - list of CStrings, each holding the full path and filename
            of an evaluation-log file that should be considered for geometrical
            calculations
        The wind speeds and wind directions will be taken from 'm_windDataBase'
        The plume heigths will be taken from 'm_plumeDataBase'
        */
    void CalculateFluxes(novac::LogContext context, const std::vector<Evaluation::CExtendedScanResult>& evalLogs);


    /** Sorts the evaluation logs in order of increasing time
        (this is mostly done since this speeds up the geometry calculations enormously) */
    void SortEvaluationLogs(std::vector<Evaluation::CExtendedScanResult>& evalLogs);

    /** Writes the calculated fluxes to the flux result file */
    void WriteFluxResult_XML(const std::list<Flux::FluxResult>& calculatedFluxes);
    void WriteFluxResult_Txt(const std::list<Flux::FluxResult>& calculatedFluxes);

    /** Takes care of uploading the result files to the FTP server */
    void UploadResultsToFTP(novac::LogContext context);

    /** Locates evaluation log files in the output directory */
    std::vector<Evaluation::CExtendedScanResult> LocateEvaluationLogFiles(novac::LogContext context, const std::string& directory) const;
};


