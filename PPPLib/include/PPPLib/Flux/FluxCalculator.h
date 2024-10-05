#pragma once

#include <PPPLib/Flux/FluxResult.h>
#include <PPPLib/Evaluation/ScanResult.h>
#include <PPPLib/PPPLib.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Meteorology/WindDataBase.h>
#include <PPPLib/Geometry/PlumeHeight.h>

#include <PPPLib/MFC/CString.h>
#include <SpectralEvaluation/Log.h>

namespace Configuration
{
class CNovacPPPConfiguration;
class CUserConfiguration;
}

namespace Flux
{
/** The class <b>CFluxCalculator</b> s used to to perform the
    flux calculation of the scans.

    This class corresponds to CPostEvaluationController when it
    comes to calculating the fluxes.

    The main function called (from the outside) is
    <b>CalculateFlux</b> which takes care of checking the spectra
    and calling the help-functions in <b>Common</b> to perform the
    actual flux calculation. This class takes care of writing the results
    to file and performing other useful things...
*/
class CFluxCalculator
{
public:
    /** Default constructor */
    CFluxCalculator(
        novac::ILogger& log,
        const Configuration::CNovacPPPConfiguration& setup,
        const Configuration::CUserConfiguration& userSettings);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------


    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Calculates the flux from the scan found in the given evaluation log file
        @param evalLogFileName - the name of the .txt-file that contains
            the result of the evaluation.
        @param windDataBase - a database with information about the wind. The
            parameters for the wind will be taken from this database. The function
            fails if no acceptable wind-field could be found
        @param plumeheight - information about the altitude of the plume. This should
            be in meters above sea level.
        @param fluxResult - will on successful calculation of the flux be filled with
            the result of the calculations.
        @return true on success.
      */
    bool CalculateFlux(
        novac::LogContext context,
        const std::string& evalLogFileName,
        const Meteorology::CWindDataBase& windDataBase,
        const Geometry::CPlumeHeight& plumeAltitude,
        CFluxResult& fluxResult);

    /** Calculates the flux using the supplied data.
        Automatically decides which algorithm to use based on the given cone angle.
        (Moved from Common.h where it previously resided in the NovacPPP) */
    static double CalculateFlux(const double* scanAngle, const double* scanAngle2, const double* column, double offset, int nDataPoints, const Meteorology::CWindField& wind, const Geometry::CPlumeHeight& relativePlumeHeight, double compass, INSTRUMENT_TYPE type, double coneAngle = 90.0, double tilt = 0.0);

private:
    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    novac::ILogger& m_log;

    const Configuration::CNovacPPPConfiguration& m_setup;

    const Configuration::CUserConfiguration& m_userSettings;

    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------

    /** Looks in the configuration of the instruments and searches
        for a configured location which is valid for the spectrometer that
        collected the given scan and is also valid at the time when the scan
        was made.
        @return 0 if successful otherwise non-zero
    */
    int GetLocation(novac::LogContext context,
        const novac::CString& serial,
        const novac::CDateTime& startTime,
        Configuration::CInstrumentLocation& instrLocation);

    /** Appends the evaluated flux to the appropriate log file.
        @param scan - the scan itself, also containing information about the evaluation and the flux.
        @return SUCCESS if operation completed sucessfully. */
    RETURN_CODE WriteFluxResult(
        novac::LogContext context,
        const Flux::CFluxResult& fluxResult,
        const Evaluation::CScanResult* result);
};
}