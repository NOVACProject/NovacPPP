#pragma once

#include <list>
#include <PPPLib/MFC/CString.h>

/// <summary>
/// Listing of the possible reasons for not inluding a scan in the final flux calculation.
/// </summary>
enum class ReasonForScanRejection
{
    SkySpectrumSaturated,
    SkySpectrumDark,
    SkySpectrumTooLongExposureTime,
    CompletenessLow,
    NoPlume,
};

/** The class <b>CPostProcessingStatistics</b> is used to keep
    track of the statistics of the processing. E.g. how many
    scans from a certain instrument are rejected due to different
    problems or how many scans have been processed... */
class CPostProcessingStatistics
{
public:

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Inserts information on a rejected scan from a certain instrument
        into the database. */
    void InsertRejection(const novac::CString& serial, ReasonForScanRejection reason);

    /** Inserts information on a accepted scan from a certain instrument into the database. */
    void InsertAcception(const novac::CString& serial);

    /** Retrieves the number of rejected full scans due to the specified reason */
    unsigned long GetRejectionNum(const novac::CString& serial, ReasonForScanRejection reason);

    /** Retrieves the number of accepted full scans */
    unsigned long GetAcceptionNum(const novac::CString& serial);

    /** Inserts the successful evaluation of a single spectrum into the statistics.
        This also increases the counter on the total amount of time used on
        evaluating spectra */
    void InsertEvaluatedSpectrum(double timeUsed);

    /** Creates a small output file containing the statistical results */
    void WriteStatToFile(const novac::CString& file);

private:
    class CInstrumentStats
    {
    public:
        novac::CString serial = "";
        unsigned long acceptedScans = 0;
        unsigned long noPlumeNum = 0;
        unsigned long lowCompletenessNum = 0;
        unsigned long darkSkySpecNum = 0;
        unsigned long saturatedSkySpecNum = 0;
        unsigned long tooLongExpTime = 0;
    };


    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    /** The statistics for each of the instrument */
    std::list <CInstrumentStats> m_instrumentStats;

    /** The number of spectra evaluated */
    unsigned long nSpectraEvaluated = 0;

    /** The total amount of time spent in the DOAS evaluations */
    double timeSpentOnEvaluations = 0;

};
