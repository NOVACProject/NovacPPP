#pragma once

#include <list>
#include <PPPLib/MFC/CString.h>

/** The class <b>CPostProcessingStatistics</b> is used to keep
    track of the statistics of the processing. E.g. how many
    scans from a certain instrument are rejected due to different
    problems or how many scans have been processed... */

class CPostProcessingStatistics
{
public:
    /** Default constructor */
    CPostProcessingStatistics(void);

    /** Default destructor */
    ~CPostProcessingStatistics(void);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    enum REASON_FOR_REJECTION
    {
        SKY_SPEC_SATURATION,
        SKY_SPEC_DARK,
        SKY_SPEC_TOO_LONG_EXPTIME,
        COMPLETENESS_LOW,
        NO_PLUME
    };

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Inserts information on a rejected scan from a certain instrument
        into the database. */
    void InsertRejection(const novac::CString& serial, const REASON_FOR_REJECTION& reason);

    /** Inserts information on a accepted scan from a certain instrument into the database. */
    void InsertAcception(const novac::CString& serial);

    /** Retrieves the number of rejected full scans due to the specified reason */
    unsigned long GetRejectionNum(const novac::CString& serial, const REASON_FOR_REJECTION& reason);

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
        CInstrumentStats();
        ~CInstrumentStats();
        novac::CString serial;
        unsigned long acceptedScans;
        unsigned long noPlumeNum;
        unsigned long lowCompletenessNum;
        unsigned long darkSkySpecNum;
        unsigned long saturatedSkySpecNum;
        unsigned long tooLongExpTime;
    };


    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    /** The statistics for each of the instrument */
    std::list <CInstrumentStats> m_instrumentStats;

    // ----------- Statistics on the performance of the program 

    /** The number of spectra evaluated */
    unsigned long nSpectraEvaluated;

    /** The total amount of time spent in the DOAS evaluations */
    double timeSpentOnEvaluations;

    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------

};
