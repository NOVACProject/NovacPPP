#include "stdafx.h"
#include "PostProcessingStatistics.h"
#include "Common/Common.h"
#include <PPPLib/MFC/CCriticalSection.h>
#include <PPPLib/MFC/CSingleLock.h>
#include <PPPLib/File/Filesystem.h>

// Include synchronization classes
// #include <afxmt.h>

novac::CCriticalSection g_processingStatCritSect; // synchronization access to the processing statistics

CPostProcessingStatistics::CInstrumentStats::CInstrumentStats()
{
    // the serial of this instrument
    serial.Format("");

    // accepted scans
    acceptedScans = 0;

    // rejected scans
    noPlumeNum = 0;
    lowCompletenessNum = 0;
    darkSkySpecNum = 0;
    saturatedSkySpecNum = 0;
    tooLongExpTime = 0;
}

CPostProcessingStatistics::CInstrumentStats::~CInstrumentStats()
{

}


CPostProcessingStatistics::CPostProcessingStatistics(void)
{
    // performance statistics
    nSpectraEvaluated = 0;
    timeSpentOnEvaluations = 0.0;
}

CPostProcessingStatistics::~CPostProcessingStatistics(void)
{
}

/** Inserts information on a rejected scan from a certain instrument
    into the database. */
void CPostProcessingStatistics::InsertRejection(const novac::CString& serial, const REASON_FOR_REJECTION& reason)
{

    novac::CSingleLock singleLock(&g_processingStatCritSect);
    singleLock.Lock();
    if (singleLock.IsLocked())
    {

        // look for the correct instrument
        std::list <CInstrumentStats>::const_iterator pos = m_instrumentStats.begin();
        while (pos != m_instrumentStats.end())
        {
            CInstrumentStats& stat = (CInstrumentStats&)*(pos++);

            // this is the instrument. Insert the new data
            if (Equals(stat.serial, serial))
            {
                switch (reason)
                {
                case SKY_SPEC_SATURATION:		++stat.saturatedSkySpecNum; return;
                case SKY_SPEC_DARK:				++stat.darkSkySpecNum; return;
                case SKY_SPEC_TOO_LONG_EXPTIME:	++stat.tooLongExpTime; return;
                case COMPLETENESS_LOW:			++stat.lowCompletenessNum; return;
                case NO_PLUME:					++stat.noPlumeNum; return;
                };
            }
        }

        // If the instrument is not in the list then insert it!!
        CInstrumentStats stat;
        stat.serial.Format(serial);
        switch (reason)
        {
        case SKY_SPEC_SATURATION:		++stat.saturatedSkySpecNum; break;
        case SKY_SPEC_DARK:				++stat.darkSkySpecNum; break;
        case SKY_SPEC_TOO_LONG_EXPTIME:	++stat.tooLongExpTime; break;
        case COMPLETENESS_LOW:			++stat.lowCompletenessNum; break;
        case NO_PLUME:					++stat.noPlumeNum; break;
        };
        m_instrumentStats.push_back(stat);
    }

    singleLock.Unlock();
}

/** Inserts information on a accepted scan from a certain instrument into the database. */
void CPostProcessingStatistics::InsertAcception(const novac::CString& serial)
{

    novac::CSingleLock singleLock(&g_processingStatCritSect);
    singleLock.Lock();
    if (singleLock.IsLocked())
    {

        // look for the correct instrument
        std::list<CInstrumentStats>::const_iterator pos = m_instrumentStats.begin();
        while (pos != m_instrumentStats.end())
        {
            CInstrumentStats& stat = (CInstrumentStats&)*(pos++);

            // this is the instrument. Insert the new data
            if (Equals(stat.serial, serial))
            {
                ++stat.acceptedScans;
                singleLock.Unlock();
                return;
            }
        }

        // If the instrument is not in the list then insert it!!
        CInstrumentStats stat;
        stat.serial.Format(serial);
        ++stat.acceptedScans;
        m_instrumentStats.push_back(stat);
    }

    singleLock.Unlock();
}

/** Retrieves the number of rejected full scans due to the specified reason */
unsigned long CPostProcessingStatistics::GetRejectionNum(const novac::CString& serial, const REASON_FOR_REJECTION& reason)
{

    // look for the correct instrument
    std::list<CInstrumentStats>::const_iterator pos = m_instrumentStats.begin();
    while (pos != m_instrumentStats.end())
    {
        CInstrumentStats& stat = (CInstrumentStats&)*(pos++);

        // this is the instrument. Retrieve the data
        if (Equals(stat.serial, serial))
        {
            switch (reason)
            {
            case SKY_SPEC_SATURATION:		return stat.saturatedSkySpecNum;
            case SKY_SPEC_DARK:				return stat.darkSkySpecNum;
            case SKY_SPEC_TOO_LONG_EXPTIME:	return stat.tooLongExpTime;
            case COMPLETENESS_LOW:			return stat.lowCompletenessNum;
            case NO_PLUME:					return stat.noPlumeNum;
            };
        }
    }

    // shouldn't happen
    return 0;
}

/** Retrieves the number of accepted full scans */
unsigned long CPostProcessingStatistics::GetAcceptionNum(const novac::CString& serial)
{

    // look for the correct instrument
    std::list<CInstrumentStats>::const_iterator pos = m_instrumentStats.begin();
    while (pos != m_instrumentStats.end())
    {
        CInstrumentStats& stat = (CInstrumentStats&)*(pos++);

        // this is the instrument. Retrieve the data
        if (Equals(stat.serial, serial))
        {
            return stat.acceptedScans;
        }
    }

    // shouldn't happen
    return 0;

}

/** Inserts the successful evaluation of a single spectrum into the statistics.
    This also increases the counter on the total amount of time used on
    evaluating spectra */
void CPostProcessingStatistics::InsertEvaluatedSpectrum(double timeUsed)
{
    ++nSpectraEvaluated;
    timeSpentOnEvaluations += timeUsed;
}

/** Creates a small output file containing the statistical results */
void CPostProcessingStatistics::WriteStatToFile(const novac::CString& file)
{
    novac::CSingleLock singleLock(&g_processingStatCritSect);
    singleLock.Lock();
    if (singleLock.IsLocked())
    {

        // open the file
        FILE* f = NULL;
        if (Filesystem::IsExistingFile(file))
        {
            f = fopen(file, "a");
        }
        else
        {
            f = fopen(file, "w");
        }
        if (f == NULL)
        {
            singleLock.Unlock();
            return;
        }

        // for each instrument processed, write the info we have on it...
        std::list <CInstrumentStats>::const_iterator pos = m_instrumentStats.begin();
        while (pos != m_instrumentStats.end())
        {
            CInstrumentStats& instr = (CInstrumentStats&)*(pos++);

            fprintf(f, "Instrument: %s\n", (const char*)instr.serial);
            fprintf(f, "\t#Accepted scans: %lu\n", instr.acceptedScans);
            fprintf(f, "\t#Rejected scans:\n");
            fprintf(f, "\t\t%lu due to too long exposure time\n", instr.tooLongExpTime);
            fprintf(f, "\t\t%lu due to saturated sky spectrum\n", instr.saturatedSkySpecNum);
            fprintf(f, "\t\t%lu due to too dark sky spectrum\n", instr.darkSkySpecNum);
            fprintf(f, "\t\t%lu due to no plume seen\n", instr.noPlumeNum);
            fprintf(f, "\t\t%lu due to too low completeness\n", instr.lowCompletenessNum);
        }

        // The timings...
        fprintf(f, "Total Number of Spectra evaluated: %lu\n", nSpectraEvaluated);
        fprintf(f, "Total Time spent on evaluating spectra: %.2lf [s] ( %.2lf mseconds / spectrum)\n", timeSpentOnEvaluations / 1000.0, timeSpentOnEvaluations / (double)nSpectraEvaluated);

        // remember to close the file
        fclose(f);
    }

    singleLock.Unlock();
}
