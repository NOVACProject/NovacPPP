#include <PPPLib/PostProcessingStatistics.h>
#include <PPPLib/MFC/CCriticalSection.h>
#include <PPPLib/MFC/CSingleLock.h>
#include <PPPLib/File/Filesystem.h>

novac::CCriticalSection g_processingStatCritSect; // synchronization access to the processing statistics

void CPostProcessingStatistics::InsertRejection(const novac::CString& serial, ReasonForScanRejection reason)
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
                case ReasonForScanRejection::SkySpectrumSaturated:       ++stat.saturatedSkySpecNum; return;
                case ReasonForScanRejection::SkySpectrumDark:             ++stat.darkSkySpecNum; return;
                case ReasonForScanRejection::SkySpectrumTooLongExposureTime: ++stat.tooLongExpTime; return;
                case ReasonForScanRejection::CompletenessLow:          ++stat.lowCompletenessNum; return;
                case ReasonForScanRejection::NoPlume:                  ++stat.noPlumeNum; return;
                };
            }
        }

        // If the instrument is not in the list then insert it!!
        CInstrumentStats stat;
        stat.serial.Format(serial);
        switch (reason)
        {
        case ReasonForScanRejection::SkySpectrumSaturated:       ++stat.saturatedSkySpecNum; break;
        case ReasonForScanRejection::SkySpectrumDark:             ++stat.darkSkySpecNum; break;
        case ReasonForScanRejection::SkySpectrumTooLongExposureTime: ++stat.tooLongExpTime; break;
        case ReasonForScanRejection::CompletenessLow:          ++stat.lowCompletenessNum; break;
        case ReasonForScanRejection::NoPlume:                  ++stat.noPlumeNum; break;
        };
        m_instrumentStats.push_back(stat);
    }

    singleLock.Unlock();
}

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

unsigned long CPostProcessingStatistics::GetRejectionNum(const novac::CString& serial, ReasonForScanRejection reason)
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
            case ReasonForScanRejection::SkySpectrumSaturated:       return stat.saturatedSkySpecNum;
            case ReasonForScanRejection::SkySpectrumDark:             return stat.darkSkySpecNum;
            case ReasonForScanRejection::SkySpectrumTooLongExposureTime: return stat.tooLongExpTime;
            case ReasonForScanRejection::CompletenessLow:          return stat.lowCompletenessNum;
            case ReasonForScanRejection::NoPlume:                  return stat.noPlumeNum;
            };
        }
    }

    // shouldn't happen
    return 0;
}

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

void CPostProcessingStatistics::InsertEvaluatedSpectrum(double timeUsed)
{
    ++nSpectraEvaluated;
    timeSpentOnEvaluations += timeUsed;
}

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
