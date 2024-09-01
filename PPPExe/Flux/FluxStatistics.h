#pragma once

#ifndef FLUXSTATISTICS_H
#define FLUXSTATISTICS_H

#include <PPPLib/MFC/CString.h>
#include <PPPLib/MFC/CList.h>
#include "FluxResult.h"

namespace Flux
{
class CFluxStatistics
{
public:
    CFluxStatistics(void);
    ~CFluxStatistics(void);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------


    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Clears all information here */
    void Clear();

    /** Attaches the supplied list of flux results to the current set
        of measured data. */
    void AttachFluxList(novac::CList<CFluxResult, CFluxResult&>& calculatedFluxes);

    /** Attaches the given flux result to the current set of
        measured data */
    void AttachFlux(const CFluxResult& result);

    /** Calculates statistics on the statistics we have here and writes
        the results to file. */
    void WriteFluxStat(const novac::CString& fileName);

private:
    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    class CMeasurementDay
    {
    public:
        CMeasurementDay();
        ~CMeasurementDay();
        novac::CDateTime day; // the date of the measurement
        novac::CList <CFluxResult> fluxList;
        static void GetHeaderLine(novac::CString& str, novac::CList <novac::CString, novac::CString&>& instruments);
        void GetStatistics(novac::CString& str, novac::CList <novac::CString, novac::CString&>& instruments);
        CMeasurementDay& operator=(const CMeasurementDay& m);
    };

    /** The list of measurement results.
        This is sorted by the 'day' in increasing order */
    novac::CList <CMeasurementDay, CMeasurementDay&> m_measurements;

    /** The list of instruments used. This is updated together with
        'm_measurements' when calling 'AttachFlux' */
    novac::CList <novac::CString, novac::CString&> m_instruments;

    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------



};
}

#endif