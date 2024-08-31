#pragma once

#include <PPPLib/Configuration/InstrumentConfiguration.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/SpectrometerId.h>

/**
    The class <b>CNovacPPPConfiguration</b> is the main configuration
        for the NovacPPP. This configuration stores an array of
        configured instruments (for each of these instruments there
        is a set of fit-windows and locations configured).

        This also defines the directories that we should use to store
        the output data and the temporary files.
*/
namespace Configuration
{
class CNovacPPPConfiguration
{
public:
    CNovacPPPConfiguration();

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The array of instruments that are configured. */
    std::vector<CInstrumentConfiguration> m_instrument;

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Returns the number of configured instruments */
    int NumberOfInstruments() const { return static_cast<int>(m_instrument.size()); }

    /** Retrieves the CInstrumentConfiguration that is connected with a given serial-number.
        @return a pointer to the found CInstrumentConfiguration.
            If none is found then return value is nullptr. */
    const CInstrumentConfiguration* GetInstrument(const novac::CString& serial) const;
    const CInstrumentConfiguration* GetInstrument(const std::string& serial) const;

    /** Retrieves the CInstrumentLocation that is valid for the given instrument and
        for the given time
        @return 0 if successful otherwise non-zero
    */
    int GetInstrumentLocation(const novac::CString& serial, const novac::CDateTime& dateAndTime, CInstrumentLocation& instrLocation) const;

    /** Retrieves the CFitWindow that is valid for the given instrument and
        for the given time
        if 'fitWindowName' is not NULL then only the fit-window with the specified
            name will be returned.
        if 'fitWindowName' is NULL then the first fit-window valid at the given time
            will be returned.
        @return 0 if successful otherwise non-zero */
    int GetFitWindow(const novac::CString& serial, int channel, const novac::CDateTime& dateAndTime, novac::CFitWindow& window, const novac::CString* fitWindowName = NULL) const;

    /** Retrieves the CDarkSettings that is valid for the given instrument and
        for the given time

        @return 0 if successful otherwise non-zero
    */
    int GetDarkCorrection(const novac::CString& serial, const novac::CDateTime& dateAndTime, CDarkSettings& settings) const;

};
}