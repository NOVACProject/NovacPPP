#pragma once

#include <PPPLib/Configuration/InstrumentConfiguration.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/SpectrometerId.h>
#include <PPPLib/Exceptions.h>

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
    const CInstrumentConfiguration* GetInstrument(const std::string& serial) const;

    /** Retrieves and returns the CInstrumentLocation that is valid for the given instrument and for the given time
    *   @throws PPPLib::NotFoundException if the instrument could not be found */
    CInstrumentLocation GetInstrumentLocation(const std::string& serial, const novac::CDateTime& dateAndTime) const;

    /** Retrieves the CFitWindow that is valid for the given instrument and for the given time
    *   if 'fitWindowName' is not NULL then only the fit-window with the specified name will be returned.
    *   if 'fitWindowName' is NULL then the first fit-window valid at the given time will be returned.
    *   @throws PPPLib::NotFoundException if the instrument could not be found */
    novac::CFitWindow GetFitWindow(const std::string& serial, int channel, const novac::CDateTime& dateAndTime, const novac::CString* fitWindowName = NULL) const;

    /** Retrieves the CDarkSettings that is valid for the given instrument and for the given time
    *   @throws PPPLib::NotFoundException if the instrument could not be found */
    CDarkSettings GetDarkCorrection(const std::string& serial, const novac::CDateTime& dateAndTime) const;

};
}