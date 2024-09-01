#pragma once

#include <PPPLib/Configuration/InstrumentType.h>
#include <SpectralEvaluation/DateTime.h>
#include <PPPLib/MFC/CString.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

/**
    The class <b>CInstrumentLocation</b> is used to store information
    about how an instrument was setup. This stores the lat and long
    of the location where the instrument was setup, the name of the
    place, the compass-direction that the instrument was poining
    in and the cone-angle of the instrument.
*/

namespace Configuration
{
class CInstrumentLocation
{
public:
    CInstrumentLocation();

    // -----------------------------------------
    // -------------- PUBLIC DATA --------------
    // -----------------------------------------

    /** The name of the location */
    novac::CString m_locationName;

    /** The volcano that was monitored from this location */
    novac::CString m_volcano;

    /** The latitude of the location */
    double  m_latitude;

    /** The longitude of the location */
    double  m_longitude;

    /** The altitude of the location, in meters above sea level */
    int     m_altitude;

    /** The compass direction of the instrument.
        this is in degrees from north, counting clockwise */
    double  m_compass;

    /** The cone-angle of the instrument during the setup here, in degrees.
        90 corresponds to a flat scanner and 60 corresponds to a 'normal'
        conical scanner */
    double  m_coneangle;

    /** The tilt of the instrument. In degrees from horizontal. This is typically 0. */
    double  m_tilt;

    /** The type of instrument that was used */
    INSTRUMENT_TYPE m_instrumentType;

    /** The type of spectrometer */
    std::string m_spectrometerModel;

    /** Time stamps, during which this location information
        is resonable/useful */
    novac::CDateTime m_validFrom, m_validTo;

    // ---------------------------------------------
    // --------------- PUBLIC METHODS --------------
    // ---------------------------------------------

    /** Clears all the information in this data structure */
    void Clear();

    /** Assignment operator */
    CInstrumentLocation& operator=(const CInstrumentLocation& l2);
};
}