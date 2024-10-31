#pragma once

#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/GPSData.h>
#include <SpectralEvaluation/NovacEnums.h>
#include <PPPLib/MFC/CString.h>

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
    novac::CString m_locationName = "";

    /** The volcano that was monitored from this location */
    novac::CString m_volcano = "";

    /** The latitude of the location */
    double m_latitude = 0.0;

    /** The longitude of the location */
    double m_longitude = 0.0;

    /** The altitude of the location, in meters above sea level */
    int m_altitude = 0;

    /** The compass direction of the instrument.
        this is in degrees from north, counting clockwise */
    double m_compass = 0.0;

    /** The cone-angle of the instrument during the setup here, in degrees.
        90 corresponds to a flat scanner and 60 corresponds to a 'normal'
        conical scanner */
    double m_coneangle = 0.0;

    /** The tilt of the instrument. In degrees from horizontal. This is typically 0. */
    double m_tilt = 0.0;

    /** The type of instrument that was used */
    novac::NovacInstrumentType m_instrumentType = novac::NovacInstrumentType::Gothenburg;

    /** The type of spectrometer */
    std::string m_spectrometerModel = "S2000";

    /** Time stamps, during which this location information
        is resonable/useful */
    novac::CDateTime m_validFrom = novac::CDateTime(0000, 00, 00, 0, 0, 0);
    novac::CDateTime m_validTo = novac::CDateTime(9999, 12, 31, 23, 59, 59);

    // ---------------------------------------------
    // --------------- PUBLIC METHODS --------------
    // ---------------------------------------------

    /** Clears all the information in this data structure */
    void Clear();

    novac::CGPSData GpsData() const
    {
        return novac::CGPSData(this->m_latitude, this->m_longitude, this->m_altitude);
    }

    /** Assignment operator */
    CInstrumentLocation& operator=(const CInstrumentLocation& l2);
};
}