#pragma once

#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Definitions.h>
#include <SpectralEvaluation/NullableValue.h>
#include <PPPLib/Meteorology/MeteorologySource.h>
#include <PPPLib/MFC/CString.h>

namespace Geometry
{
class CGeometryResult
{
public:
    /** Assignment operator */
    CGeometryResult& operator=(const CGeometryResult& gr);

    /** How this calculation was made. Can be either one of
        GeometryCalculationTwoInstruments (default) or
        GeometryCalculationSingleInstrument */
    Meteorology::MeteorologySource m_calculationType = Meteorology::MeteorologySource::GeometryCalculationTwoInstruments;

    /** The average of the starting-time of the
            two scans that were combined to generate the result
            (seconds since midnight, UTC) */
    novac::CDateTime m_averageStartTime;

    /** The difference in start-time between the two
        scans that were combined to make this measurement.
        In seconds. */
    int	m_startTimeDifference = 0;

    /** The two instruments that were used to derive this
        geometry result. */
    std::string m_instrumentSerial1 = "";
    std::string m_instrumentSerial2 = "";

    /** The plume centre-angles for the two scans that were combined */
    novac::Nullable<double> m_plumeCentre1;
    novac::Nullable<double> m_plumeCentre2;

    /** The estimated error in the plume-centre angles
        for each of the two scans that were combined */
    novac::Nullable<double> m_plumeCentreError1;
    novac::Nullable<double> m_plumeCentreError2;

    /** The calculated plume height (in meters above sea level) */
    double m_plumeAltitude = 0.0;

    /** The estimated error in plume height (in meters) */
    double m_plumeAltitudeError = 0.0;

    /** The calculated wind-direction (degrees from north) */
    double m_windDirection = 0.0;

    /** The estimated error in the calculated wind-direction (degrees) */
    double m_windDirectionError = 0.0;
};
}