#pragma once

#include <PPPLib/Meteorology/MeteorologySource.h>
#include <SpectralEvaluation/DateTime.h>

namespace Geometry
{

/** The struct PlumeHeight is intended to hold information about the
    height of the plume at a given time.
    It should be used as a container for this kind of data. */
struct PlumeHeight
{
public:

    PlumeHeight() = default;
    PlumeHeight(const PlumeHeight& other) = default;
    PlumeHeight(PlumeHeight&& other) = default;

    PlumeHeight& operator=(const PlumeHeight& ph) = default;
    PlumeHeight& operator=(PlumeHeight&& ph) = default;

    /** The altitude of the plume. In meters above sea level. */
    double m_plumeAltitude = 1000.0;

    /** The uncertainty of the altitude of the plume.
            In meters (above sea level). */
    double m_plumeAltitudeError = 1000.0;

    /** The source of our knowledge of this altitude. */
    Meteorology::MeteorologySource m_plumeAltitudeSource = Meteorology::MeteorologySource::Default;

    /** The time range over which this information of the plume is valid */
    novac::CDateTime m_validFrom = novac::CDateTime(0, 0, 0, 0, 0, 0);
    novac::CDateTime m_validTo = novac::CDateTime(9999, 12, 31, 23, 59, 59);
};
}
