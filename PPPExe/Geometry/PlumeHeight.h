#pragma once

#include <PPPLib/Meteorology/MeteorologySource.h>
#include <SpectralEvaluation/DateTime.h>

#ifndef PLUMEHEIGHT_H
#define PLUMEHEIGHT_H


namespace Geometry {

    /** The class <b>CPlumeHeight</b> is intended to hold information about the
        height of the plume at a given time.
        It should be used as a container for this kind of data. */
    class CPlumeHeight
    {
    public:
        /** Default constructor */
        CPlumeHeight();
        ~CPlumeHeight() = default;

        /** assignment operator */
        CPlumeHeight &operator=(const CPlumeHeight &ph);

        /** The altitude of the plume. In meters above sea level. */
        double  m_plumeAltitude;

        /** The uncertainty of the altitude of the plume.
                In meters (above sea level). */
        double  m_plumeAltitudeError;

        /** The source of our knowledge of this altitude. */
        Meteorology::MET_SOURCE m_plumeAltitudeSource;

        /** The time range over which this information of the plume is valid */
        novac::CDateTime   m_validFrom;
        novac::CDateTime   m_validTo;
    };
}

#endif