#include <PPPLib/Geometry/PlumeHeight.h>

namespace Geometry
{

CPlumeHeight& CPlumeHeight::operator=(const CPlumeHeight& ph2)
{
    this->m_plumeAltitude = ph2.m_plumeAltitude;
    this->m_plumeAltitudeError = ph2.m_plumeAltitudeError;
    this->m_plumeAltitudeSource = ph2.m_plumeAltitudeSource;

    this->m_validFrom = ph2.m_validFrom;
    this->m_validTo = ph2.m_validTo;

    return *this;
}
}
