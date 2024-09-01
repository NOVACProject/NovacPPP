#include <PPPLib/Configuration/InstrumentLocation.h>

using namespace novac;

namespace Configuration
{
CInstrumentLocation::CInstrumentLocation()
{
    Clear();
}

void CInstrumentLocation::Clear()
{
    this->m_altitude = 0;
    this->m_compass = 0.0;
    this->m_coneangle = 0.0;
    this->m_latitude = 0.0;
    this->m_locationName = "";
    this->m_longitude = 0.0;
    this->m_tilt = 0.0;
    this->m_validFrom = CDateTime(0000, 00, 00, 0, 0, 0);
    this->m_validTo = CDateTime(9999, 12, 31, 23, 59, 59);
    this->m_volcano = "";
    this->m_instrumentType = INSTRUMENT_TYPE::INSTR_GOTHENBURG;
    this->m_spectrometerModel = "S2000";
}

CInstrumentLocation& CInstrumentLocation::operator =(const CInstrumentLocation& l2)
{
    this->m_altitude = l2.m_altitude;
    this->m_compass = l2.m_compass;
    this->m_coneangle = l2.m_coneangle;
    this->m_latitude = l2.m_latitude;
    this->m_locationName.Format(l2.m_locationName);
    this->m_longitude = l2.m_longitude;
    this->m_tilt = l2.m_tilt;
    this->m_validFrom = l2.m_validFrom;
    this->m_validTo = l2.m_validTo;
    this->m_volcano.Format(l2.m_volcano);
    this->m_instrumentType = l2.m_instrumentType;
    this->m_spectrometerModel = l2.m_spectrometerModel;

    return *this;
}

}
