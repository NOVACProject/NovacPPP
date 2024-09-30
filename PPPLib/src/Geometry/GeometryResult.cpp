#include <PPPLib/Geometry/GeometryResult.h>

namespace Geometry
{

CGeometryResult& CGeometryResult::operator=(const CGeometryResult& gr)
{
    this->m_averageStartTime = gr.m_averageStartTime;
    this->m_startTimeDifference = gr.m_startTimeDifference;
    m_plumeAltitude = gr.m_plumeAltitude;
    m_plumeAltitudeError = gr.m_plumeAltitudeError;
    m_windDirection = gr.m_windDirection;
    m_windDirectionError = gr.m_windDirectionError;

    m_instr1.Format(gr.m_instr1);
    m_instr2.Format(gr.m_instr2);

    m_plumeCentre1 = gr.m_plumeCentre1;
    m_plumeCentreError1 = gr.m_plumeCentreError1;
    m_plumeCentre2 = gr.m_plumeCentre2;
    m_plumeCentreError2 = gr.m_plumeCentreError2;

    m_calculationType = gr.m_calculationType;

    return *this;
}
}
