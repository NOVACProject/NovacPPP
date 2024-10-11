#include <PPPLib/Flux/FluxResult.h>

namespace Flux
{

CFluxResult& CFluxResult::operator=(const CFluxResult& res)
{
    m_flux = res.m_flux;
    m_fluxQualityFlag = res.m_fluxQualityFlag;

    m_fluxError_Wind = res.m_fluxError_Wind;
    m_fluxError_PlumeHeight = res.m_fluxError_PlumeHeight;

    m_windField = res.m_windField;
    m_plumeHeight = res.m_plumeHeight;

    m_numGoodSpectra = res.m_numGoodSpectra;

    m_coneAngle = res.m_coneAngle;
    m_tilt = res.m_tilt;
    m_compass = res.m_compass;
    m_volcano = res.m_volcano;

    m_startTime = res.m_startTime;
    m_stopTime = res.m_stopTime;

    m_instrument = res.m_instrument;

    m_instrumentType = res.m_instrumentType;
    m_scanOffset = res.m_scanOffset;
    m_completeness = res.m_completeness;
    m_plumeCentre[0] = res.m_plumeCentre[0];
    m_plumeCentre[1] = res.m_plumeCentre[1];


    return *this;
}

}
