#pragma once

#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/NovacEnums.h>
#include <PPPLib/Meteorology/MeteorologySource.h>
#include <PPPLib/Meteorology/WindField.h>
#include <PPPLib/Geometry/PlumeHeight.h>

enum class FluxQuality
{
    Green,
    Yellow,
    Red
};

namespace Flux
{

/** The data structure FluxResult stores the results from a flux-calculation of a scan.
    The struct holds the values of all the parameters used in the calculation (except for the measurement itself)
    and the result of the measurement. */
struct FluxResult
{
public:
    /** The calculated flux, in kg/s */
    double m_flux = 0.0;

    /** The quality of this flux measurement */
    FluxQuality m_fluxQualityFlag = FluxQuality::Green;

    /** The error in flux due to the uncertainty in
        wind speed and wind direction */
    double m_fluxError_Wind = 0.0;

    /** The error in flux due to the uncertainty in
        plume height */
    double m_fluxError_PlumeHeight = 0.0;

    /** The wind field that was used to calculate this flux */
    Meteorology::WindField m_windField;

    /** The plume height that was used to calculate this flux */
    Geometry::PlumeHeight m_plumeHeight;

    /** The number of good spectra in this measurement. This is also
        the number of column-values that were used to calculate the
        flux */
    unsigned int m_numGoodSpectra = 0U;

    /** The cone-angle of the scanner that collected this scan */
    double m_coneAngle = 0.0;

    /** The tilt of the scanner that collected this scan */
    double m_tilt = 0.0;

    /** The compass-direction of the scanner that collected this scan */
    double m_compass = 0.0;

    /** The date and time (UTC) when the measurement was started */
    novac::CDateTime m_startTime;

    /** The date and time (UTC) when the measurement was finished */
    novac::CDateTime m_stopTime;

    /** The instrument that collected this scan */
    std::string m_instrument;

    /** The type of the instrument that collected this scan */
    novac::NovacInstrumentType m_instrumentType = novac::NovacInstrumentType::Gothenburg;

    /** The calculated offset of the scan */
    double m_scanOffset = 0.0;

    /** The calculated plume completeness */
    double m_completeness = 0.0;

    /** The calculated plume centre position */
    double m_plumeCentre[2] = { 0.0, 0.0 };

};
}