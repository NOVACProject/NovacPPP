#include <PPPLib/Meteorology/WindField.h>

namespace Meteorology
{

// global function that converts a MeteorologySource item to string
void MetSourceToString(const MeteorologySource src, novac::CString& str)
{
    if (MeteorologySource::User == src)
        str.Format("user");
    else if (MeteorologySource::Default == src)
        str.Format("default");
    else if (MeteorologySource::EcmwfForecast == src)
        str.Format("ecmwf_forecast");
    else if (MeteorologySource::EcmwfAnalysis == src)
        str.Format("ecmwf_analysis");
    else if (MeteorologySource::DualBeamMeasurement == src)
        str.Format("dual_beam_measurement");
    else if (MeteorologySource::ModelledWrf == src)
        str.Format("model_wrf");
    else if (MeteorologySource::GeometryCalculationTwoInstruments == src)
        str.Format("geometry_calc");
    else if (MeteorologySource::GeometryCalculationSingleInstrument == src)
        str.Format("geometry_calc_single_instr");
    else if (MeteorologySource::NoaaGdas == src)
        str.Format("noaa_gdas");
    else if (MeteorologySource::NoaaFnl == src)
        str.Format("noaa_fnl");
    else if (MeteorologySource::None == src)
        str.Format("none");
    else
        str.Format("unknown");
}

MeteorologySource StringToMetSource(const novac::CString& str)
{
    novac::CString trimmedStr(str);
    trimmedStr.Trim(); // remove blanks in the beginning and in the end

    if (Equals(trimmedStr, "user"))
    {
        return MeteorologySource::User;

    }
    else if (Equals(trimmedStr, "none"))
    {
        return MeteorologySource::None;

    }
    else if (Equals(trimmedStr, "default"))
    {
        return MeteorologySource::Default;

    }
    else if (Equals(trimmedStr, "default"))
    {
        return MeteorologySource::User;

    }
    else if (Equals(trimmedStr, "ecmwf_forecast"))
    {
        return MeteorologySource::EcmwfForecast;

    }
    else if (Equals(trimmedStr, "ecmwf_analysis"))
    {
        return MeteorologySource::EcmwfAnalysis;

    }
    else if (Equals(trimmedStr, "ecmwf_postanalysis"))
    {
        return MeteorologySource::EcmwfAnalysis;

    }
    else if (Equals(trimmedStr, "dual_beam_measurement"))
    {
        return MeteorologySource::DualBeamMeasurement;

    }
    else if (Equals(trimmedStr, "model_wrf"))
    {
        return MeteorologySource::ModelledWrf;

    }
    else if (Equals(trimmedStr, "noaa_gdas"))
    {
        return MeteorologySource::NoaaGdas;

    }
    else if (Equals(trimmedStr, "noaa_fnl"))
    {
        return MeteorologySource::NoaaFnl;

    }
    else if (Equals(trimmedStr, "geometry_calc"))
    {
        return MeteorologySource::GeometryCalculationTwoInstruments;

    }
    else if (Equals(trimmedStr, "geometry_calc_single_instr"))
    {
        return MeteorologySource::GeometryCalculationSingleInstrument;

    }
    else
    {
        return MeteorologySource::None;
    }
}

int GetSourceQuality(MeteorologySource src)
{
    MeteorologySource list[] = {
        MeteorologySource::DualBeamMeasurement,
        MeteorologySource::GeometryCalculationTwoInstruments,
        MeteorologySource::GeometryCalculationSingleInstrument,
        MeteorologySource::ModelledWrf,
        MeteorologySource::EcmwfAnalysis,
        MeteorologySource::NoaaGdas,
        MeteorologySource::EcmwfForecast,
        MeteorologySource::NoaaFnl,
        MeteorologySource::User,
        MeteorologySource::Default,
        MeteorologySource::None
    };

    const int numberOfSources = static_cast<int>(MeteorologySource::NumberOfSources);

    for (int k = 0; k < numberOfSources; ++k)
    {
        if (src == list[k])
        {
            return (numberOfSources - k - 1);
        }
    }
    return -1; // shouldn't happen
}

WindField::WindField(double windSpeed, MeteorologySource windSpeedSrc, double windDir, MeteorologySource windDirSrc, const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double lat, double lon, double alt)
    : m_windSpeed(windSpeed), m_windSpeedError(0.0), m_windSpeedSource(windSpeedSrc),
    m_windDirection(windDir), m_windDirectionError(0.0), m_windDirectionSource(windDirSrc),
    m_validFrom(validFrom), m_validTo(validTo), m_location(lat, lon, alt)
{}

WindField::WindField(double windSpeed, double windSpeedErr, MeteorologySource windSpeedSrc, double windDir, double windDirErr, MeteorologySource windDirSrc, const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double lat, double lon, double alt)
    : m_windSpeed(windSpeed), m_windSpeedError(windSpeedErr), m_windSpeedSource(windSpeedSrc), 
    m_windDirection(windDir), m_windDirectionError(windDirErr), m_windDirectionSource(windDirSrc),
    m_validFrom(validFrom), m_validTo(validTo), m_location(lat, lon, alt)
{
}

void WindField::SetWindSpeed(double ws, MeteorologySource source)
{
    this->m_windSpeed = ws;
    this->m_windSpeedSource = source;
}

void WindField::SetWindDirection(double wd, MeteorologySource source)
{
    this->m_windDirection = wd;
    this->m_windDirectionSource = source;

}

void WindField::SetValidTimeFrame(const novac::CDateTime& from, const novac::CDateTime& to)
{
    this->m_validFrom = from;
    this->m_validTo = to;
}

/** Gets the wind-speed */
double WindField::GetWindSpeed() const
{
    return this->m_windSpeed;
}

/** Gets the source of the wind-speed */
MeteorologySource WindField::GetWindSpeedSource() const
{
    return this->m_windSpeedSource;
}

/** Gets the source of the wind-speed */
void WindField::GetWindSpeedSource(novac::CString& str) const
{
    return MetSourceToString(m_windSpeedSource, str);
}

/** Gets the estimate for the total error in the wind-speed */
double WindField::GetWindSpeedError() const
{
    return this->m_windSpeedError;
}

/** Sets the estimate for the total error in the wind-speed */
void WindField::SetWindSpeedError(double err)
{
    this->m_windSpeedError = err;
}


/** Gets the wind-direction */
double WindField::GetWindDirection() const
{
    return this->m_windDirection;
}

/** Gets the source of the wind-direction */
MeteorologySource WindField::GetWindDirectionSource() const
{
    return this->m_windDirectionSource;
}

/** Gets the source of the wind-direction */
void WindField::GetWindDirectionSource(novac::CString& str) const
{
    return MetSourceToString(m_windDirectionSource, str);
}

/** Gets the estimate for the total error in the wind-direction */
double WindField::GetWindDirectionError() const
{
    return this->m_windDirectionError;
}

/** Sets the estimate for the total error in the wind-direction */
void WindField::SetWindDirectionError(double err)
{
    this->m_windDirectionError = err;
}

/** Gets the time and date for which this wind-field is valid */
void WindField::GetValidTimeFrame(novac::CDateTime& from, novac::CDateTime& to) const
{
    from = this->m_validFrom;
    to = this->m_validTo;
}

/** Gets the position that this wind-field is valid for */
void WindField::GetValidPosition(double& lat, double& lon, double& alt) const
{
    lat = this->m_location.m_latitude;
    lon = this->m_location.m_longitude;
    alt = this->m_location.m_altitude;
}

} // namespace Meteorology