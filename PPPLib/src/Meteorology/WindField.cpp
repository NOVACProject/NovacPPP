#include <PPPLib/Meteorology/WindField.h>

using namespace Meteorology;

// global function that converts a MET_SOURCE item to string
void Meteorology::MetSourceToString(const MET_SOURCE src, novac::CString& str)
{
    if (MET_USER == src)
        str.Format("user");
    else if (MET_DEFAULT == src)
        str.Format("default");
    else if (MET_ECMWF_FORECAST == src)
        str.Format("ecmwf_forecast");
    else if (MET_ECMWF_ANALYSIS == src)
        str.Format("ecmwf_analysis");
    else if (MET_DUAL_BEAM_MEASUREMENT == src)
        str.Format("dual_beam_measurement");
    else if (MET_MODEL_WRF == src)
        str.Format("model_wrf");
    else if (MET_GEOMETRY_CALCULATION == src)
        str.Format("geometry_calc");
    else if (MET_GEOMETRY_CALCULATION_SINGLE_INSTR == src)
        str.Format("geometry_calc_single_instr");
    else if (MET_NOAA_GDAS == src)
        str.Format("noaa_gdas");
    else if (MET_NOAA_FNL == src)
        str.Format("noaa_fnl");
    else if (MET_NONE == src)
        str.Format("none");
    else
        str.Format("unknown");
}

MET_SOURCE Meteorology::StringToMetSource(const novac::CString& str)
{
    novac::CString trimmedStr(str);
    trimmedStr.Trim(); // remove blanks in the beginning and in the end

    if (Equals(trimmedStr, "user"))
    {
        return MET_USER;

    }
    else if (Equals(trimmedStr, "none"))
    {
        return MET_NONE;

    }
    else if (Equals(trimmedStr, "default"))
    {
        return MET_DEFAULT;

    }
    else if (Equals(trimmedStr, "default"))
    {
        return MET_USER;

    }
    else if (Equals(trimmedStr, "ecmwf_forecast"))
    {
        return MET_ECMWF_FORECAST;

    }
    else if (Equals(trimmedStr, "ecmwf_analysis"))
    {
        return MET_ECMWF_ANALYSIS;

    }
    else if (Equals(trimmedStr, "ecmwf_postanalysis"))
    {
        return MET_ECMWF_ANALYSIS;

    }
    else if (Equals(trimmedStr, "dual_beam_measurement"))
    {
        return MET_DUAL_BEAM_MEASUREMENT;

    }
    else if (Equals(trimmedStr, "model_wrf"))
    {
        return MET_MODEL_WRF;

    }
    else if (Equals(trimmedStr, "noaa_gdas"))
    {
        return MET_NOAA_GDAS;

    }
    else if (Equals(trimmedStr, "noaa_fnl"))
    {
        return MET_NOAA_FNL;

    }
    else if (Equals(trimmedStr, "geometry_calc"))
    {
        return MET_GEOMETRY_CALCULATION;

    }
    else if (Equals(trimmedStr, "geometry_calc_single_instr"))
    {
        return MET_GEOMETRY_CALCULATION_SINGLE_INSTR;

    }
    else
    {
        return MET_NONE;
    }
}

/** Retrieves the judged quality of a given source */
int Meteorology::GetSourceQuality(MET_SOURCE src)
{
    int list[] = {
        MET_DUAL_BEAM_MEASUREMENT,
        MET_GEOMETRY_CALCULATION,
        MET_GEOMETRY_CALCULATION_SINGLE_INSTR,
        MET_MODEL_WRF,
        MET_ECMWF_ANALYSIS,
        MET_NOAA_GDAS,
        MET_ECMWF_FORECAST,
        MET_NOAA_FNL,
        MET_USER,
        MET_DEFAULT,
        MET_NONE
    };
    for (int k = 0; k < MET_NUMBER_OF_DEFINED_SOURCES; ++k)
    {
        if (src == list[k])
        {
            return (MET_NUMBER_OF_DEFINED_SOURCES - k - 1);
        }
    }
    return -1; // shouldn't happen
}

WindField::WindField(double windSpeed, MET_SOURCE windSpeedSrc, double windDir, MET_SOURCE windDirSrc, const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double lat, double lon, double alt)
    : m_windSpeed(windSpeed), m_windSpeedSource(windSpeedSrc), m_windSpeedError(0.0),
    m_windDirection(windDir), m_windDirectionSource(windDirSrc), m_windDirectionError(0.0),
    m_validFrom(validFrom), m_validTo(validTo), m_location(lat, lon, alt)
{}

WindField::WindField(double windSpeed, double windSpeedErr, MET_SOURCE windSpeedSrc, double windDir, double windDirErr, MET_SOURCE windDirSrc, const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double lat, double lon, double alt)
    : m_windSpeed(windSpeed), m_windSpeedSource(windSpeedSrc), m_windSpeedError(windSpeedErr),
    m_windDirection(windDir), m_windDirectionSource(windDirSrc), m_windDirectionError(windDirErr),
    m_validFrom(validFrom), m_validTo(validTo), m_location(lat, lon, alt)
{
}

WindField& WindField::operator=(const WindField& wf2)
{
    this->m_windDirection = wf2.m_windDirection;
    this->m_windDirectionError = wf2.m_windDirectionError;
    this->m_windDirectionSource = wf2.m_windDirectionSource;

    this->m_windSpeed = wf2.m_windSpeed;
    this->m_windSpeedError = wf2.m_windSpeedError;
    this->m_windSpeedSource = wf2.m_windSpeedSource;

    this->m_validFrom = wf2.m_validFrom;
    this->m_validTo = wf2.m_validTo;

    this->m_location = wf2.m_location;

    return *this;
}

void WindField::SetWindSpeed(double ws, MET_SOURCE source)
{
    this->m_windSpeed = ws;
    this->m_windSpeedSource = source;
}

void WindField::SetWindDirection(double wd, MET_SOURCE source)
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
MET_SOURCE WindField::GetWindSpeedSource() const
{
    return this->m_windSpeedSource;
}

/** Gets the source of the wind-speed */
void WindField::GetWindSpeedSource(novac::CString& str) const
{
    return Meteorology::MetSourceToString(m_windSpeedSource, str);
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
MET_SOURCE WindField::GetWindDirectionSource() const
{
    return this->m_windDirectionSource;
}

/** Gets the source of the wind-direction */
void WindField::GetWindDirectionSource(novac::CString& str) const
{
    return Meteorology::MetSourceToString(m_windDirectionSource, str);
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
