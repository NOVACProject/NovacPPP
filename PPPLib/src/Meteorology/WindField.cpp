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

CWindField::CWindField(void)
{
    this->m_windDirection = 0.0;
    this->m_windDirectionSource = MET_DEFAULT;
    this->m_windDirectionError = 360.0;

    this->m_windSpeed = 10.0;
    this->m_windSpeedSource = MET_DEFAULT;
    this->m_windSpeedError = 10.0;

    m_validFrom = novac::CDateTime(2005, 10, 01, 00, 00, 00);
    m_validTo = novac::CDateTime(9999, 12, 31, 23, 59, 59);

    m_location = novac::CGPSData(0.0, 0.0, 0.0);
}

CWindField::CWindField(double windSpeed, MET_SOURCE windSpeedSrc, double windDir, MET_SOURCE windDirSrc, const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double lat, double lon, double alt)
{
    this->m_windSpeed = windSpeed;
    this->m_windSpeedSource = windSpeedSrc;
    this->m_windDirection = windDir;
    this->m_windDirectionSource = windDirSrc;
    this->m_validFrom = validFrom;
    this->m_validTo = validTo;
    this->m_location.m_latitude = lat;
    this->m_location.m_longitude = lon;
    this->m_location.m_altitude = alt;

    this->m_windSpeedError = 0.0;
    this->m_windDirectionError = 0.0;
}

CWindField::CWindField(double windSpeed, double windSpeedErr, MET_SOURCE windSpeedSrc, double windDir, double windDirErr, MET_SOURCE windDirSrc, const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double lat, double lon, double alt)
{
    this->m_windSpeed = windSpeed;
    this->m_windSpeedSource = windSpeedSrc;
    this->m_windSpeedError = windSpeedErr;

    this->m_windDirection = windDir;
    this->m_windDirectionSource = windDirSrc;
    this->m_windDirectionError = windDirErr;

    this->m_validFrom = validFrom;
    this->m_validTo = validTo;
    this->m_location.m_latitude = lat;
    this->m_location.m_longitude = lon;
    this->m_location.m_altitude = alt;
}


CWindField::~CWindField(void)
{
}

/** assignment operator */
CWindField& CWindField::operator=(const CWindField& wf2)
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

/** Sets the wind-speed */
void CWindField::SetWindSpeed(double ws, MET_SOURCE source)
{
    this->m_windSpeed = ws;
    this->m_windSpeedSource = source;
}

/** Sets the wind-direction */
void CWindField::SetWindDirection(double wd, MET_SOURCE source)
{
    this->m_windDirection = wd;
    this->m_windDirectionSource = source;

}

/** Sets the time and/or date the wind-field is valid for */
void CWindField::SetValidTimeFrame(const novac::CDateTime& from, const novac::CDateTime& to)
{
    this->m_validFrom = from;
    this->m_validTo = to;
}

/** Gets the wind-speed */
double CWindField::GetWindSpeed() const
{
    return this->m_windSpeed;
}

/** Gets the source of the wind-speed */
MET_SOURCE CWindField::GetWindSpeedSource() const
{
    return this->m_windSpeedSource;
}

/** Gets the source of the wind-speed */
void CWindField::GetWindSpeedSource(novac::CString& str) const
{
    return Meteorology::MetSourceToString(m_windSpeedSource, str);
}

/** Gets the estimate for the total error in the wind-speed */
double CWindField::GetWindSpeedError() const
{
    return this->m_windSpeedError;
}

/** Sets the estimate for the total error in the wind-speed */
void CWindField::SetWindSpeedError(double err)
{
    this->m_windSpeedError = err;
}


/** Gets the wind-direction */
double CWindField::GetWindDirection() const
{
    return this->m_windDirection;
}

/** Gets the source of the wind-direction */
MET_SOURCE CWindField::GetWindDirectionSource() const
{
    return this->m_windDirectionSource;
}

/** Gets the source of the wind-direction */
void CWindField::GetWindDirectionSource(novac::CString& str) const
{
    return Meteorology::MetSourceToString(m_windDirectionSource, str);
}

/** Gets the estimate for the total error in the wind-direction */
double CWindField::GetWindDirectionError() const
{
    return this->m_windDirectionError;
}

/** Sets the estimate for the total error in the wind-direction */
void CWindField::SetWindDirectionError(double err)
{
    this->m_windDirectionError = err;
}

/** Gets the time and date for which this wind-field is valid */
void CWindField::GetValidTimeFrame(novac::CDateTime& from, novac::CDateTime& to) const
{
    from = this->m_validFrom;
    to = this->m_validTo;
}

/** Gets the position that this wind-field is valid for */
void CWindField::GetValidPosition(double& lat, double& lon, double& alt) const
{
    lat = this->m_location.m_latitude;
    lon = this->m_location.m_longitude;
    alt = this->m_location.m_altitude;
}
/** Gets the position that this wind-field is valid for */
void CWindField::GetValidPosition(float& lat, float& lon, float& alt) const
{
    lat = (float)this->m_location.m_latitude;
    lon = (float)this->m_location.m_longitude;
    alt = (float)this->m_location.m_altitude;
}

/** Sets the position that this wind-field is valid for */
void CWindField::SetValidPosition(const double& lat, const double& lon, const double& alt)
{
    this->m_location.m_latitude = lat;
    this->m_location.m_longitude = lon;
    this->m_location.m_altitude = alt;
}
