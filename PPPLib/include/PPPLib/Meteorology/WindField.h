#pragma once

#include <PPPLib/MFC/CString.h>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/GPSData.h>
#include <PPPLib/Meteorology/MeteorologySource.h>

#ifndef WINDFIELD_H
#define WINDFIELD_H

namespace Meteorology {

    /** global function that converts a MET_SOURCE item to string.
        This is used e.g. to write the source to file.
     */
    void MetSourceToString(const MET_SOURCE src, novac::CString& str);

    /** Global function that converts a string to a MET_SOURCE.
        This is basically the inverse function of 'MetSourceToString'.
        If the string does not evaluate to a known MET_SOURCE then
        MET_NONE is returned.
        This function is used e.g. to parse a log-file.
    */
    MET_SOURCE StringToMetSource(const novac::CString& str);

    /** Global function that retrieves the judged quality of a given source.
        Not all sources are equally good. A source with a higher 'quality'
        is judged to be better than a source with lower 'quality'.
        For instance will the source 'MET_GEOMETRY_CALCULATION' give a higher
        'quality' than 'MET_ECMWF_FORECAST' or 'MET_DEFAULT'.
        @return - an integer >= 0 and < MET_NUMBER_OF_DEFINED_SOURCES

     */
    int GetSourceQuality(MET_SOURCE src);


    /** The class <b>CWindField</b> is intended to hold information about the
        wind speed and direction at a given location and at a given point in time.
        It should be used as a container for this kind of data.
    */
    class CWindField
    {
    public:
        /** Default constructor */
        CWindField(void);

        /** Constructor */
        CWindField(double windSpeed, MET_SOURCE windSpeedSrc, double windDir, MET_SOURCE windDirSrc, const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double lat, double lon, double alt);

        /** Constructor */
        CWindField(double windSpeed, double windSpeedErr, MET_SOURCE windSpeedSrc, double windDir, double windDirErr, MET_SOURCE windDirSrc, const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double lat, double lon, double alt);

        /** Default destructor */
        ~CWindField(void);

        /** assignment operator */
        CWindField& operator=(const CWindField& wf2);

        /** Sets the wind-speed */
        void SetWindSpeed(double ws, MET_SOURCE source);

        /** Sets the wind-direction */
        void SetWindDirection(double wd, MET_SOURCE source);

        /** Sets the time-frame the wind-field is valid for */
        void SetValidTimeFrame(const novac::CDateTime& from, const novac::CDateTime& to);

        /** Gets the wind-speed */
        double GetWindSpeed() const;

        /** Gets the estimate for the total error in the wind-speed */
        double	GetWindSpeedError() const;

        /** Sets the estimate for the total error in the wind-speed */
        void	SetWindSpeedError(double err);

        /** Gets the source of the wind-speed */
        MET_SOURCE GetWindSpeedSource() const;

        /** Gets the source of the wind-speed */
        void GetWindSpeedSource(novac::CString& str) const;

        /** Gets the wind-direction */
        double GetWindDirection() const;

        /** Gets the source of the wind-direction */
        MET_SOURCE GetWindDirectionSource() const;

        /** Gets the source of the wind-direction */
        void GetWindDirectionSource(novac::CString& str) const;

        /** Gets the estimate for the total error in the wind-direction */
        double	GetWindDirectionError() const;

        /** Sets the estimate for the total error in the wind-direction */
        void	SetWindDirectionError(double err);

        /** Gets the time and date for which this wind-field is valid */
        void GetValidTimeFrame(novac::CDateTime& from, novac::CDateTime& to) const;

        /** Gets the position that this wind-field is valid for */
        void GetValidPosition(double& lat, double& lon, double& alt) const;
        void GetValidPosition(float& lat, float& lon, float& alt) const;

        /** Sets the position that this wind-field is valid for */
        void SetValidPosition(const double& lat, const double& lon, const double& alt);

    protected:

        /** The speed of the wind */
        double			m_windSpeed;
        double			m_windSpeedError;
        MET_SOURCE		m_windSpeedSource;

        /** The direction of the wind */
        double			m_windDirection;
        double			m_windDirectionError;
        MET_SOURCE		m_windDirectionSource;

        /** The time frame during which this piece of wind information
            is valid. */
        novac::CDateTime		m_validFrom;
        novac::CDateTime		m_validTo;

        /** The place on earth where this wind-field comes from */
        novac::CGPSData		m_location;
    };
}
#endif