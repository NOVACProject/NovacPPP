#pragma once

#include <list>
#include <vector>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Definitions.h>
#include <PPPLib/Meteorology/MeteorologySource.h>
#include <PPPLib/Meteorology/WindField.h>

namespace Meteorology
{
enum class InterpolationMethod
{
    Exact,
    NearestNeighbour,
    Bilinear,
};

/** An instance of the class <b>CWindDataBase</b> can be used to
    keep track of the wind field in a given region, typically
    for a single volcano.

    The CWindDataBase can store the variation of the wind-speed and
    wind direction with time and for several positions and altitudes.

*/
class CWindDataBase
{
public:

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The name of this database. Typically this is the name of the volcano for
        which the database is valid. Used to identify the database */
    std::string m_dataBaseName;

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Retrieves the wind field at a given time and at a given location.
        The data point with the highest quality that is valid at the given time
            and location will be returned.

        @param time - the time for which the wind field should be retrieved
        @param location - the position and altitude for which the wind field
            should be retrieved.
        @param method - specifies how the wind-field should be interpolated (in space) from
            the data in the database.
        @param windField - will on successful return be filled with the information
            about the wind field at the requested time and location.
        @return true if the wind field could be retrieved, otherwise false.
     */
    bool GetWindField(const novac::CDateTime& time, const novac::CGPSData& location, InterpolationMethod method, WindField& windField) const;

    /** Inserts a wind field into the database */
    void InsertWindField(const WindField& windField);

    /** Inserts a wind-direction into the database.
        @param validFrom - the time from which the wind-direction is judged to be ok
        @param validTo - the time until which the wind-direction is judged to be ok
        @param wd - the actual wind-direction, in degrees from north, positive angles clock-wise
        @param wd_err - the judged error in the wind-direction. Also in degrees
        @param wd_src - the source of this piece of information
        @param location - the point in space where this piece of information is valid.
            If location is NULL then the wind-direction is assumed to be valid everywhere.
            The altitude in 'location ' will be ignored. */
    void InsertWindDirection(const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double wd, double wd_err, MeteorologySource wd_src, const novac::CGPSData* location);

    /** Inserts a wind-speed into the database.
        @param validFrom - the time from which the wind-direction is judged to be ok
        @param validTo - the time until which the wind-direction is judged to be ok
        @param windSpeed - the actual wind-speed, in meters per second.
        @param windSpeed_err - the judged error in the wind-speed. Also in meters per second.
        @param ws_src - the source of this piece of information
        @param location - the point in space where this piece of information is valid.
            If location is NULL then the wind-speed is assumed to be valid everywhere.
            The altitude in 'location ' will be ignored. */
    void InsertWindSpeed(const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double windSpeed, double windSpeed_err, MeteorologySource wd_src, const novac::CGPSData* location);

    /** Writes the contents of this database to file.
        @return 0 on success. */
    int WriteToFile(const novac::CString& fileName) const;

    /** Retrieves the size of the database */
    int GetDataBaseSize() const;

private:

    // the structure winddata is used to hold the information about the wind for a
    //  single location and a single altitude
    struct WindData
    {
    public:
        WindData() = default;
        WindData(const WindData& other) = default;
        WindData& operator=(const WindData& other) = default;

        bool HasWindSpeed() const { return ws != NOT_A_NUMBER; }
        bool HasWindDirection() const { return wd != NOT_A_NUMBER; }

        /** This is the position where this wind field is valid
            If equal to -1 then this is valid everywhere */
        int location = -1;

        double ws = 10.0;   // the wind-speed, in meters/second
        double ws_err = 10.0;   // the absolute error in wind-speed, in meters/second
        MeteorologySource ws_src = MeteorologySource::Default;    // the source for the wind-speed

        double wd = 0.0;    // the wind-direction, in degrees from north, positive clock-wise
        double wd_err = 360.0;  // the absolute error in wind-direction, in degrees
        MeteorologySource wd_src = MeteorologySource::Default;    // the source for the wind-direction
    };

    /** This is used to organise all data that is valid within a certain
        time frame. */
    struct WindInTime
    {
    public:
        WindInTime();
        WindInTime(const WindInTime& w);
        WindInTime& operator=(const WindInTime& w);

        novac::CDateTime validFrom = novac::CDateTime(0, 0, 0, 0, 0, 0);  // this wind-data is valid from this day and time
        novac::CDateTime validTo = novac::CDateTime(9999, 12, 31, 23, 59, 59);    // this wind-data is valid until this day and time
        std::list <WindData> windData; // the list of wind-datas
    };

    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    /** This is the database of wind information. Each item in the list
        holds the information of the wind for a single time frame.

        Each WindInTime object in the list MUST have an unique time frame.
        The time frames MUST be disjoint!!!!
        */
    std::list <WindInTime> m_dataBase;

    /** These are all the positions that we have in our database */
    std::vector <novac::CGPSData> m_locations;


    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------

    /** Returns the location that belongs to the given location index.
        index == -1 represents a wind field which is valid everywhere. */
    const novac::CGPSData& GetLocation(int index) const;

    /** Retrieves the index of a given location in our list of locations.
        If not found in the list, this will return -1 */
    int GetLocationIndex(double lat, double lon, double alt) const;
    int GetLocationIndex(const novac::CGPSData& gps) const;

    /** Inserts a location into the array of locations.
        @return the index of the newly inserted location.
        If this location already exists in the database then it's index
        will be returned and nothing inserted in the list.
         */
    int InsertLocation(double lat, double lon, double alt);
    int InsertLocation(const novac::CGPSData& gps);

    /** The implementations of the different (spatial) interpolation methods.
    */
    bool GetWindField_Exact(const novac::CDateTime& time, const novac::CGPSData& location, WindField& windField) const;
    bool GetWindField_Nearest(const novac::CDateTime& time, const novac::CGPSData& location, WindField& windField) const;
    bool GetWindField_Bilinear(const novac::CDateTime& time, const novac::CGPSData& location, WindField& windField) const;


};
}