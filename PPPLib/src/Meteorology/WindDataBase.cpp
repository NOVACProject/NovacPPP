#include <PPPLib/Meteorology/WindDataBase.h>
#include <PPPLib/Gps.h>
#include <math.h>

namespace Meteorology
{

// ----------- THE SUB-CLASS WindInTime --------------
CWindDataBase::WindInTime::WindInTime()
{}
CWindDataBase::WindInTime::WindInTime(const WindInTime& other)
{
    this->validFrom = other.validFrom;
    this->validTo = other.validTo;
    std::list<WindData>::const_iterator p = other.windData.begin();
    while (p != other.windData.end())
    {
        this->windData.push_back((WindData&)*(p++));
    }
}
CWindDataBase::WindInTime& CWindDataBase::WindInTime::operator=(const WindInTime& other)
{
    this->validFrom = other.validFrom;
    this->validTo = other.validTo;
    std::list<WindData>::const_iterator p = other.windData.begin();
    while (p != other.windData.end())
    {
        this->windData.push_back((WindData&)*(p++));
    }
    return *this;
}


// --------- THE CLASS CWindDataBase ----------


bool CWindDataBase::GetWindField(const novac::CDateTime& time, const novac::CGPSData& location, InterpolationMethod method, WindField& windField) const
{
    if (InterpolationMethod::Exact == method)
    {
        return GetWindField_Exact(time, location, windField);
    }
    else if (InterpolationMethod::NearestNeighbour == method)
    {
        return GetWindField_Nearest(time, location, windField);
    }
    else if (InterpolationMethod::Bilinear == method)
    {
        return GetWindField_Bilinear(time, location, windField);
    }

    return false; // nothing found in the database
}

void CWindDataBase::InsertWindField(const WindField& windField)
{
    bool foundMatchingTimeFrame = false;
    novac::CDateTime startTime, endTime;
    novac::CGPSData position;

    // Get the time-frame from the wind-field
    windField.GetValidTimeFrame(startTime, endTime);

    // Loop through the database and see if there is alreay an item with this time-frame
    std::list<WindInTime>::const_iterator pos = m_dataBase.begin();
    while (pos != m_dataBase.end())
    {
        WindInTime& t = (WindInTime&)*pos;

        // check if the timeframe for this item matches the wind-field to insert
        if (t.validFrom == startTime && t.validTo == endTime)
        {
            foundMatchingTimeFrame = true;
            break; // jump out of the while-loop
        }

        // increase the iterator...
        ++pos;
    }

    // create a new data object to insert
    WindData data;
    data.ws = windField.GetWindSpeed();
    data.ws_err = windField.GetWindSpeedError();
    data.ws_src = windField.GetWindSpeedSource();
    data.wd = windField.GetWindDirection();
    data.wd_err = windField.GetWindDirectionError();
    data.wd_src = windField.GetWindDirectionSource();
    windField.GetValidPosition(position.m_latitude, position.m_longitude, position.m_altitude);
    if (position.m_latitude == NOT_A_NUMBER && position.m_longitude == NOT_A_NUMBER)
    {
        data.location = -1;
    }
    else
    {
        data.location = InsertLocation(position);
    }

    // if we found a matching time frame then insert the wind at that location in the database
    if (foundMatchingTimeFrame)
    {
        WindInTime& t = (WindInTime&)*pos;

        // insert the new data object into the database and return
        t.windData.push_back(data);

        return;
    }
    else
    {
        // it's not found in the database. Insert it as a new item.
        WindInTime t;
        t.validFrom = startTime;
        t.validTo = endTime;
        t.windData.push_back(data);

        this->m_dataBase.push_back(t);

        return;
    }
}

/** Inserts a wind-direction into the database */
void CWindDataBase::InsertWindDirection(const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double wd, double wd_err, MeteorologySource wd_src, const novac::CGPSData* location)
{
    WindField windField;

    if (location != NULL)
    {
        windField = WindField(NOT_A_NUMBER, NOT_A_NUMBER, MeteorologySource::None, wd, wd_err, wd_src, validFrom, validTo, location->m_latitude, location->m_longitude, location->m_altitude);
        this->InsertWindField(windField);
    }
    else
    {
        windField = WindField(NOT_A_NUMBER, NOT_A_NUMBER, MeteorologySource::None, wd, wd_err, wd_src, validFrom, validTo, NOT_A_NUMBER, NOT_A_NUMBER, NOT_A_NUMBER);
        this->InsertWindField(windField);
    }
}

/** Inserts a wind-direction into the database */
void CWindDataBase::InsertWindSpeed(const novac::CDateTime& validFrom, const novac::CDateTime& validTo, double ws, double ws_err, MeteorologySource ws_src, const novac::CGPSData* location)
{
    WindField windField;

    if (location != NULL)
    {
        windField = WindField(ws, ws_err, ws_src, NOT_A_NUMBER, NOT_A_NUMBER, MeteorologySource::None, validFrom, validTo, location->m_latitude, location->m_longitude, location->m_altitude);
        this->InsertWindField(windField);
    }
    else
    {
        windField = WindField(ws, ws_err, ws_src, NOT_A_NUMBER, NOT_A_NUMBER, MeteorologySource::None, validFrom, validTo, NOT_A_NUMBER, NOT_A_NUMBER, NOT_A_NUMBER);
        this->InsertWindField(windField);
    }
}

int CWindDataBase::WriteToFile(const novac::CString& fileName) const
{
    novac::CString indent, sourceStr;

    // open the file
    FILE* f = fopen(fileName, "w");
    if (f == NULL)
    {
        return 1;
    }

    // write the header lines and the start of the file
    fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(f, "<!-- This file defines the wind field for a given volcano. To be used\n for the calculation of fluxes in the NOVAC Post Processing Program -->\n\n");
    fprintf(f, "<Wind volcano=\"%s\">\n", m_dataBaseName.c_str());
    indent.Format("\t");

    // loop through the list of "WindInTime's" and write them to file 
    std::list<WindInTime>::const_iterator pos = m_dataBase.begin();
    while (pos != m_dataBase.end())
    {
        // write the start of the <windfield> section
        fprintf(f, "%s<windfield>\n", (const char*)indent);

        // get the next element in the list
        const WindInTime& time = (const WindInTime&)*(pos++);

        // make sure that there's at least one item in this list...
        std::list<WindData>::const_iterator wPos = time.windData.begin();
        if (wPos != time.windData.end())
        {
            const novac::CDateTime& from = time.validFrom;
            const novac::CDateTime& to = time.validTo;
            const WindData& data = (const WindData&)*(wPos);

            if (data.wd == NOT_A_NUMBER)
            {
                Meteorology::MetSourceToString(data.ws_src, sourceStr);
            }
            else
            {
                Meteorology::MetSourceToString(data.wd_src, sourceStr);
            }

            // write this one to string
            fprintf(f, "\t%s<source>%s</source>\n", (const char*)indent, (const char*)sourceStr);
            fprintf(f, "\t%s<altitude>%.1f</altitude>\n", (const char*)indent, GetLocation(data.location).m_altitude);
            fprintf(f, "\t%s<valid_from>%04d.%02d.%02dT%02d:%02d:%02d</valid_from>\n", (const char*)indent, from.year, from.month, from.day, from.hour, from.minute, from.second);
            fprintf(f, "\t%s<valid_to>%04d.%02d.%02dT%02d:%02d:%02d</valid_to>\n", (const char*)indent, to.year, to.month, to.day, to.hour, to.minute, to.second);

            std::list<WindData>::const_iterator pos2 = time.windData.begin();
            // loop through each item in the list and write it down
            while (pos2 != time.windData.end())
            {
                const WindData& data2 = (const WindData&)*(pos2++);
                const novac::CGPSData& dataPos2 = GetLocation(data2.location);
                fprintf(f, "\t%s<item lat=\"%.2f\" lon=\"%.2f\" ws=\"%.2f\" wse=\"%.2f\" wd=\"%.2f\" wde=\"%.2f\"/>\n", (const char*)indent, dataPos2.m_latitude, dataPos2.m_longitude, data2.ws, data2.ws_err, data2.wd, data2.wd_err);
            }
        }

        // write the end of the <windfield> section
        fprintf(f, "%s</windfield>\n", (const char*)indent);
    }

    fprintf(f, "</Wind>\n");

    // remember to close the file when we're done
    fclose(f);

    return 0;
}


/** Returns the location that belongs to the given location index. */
const novac::CGPSData& CWindDataBase::GetLocation(int index) const
{
    static const novac::CGPSData nullPos = novac::CGPSData(-1, -1, -1);
    if ((index < 0) || (index >= (int)m_locations.size()))
        return nullPos;
    else
        return m_locations.at(index);
}

/** Retrieves the index of a given location in our list of locations.
    If not found in the list, this will return -1 */
int CWindDataBase::GetLocationIndex(double lat, double lon, double alt) const
{
    return GetLocationIndex(novac::CGPSData(lat, lon, alt));
}

int CWindDataBase::GetLocationIndex(const novac::CGPSData& gps) const
{
    int N = (int)m_locations.size();
    for (int k = 0; k < N; ++k)
    {
        const novac::CGPSData& gps2 = m_locations.at(k);
        if (gps == gps2)
        {
            return k;
        }
    }
    // not found in the list
    return -1;
}

/** Inserts a location into the array of locations.
    @return the index of the newly inserted location */
int CWindDataBase::InsertLocation(double lat, double lon, double alt)
{
    return InsertLocation(novac::CGPSData(lat, lon, alt));
}
int CWindDataBase::InsertLocation(const novac::CGPSData& gps)
{
    // check if this position is already in the list. If so then don't insert it...
    int locationIndex = GetLocationIndex(gps);
    if (locationIndex >= 0)
        return locationIndex;

    // this position does not already exist, add it...
    int N = (int)m_locations.size();
    m_locations.push_back(novac::CGPSData(gps));
    return N;
}


bool CWindDataBase::GetWindField_Exact(const novac::CDateTime& time, const novac::CGPSData& location, WindField& windField) const
{
    int bestWs_Quality = -1; // the best quality data of wind-speed that we found
    int bestWd_Quality = -1; // the best quality data of wind-speed that we found
    WindField tempWindField = WindField(NOT_A_NUMBER, MeteorologySource::None, NOT_A_NUMBER, MeteorologySource::None, novac::CDateTime(), novac::CDateTime(9999, 12, 31, 23, 59, 59), location.m_latitude, location.m_longitude, location.m_altitude);
    novac::CDateTime validFrom = novac::CDateTime(0, 0, 0, 0, 0, 0);
    novac::CDateTime validTo = novac::CDateTime(9999, 12, 31, 23, 59, 59);
    int ws_Average = 0; // how many data-points is the wind-speed an average of...
    int wd_Average = 0; // how many data-points is the wind-direction an average of...

    // Get the location index for this location
    int locationIndex = GetLocationIndex(location);
    if (locationIndex == -1)
    {
        // this point does not exist in our database...
        return false;
    }

    // search through the database to see if we can find any item that is valid for this time
    // Loop through the database and see if there is alreay an item with this time-frame
    std::list<WindInTime>::const_iterator pos_t = m_dataBase.begin();
    while (pos_t != m_dataBase.end())
    {
        const WindInTime& t = (const WindInTime&)*(pos_t++);

        // check if the given time matches this interval
        if ((t.validFrom > time) || (time > t.validTo))
            continue;

        // loop through all the data points at this time-step to extract the data point
        //	with the highest quality at this time
        std::list<WindData>::const_iterator pos = t.windData.begin();
        while (pos != t.windData.end())
        {
            const WindData& data = (const WindData&)*(pos++);

            // if this is not the right spot...
            if ((data.location != -1) && (data.location != locationIndex))
                continue;

            int ws_quality = GetSourceQuality(data.ws_src);
            int wd_quality = GetSourceQuality(data.wd_src);

            // ------- The wind-speed ---------
            if (ws_quality > bestWs_Quality)
            {
                // we found a better source than we already have
                //	replace the information that we have with the new one.
                bestWs_Quality = ws_quality;
                tempWindField.SetWindSpeed(data.ws, data.ws_src);
                tempWindField.SetWindSpeedError(data.ws_err * data.ws_err);
                ws_Average = 1;
                if (t.validFrom > validFrom)
                {
                    validFrom = t.validFrom;
                }
                if (t.validTo < validTo)
                {
                    validTo = t.validTo;
                }
            }
            else if (ws_quality == bestWs_Quality)
            {
                // we found data with the same quality as we already have
                //	make the information an average of the old and the new information
                tempWindField.SetWindSpeed(tempWindField.GetWindSpeed() + data.ws, data.ws_src);
                tempWindField.SetWindSpeedError(data.ws_err * data.ws_err + tempWindField.GetWindSpeedError());
                ++ws_Average;
                if (t.validFrom > validFrom)
                {
                    validFrom = t.validFrom;
                }
                if (t.validTo < validTo)
                {
                    validTo = t.validTo;
                }
            }

            // ------- The wind direction ------
            if (wd_quality > bestWd_Quality)
            {
                // we found a better source than we already have
                //	replace the information that we have with the new one.
                bestWd_Quality = wd_quality;
                tempWindField.SetWindDirection(data.wd, data.wd_src);
                tempWindField.SetWindDirectionError(data.wd_err * data.wd_err);
                wd_Average = 1;
                if (t.validFrom > validFrom)
                {
                    validFrom = t.validFrom;
                }
                if (t.validTo < validTo)
                {
                    validTo = t.validTo;
                }
            }
            else if (wd_quality == bestWd_Quality)
            {
                // we found data with the same quality as we already have
                //	make the information an average of the old and the new information
                tempWindField.SetWindDirection(tempWindField.GetWindDirection() + data.wd, data.wd_src);
                tempWindField.SetWindDirectionError(data.wd_err * data.wd_err + tempWindField.GetWindDirectionError());
                ++wd_Average;
                if (t.validFrom > validFrom)
                {
                    validFrom = t.validFrom;
                }
                if (t.validTo < validTo)
                {
                    validTo = t.validTo;
                }
            }
        }
    }

    if (bestWs_Quality <= GetSourceQuality(MeteorologySource::None) || bestWd_Quality <= GetSourceQuality(MeteorologySource::None))
        return false; // no matching location found.
    else
    {
        // make the wind-speeds and wind-direction averages...
        double avgWindSpeed = tempWindField.GetWindSpeed() / ws_Average;
        double avgWindDir = tempWindField.GetWindDirection() / wd_Average;
        double wsErr = sqrt(tempWindField.GetWindSpeedError());
        double wdErr = sqrt(tempWindField.GetWindDirectionError());

        windField.SetWindSpeed(avgWindSpeed, tempWindField.GetWindSpeedSource());
        windField.SetWindDirection(avgWindDir, tempWindField.GetWindDirectionSource());
        windField.SetWindSpeedError(wsErr);
        windField.SetWindDirectionError(wdErr);
        windField.SetValidTimeFrame(validFrom, validTo);
        windField.m_location = location;
        return true;
    }
}

// This function takes the wind-field in the nearest datapoint in the database
bool CWindDataBase::GetWindField_Nearest(const novac::CDateTime& time, const novac::CGPSData& location, WindField& windField) const
{
    double smallestDistance = 1e99; // the smallest distance from a point in the database to 'location'
    int closestPoint = -1; // the index of the closest location

    // loop through all locations to see which one is the closest
    for (unsigned int k = 0; k < m_locations.size(); ++k)
    {
        const novac::CGPSData& pos = GetLocation(k);

        // compare the position with the given one
        double dist = Gps::GpsMath::GPSDistance(location.m_latitude, location.m_longitude, pos.m_latitude, pos.m_longitude);
        if (dist < smallestDistance)
        {
            closestPoint = k;
            smallestDistance = dist;
        }
    }
    if (closestPoint == -1)
    {
        return false; // no point found.
    }

    // return the wind field at the closest point
    const novac::CGPSData& closestGPSPoint = GetLocation(closestPoint);
    return GetWindField_Exact(time, closestGPSPoint, windField);
}

// This function calculates the wind-field as a bi-linear interpolation of
//	the wind-field in the nearest four datapoints in the database
//	This assumes that the grid is regular
bool CWindDataBase::GetWindField_Bilinear(const novac::CDateTime& /*time*/, const novac::CGPSData& /*location*/, WindField& /*windField*/) const
{

    //// First make a reasonability check to make sure that the database is ok with this
    //if(time.windData.GetCount() < 4)
    //	return false;

    //// 1. -------- Find four data-points that encloses the given location ---------
    //WindField wf11, wf12, wf21, wf22;
    //double dist11 = 1e99, dist12 = 1e99, dist21 = 1e99, dist22 = 1e99;

    //// loop through all locations to find the four closest ones
    //POSITION pos = time.windData.GetHeadPosition();
    //while(pos != NULL){
    //	const WindData &data = time.windData.GetNext(pos);
    //	const novac::CGPSData &dataPos = GetLocation(data.location);

    //	// compare the position with the given one
    //	double dist = Common::GPSDistance(location.m_latitude, location.m_longitude, dataPos.m_latitude, dataPos.m_longitude);
    //	
    //	// Sort out if the point is in the south-west (11),
    //	//	south-east(21), north-west(12) or the north-east(22)
    //	if(dataPos.m_latitude <= location.m_latitude){
    //		// wf.1
    //		if(dataPos.m_longitude <= location.m_longitude){
    //			// wf11
    //			if(dist < dist11){
    //				wf11	= WindField(data.ws, data.ws_err, data.ws_src, data.wd, data.wd_err, data.wd_src, time.validFrom, time.validTo, dataPos.m_latitude, dataPos.m_longitude, dataPos.m_altitude);
    //				dist11	= dist;
    //			}
    //		}else{
    //			// wf21
    //			if(dist < dist21){
    //				wf21	= WindField(data.ws, data.ws_err, data.ws_src, data.wd, data.wd_err, data.wd_src, time.validFrom, time.validTo, dataPos.m_latitude, dataPos.m_longitude, dataPos.m_altitude);
    //				dist21	= dist;
    //			}
    //		}
    //	}else{
    //		// wf.2
    //		if(dataPos.m_longitude <= location.m_longitude){
    //			// wf12
    //			if(dist < dist12){
    //				wf12	= WindField(data.ws, data.ws_err, data.ws_src, data.wd, data.wd_err, data.wd_src, time.validFrom, time.validTo, dataPos.m_latitude, dataPos.m_longitude, dataPos.m_altitude);
    //				dist12	= dist;
    //			}
    //		}else{
    //			// wf22
    //			if(dist < dist22){
    //				wf22	= WindField(data.ws, data.ws_err, data.ws_src, data.wd, data.wd_err, data.wd_src, time.validFrom, time.validTo, dataPos.m_latitude, dataPos.m_longitude, dataPos.m_altitude);
    //				dist22	= dist;
    //			}
    //		}
    //	}//end if
    //} // endwhile
    //
    //// 2.  --------------- Check the data ----------------

    //// this algorithm assumes that the grid is regular. Make sure that this is the case.
    ////  start by extracting the four corners of the rectangle
    //double x1, x2, y1, y2, trams;
    //wf11.GetValidPosition(y1, x1, trams);
    //wf22.GetValidPosition(y2, x2, trams);
    //
    //// now check that each data point has it's lat equal to x1 or x2 and it's lon equal to y1 or y2.
    //double tmpLat, tmpLon, tmpAlt;
    //wf12.GetValidPosition(tmpLat, tmpLon, tmpAlt);
    //if(fabs(tmpLat - y2) > 1e-4 || fabs(tmpLon - x1) > 1e-4)
    //	return false;
    //wf21.GetValidPosition(tmpLat, tmpLon, tmpAlt);
    //if(fabs(tmpLat - y1) > 1e-4 || fabs(tmpLon - x2) > 1e-4)
    //	return false;
    //
    //// 3. -------- Make the interpolation ---------
    //
    //// to make this, we need to extract the u and v components and interpolate them
    ////	separately
    //double u11, u12, u21, u22;
    //double v11, v12, v21, v22;
    //int nFound = 0;
    //u11 = wf11.GetWindSpeed() * cos(DEGREETORAD * wf11.GetWindDirection());
    //v11 = wf11.GetWindSpeed() * sin(DEGREETORAD * wf11.GetWindDirection());
    //u12 = wf12.GetWindSpeed() * cos(DEGREETORAD * wf12.GetWindDirection());
    //v12 = wf12.GetWindSpeed() * sin(DEGREETORAD * wf12.GetWindDirection());
    //u21 = wf21.GetWindSpeed() * cos(DEGREETORAD * wf21.GetWindDirection());
    //v21 = wf21.GetWindSpeed() * sin(DEGREETORAD * wf21.GetWindDirection());
    //u22 = wf22.GetWindSpeed() * cos(DEGREETORAD * wf22.GetWindDirection());
    //v22 = wf22.GetWindSpeed() * sin(DEGREETORAD * wf22.GetWindDirection());

    //// the lat and long of the point we want to extract
    //double x = location.m_longitude;
    //double y = location.m_latitude;

    ////  The interpolated u and v components (from http://en.wikipedia.org/wiki/Bilinear_interpolation)
    //double u_interp = u11 * (x2 - x)*(y2 - y) + u21 * (x - x1)*(y2 - y) + 
    //                  u12 * (x2 - x)*(y - y1) + u22 * (x - x1)*(y - y1);
    //u_interp /= (x2 - x1)*(y2 - y1);
    //double v_interp = v11 * (x2 - x)*(y2 - y) + v21 * (x - x1)*(y2 - y) + 
    //                  v12 * (x2 - x)*(y - y1) + v22 * (x - x1)*(y - y1);
    //v_interp /= (x2 - x1)*(y2 - y1);

    //// Finally put together the u and v to a wind speed and direction
    //double wd = RADTODEGREE * atan2(v_interp, u_interp);
    //double ws = sqrt(u_interp * u_interp + v_interp * v_interp);	

    //// extract the position
    //pos = time.windData.GetHeadPosition();
    //const WindData &data = time.windData.GetAt(pos);
    //const novac::CGPSData &dataPos = GetLocation(data.location);

    //// Construct the wind field
    //windField = WindField(ws, data.ws_err, data.ws_src, wd, data.wd_err, data.ws_src, time.validFrom, time.validTo, dataPos.m_latitude, dataPos.m_longitude, dataPos.m_altitude);

    //return true;
    return false;
}

int CWindDataBase::GetDataBaseSize() const
{
    return (int)m_dataBase.size();
}

}
