#include "Common.h"

#include <algorithm>
#include <iostream>
#include <time.h>

#include <Poco/File.h>
#include <Poco/DirectoryIterator.h>
#include <Poco/DateTime.h>
#include <Poco/Message.h>
#include <Poco/Logger.h>

// include the global settings
#include <PPPLib/VolcanoInfo.h>

#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Flux/Flux.h>

#include "../Geometry/PlumeHeight.h"
#include <PPPLib/Meteorology/WindField.h>

extern novac::CVolcanoInfo g_volcanoes; // <-- the list of volcanoes

extern std::string s_exePath;
extern std::string s_exeFileName;

#undef min
#undef max

void UpdateMessage(const novac::CString& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message.std_str());
}

void ShowMessage(const novac::CString& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message.std_str());
}
void ShowMessage(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message);
}
void ShowMessage(const novac::CString& message, novac::CString connectionID)
{
    novac::CString msg;

    msg.Format("<%s> : %s", (const char*)connectionID, (const char*)message);

    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(msg.std_str());
}

void ShowMessage(const char message[])
{
    novac::CString msg;
    msg.Format("%s", message);
    ShowMessage(msg);
}

void ShowError(const novac::CString& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.fatal(message.std_str());
}
void ShowError(const char message[])
{
    novac::CString msg;
    msg.Format("%s", message);
    ShowError(msg);
}

void PocoLogger::Information(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message);
}

void PocoLogger::Error(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.fatal(message);
}

Common::Common()
    :m_exePath(s_exePath), m_exeFileName(s_exeFileName)
{

}

/** Calculate the distance (in meters) between the two points (lat1, lon1) and
    (lat2, lon2). All latitudes and longitudes should be in degrees. */
double Common::GPSDistance(double lat1, double lon1, double lat2, double lon2)
{
    const double R_Earth = 6367000; // radius of the earth
    lat1 = lat1 * DEGREETORAD;
    lat2 = lat2 * DEGREETORAD;
    lon1 = lon1 * DEGREETORAD;
    lon2 = lon2 * DEGREETORAD;

    double dLon = lon2 - lon1;
    double dLat = lat2 - lat1;

    if ((dLon == 0) && (dLat == 0))
        return 0;

    double a = std::pow((std::sin(dLat / 2.0)), 2.0) + std::cos(lat1) * std::cos(lat2) * std::pow((std::sin(dLon / 2.0)), 2.0);
    double c = 2 * std::asin(std::min(1.0, std::sqrt(a)));
    double distance = R_Earth * c;

    return distance;
}

/**count the angle from wind to north,also the plume direction compared with north
* the direction is from plume center to the source of the plume
* return degree value
*@lat1 - the latitude of beginning point or plume source, rad
*@lon1 - the longitude of beginning point orplume source,rad
*@lat2   - the latitude of ending point or plume center,rad
*@lon2   - the longitude of ending point or plume center,rad
*/
double Common::GPSBearing(double lat1, double lon1, double lat2, double lon2)
{
    lat1 = lat1 * DEGREETORAD;
    lat2 = lat2 * DEGREETORAD;
    lon1 = lon1 * DEGREETORAD;
    lon2 = lon2 * DEGREETORAD;
    double tmpAngle;
    double dLat = lat1 - lat2;
    double dLon = lon1 - lon2;

    if ((dLon == 0) && (dLat == 0))
        return 0;

    tmpAngle = atan2(-sin(dLon) * cos(lat2),
        cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon));

    /*  	tmpAngle = atan2(lon1*cos(lat1)-lon2*cos(lat2), lat1-lat2); */

    if (tmpAngle < 0)
    {
        tmpAngle = TWO_PI + tmpAngle;
    }

    tmpAngle = RADTODEGREE * tmpAngle;
    return tmpAngle;
}

/** This function calculates the latitude and longitude for a point
        which is the distance 'dist' m and bearing 'az' degrees from
        the point defied by 'lat1' and 'lon1' */
void Common::CalculateDestination(double lat1, double lon1, double dist, double az, double& lat2, double& lon2)
{
    const double R_Earth = 6367000; // radius of the earth

    double dR = dist / R_Earth;

    // convert to radians
    lat1 = lat1 * DEGREETORAD;
    lon1 = lon1 * DEGREETORAD;
    az = az * DEGREETORAD;

    // calculate the second point
    lat2 = asin(sin(lat1) * cos(dR) + cos(lat1) * sin(dR) * cos(az));

    lon2 = lon1 + atan2(sin(az) * sin(dR) * cos(lat1), cos(dR) - sin(lat1) * sin(lat2));

    // convert back to degrees
    lat2 = lat2 * RADTODEGREE;
    lon2 = lon2 * RADTODEGREE;
}

/** Retrieves the solar zenith angle (SZA) and the solar azimuth angle (SAZ)
        for the site specified by (lat, lon) and for the time given in gmtTime.
        Note that the returned angles are in degrees and that the specified
        time _must_ be GMT-time. */
RETURN_CODE Common::GetSunPosition(const novac::CDateTime& gmtTime, double lat, double lon, double& SZA, double& SAZ)
{
    SZA = SAZ = 0; // reset the numbers

    // Get the julian day
    double D = JulianDay(gmtTime) - 2451545.0;

    // Get the Equatorial coordinates...
    double	RA; //	the right ascension (deg)
    double	dec; // the declination	(deg)
    double	EqT;	// the equation of time (hours)
    EquatorialCoordinates(D, RA, dec, EqT);

    // Get the hour angle
    double fractionalHour = (double)gmtTime.hour + gmtTime.minute / 60.0 + gmtTime.second / 3600.0;
    double H = GetHourAngle(fractionalHour, lon, EqT);

    // Get the horizontal coordinates
    double	elev, sAzim; // The elevation and azimuth (towards south);
    HorizontalCoordinates(lat, H, dec, elev, sAzim);

    // Convert the elevation into sza
    SZA = 90.0 - elev;

    // Convert the azimuth to a value counted from the north and 
    SAZ = fmod(180.0 + sAzim, 360.0);

    return RETURN_CODE::SUCCESS;
}

double Common::CalculateFlux(const double* scanAngle, const double* scanAngle2, const double* column, double offset, int nDataPoints, const Meteorology::CWindField& wind, const Geometry::CPlumeHeight& relativePlumeHeight, double compass, INSTRUMENT_TYPE type, double coneAngle, double tilt)
{
    double windSpeed = wind.GetWindSpeed();
    double windDirection = wind.GetWindDirection();
    double plumeHeight = relativePlumeHeight.m_plumeAltitude;

    if (type == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
    {
        return CalculateFluxHeidelbergScanner(scanAngle, scanAngle2, column, offset, nDataPoints, windSpeed, windDirection, plumeHeight, compass);
    }
    else if (type == INSTRUMENT_TYPE::INSTR_GOTHENBURG)
    {
        // In the NovacPPP, the gas factor isn't used. However the flux-calculation formula, shared with the NovacProgram, requires the gas factor.
        //  This compensation factor is used to compensate for how the gas factor is weighted into the calculation...
        const double gasFactorCompensation = 1e6;
        if (fabs(coneAngle - 90.0) < 1.0)
        {
            return CalculateFluxFlatScanner(scanAngle, column, offset, nDataPoints, windSpeed, windDirection, plumeHeight, compass, gasFactorCompensation);
        }
        else
        {
            return CalculateFluxConicalScanner(scanAngle, column, offset, nDataPoints, windSpeed, windDirection, plumeHeight, compass, coneAngle, tilt, gasFactorCompensation);
        }
    }
    else
    {
        return 0.0; // unsupported instrument-type
    }
}


void Common::GuessSpecieName(const novac::CString& fileName, novac::CString& specie)
{
    specie.Format("");
    novac::CString spc[] = { "SO2", "NO2", "O3", "O4", "HCHO", "RING", "H2O", "CLO", "BRO", "CHOCHO", "Glyoxal", "Formaldehyde", "HONO", "NO3" };
    int nSpecies = 12;

    int index = fileName.ReverseFind('\\');
    if (index == 0)
        return;

    novac::CString fil;
    fil.Format("%s", (const char*)fileName.Right((int)strlen(fileName) - index - 1));
    fil.MakeUpper();

    for (int i = 0; i < nSpecies; ++i)
    {
        if (strstr(fil, spc[i]))
        {
            specie.Format("%s", (const char*)spc[i]);
            return;
        }
    }

    // nothing found
    return;
}

void Common::GetFileName(novac::CString& fileName)
{
    // look for slashes in the path
    int position = std::max(fileName.ReverseFind('\\'), fileName.ReverseFind('/'));
    int length = fileName.GetLength();
    fileName = fileName.Right(length - position - 1);
}

/** Take out the directory from a long path name.
    @param fileName - the complete path of the file */
void Common::GetDirectory(novac::CString& fileName)
{
    int position = fileName.ReverseFind('\\');
    if (position >= 0)
    {
        fileName = fileName.Left(position + 1);
    }
}

void Common::CopyFile(const novac::CString& oldName, const novac::CString& newName)
{
    Poco::File oldFile(oldName.std_str());

    oldFile.copyTo(newName.std_str());
}

long Common::RetrieveFileSize(novac::CString& fileName)
{
    Poco::File file(fileName.std_str());
    return file.getSize();
}


/** Compares two files to see if their contents are the same */
bool Common::AreIdenticalFiles(const novac::CString& fileName1, const novac::CString& fileName2)
{
    if (Equals(fileName1, fileName2))
        return true; // a file is always identical to itself

    FILE* f1 = fopen(fileName1, "r");
    if (f1 == NULL)
        return false;

    FILE* f2 = fopen(fileName2, "r");
    if (f2 == NULL)
        return false;

    while (1)
    {
        int c1 = fgetc(f1);
        int c2 = fgetc(f2);

        if (c1 == EOF)
        {
            fclose(f1); fclose(f2);
            if (c2 == EOF)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        if (c1 != c2)
            return false;
    }


    fclose(f1); fclose(f2);

    // should never reach this point...
    return false;
}

/** If there's a file with the given input name, then it will be renamed to
    PATH\\FILENAME_creationDate_creationTime.FILEENDING */
bool Common::ArchiveFile(const novac::CString& fileName)
{
    novac::CString newFileName, errorMsg;

    // Search for the file
    Poco::File oldFile(fileName.std_str());

    if (!oldFile.exists())
    {
        return false; // file does not exist
    }

    // Get the time the file was created
    // TODO: Is this the local time or the UTC time stamp??
    Poco::DateTime creationTime = Poco::DateTime(oldFile.created());

    // build the new file-name
    int lastDot = fileName.ReverseFind('.');
    if (lastDot == -1)
    {
        newFileName.Format("%s_%04d%02d%02d_%02d%02d", (const char*)fileName,
            creationTime.year(), creationTime.month(), creationTime.day(), creationTime.hour(), creationTime.minute());
    }
    else
    {
        newFileName.Format("%s_%04d%02d%02d_%02d%02d%s", (const char*)fileName.Left(lastDot),
            creationTime.year(), creationTime.month(), creationTime.day(), creationTime.hour(), creationTime.minute(), (const char*)fileName.Right(fileName.GetLength() - lastDot));
    }

    // move the file
    oldFile.moveTo(newFileName.std_str());

    return true;
}



/*EQUATORIAL COORDINATES:RIGHT ASCENSION AND DECLINATION*/
void Common::EquatorialCoordinates(double D, double& RA, double& dec, double& EQT)
{
    double g_deg, q_deg, L_deg;				/*ANGLES IN DEGREES*/
    double g_rad, q_rad, L_rad;				/*ANGLES IN	RADIANS*/
    double R;												/*DISTANCE SUN-EARTH IN A.U*/
    double obliq_deg, obliq_rad;			/*OBLIQUITY OF THE ECLIPTIC*/
    double RA_rad, dec_rad;					/*EQUATORIAL COORDINATES IN RADIANS*/

    g_deg = fmod(357.529 + 0.98560028 * D, 360.0);
    g_rad = g_deg * DEGREETORAD;
    q_deg = fmod(280.459 + 0.98564736 * D, 360.0);
    q_rad = q_deg * DEGREETORAD;

    L_deg = q_deg + 1.915 * sin(g_rad) + 0.02 * sin(2 * g_rad);
    L_rad = L_deg * DEGREETORAD;

    // The distance between the sun and the earth (in Astronomical Units)
    R = 1.00014 - 0.01671 * cos(g_rad) - 0.00014 * cos(2 * g_rad);

    // The obliquity of the earth's orbit:
    obliq_deg = 23.439 - 0.00000036 * D;
    obliq_rad = obliq_deg * DEGREETORAD;

    // The right ascension (RA)
    RA_rad = atan(cos(obliq_rad) * sin(L_rad) / cos(L_rad));
    if (RA_rad < 0)
        RA_rad = TWO_PI + RA_rad;

    if (fabs(RA_rad - L_rad) > 1.570796)
        RA_rad = M_PI + RA_rad;

    dec_rad = asin(sin(obliq_rad) * sin(L_rad));
    RA = fmod(RA_rad * RADTODEGREE, 360.0);		// The right ascension

    // The declination
    dec = dec_rad * RADTODEGREE;

    // The Equation of Time
    EQT = q_deg / 15.0 - RA / 15.0;
}

void Common::HorizontalCoordinates(double lat, double H, double dec, double& elev, double& azim)
{
    const double H_rad = H * DEGREETORAD;
    const double lat_rad = lat * DEGREETORAD;
    const double dec_rad = dec * DEGREETORAD;

    double azim_rad;

    // The elevation angle
    const double elev_rad = asin(cos(H_rad) * cos(dec_rad) * cos(lat_rad) + sin(dec_rad) * sin(lat_rad));

    // The cosine of the azimuth - angle
    const double cazim_rad = (cos(H_rad) * cos(dec_rad) * sin(lat_rad) - sin(dec_rad) * cos(lat_rad)) / cos(elev_rad);

    // The sine of the azimuth - angle
    double sazim_rad = (sin(H_rad) * cos(dec_rad)) / cos(elev_rad);

    if (cazim_rad > 0 && sazim_rad > 0)
        azim_rad = asin(sazim_rad);				// azim is in the range 0 - 90 degrees
    else if (cazim_rad < 0 && sazim_rad > 0)
        azim_rad = M_PI - asin(sazim_rad);		// azim is in the range 90 - 180 degrees
    else if (cazim_rad < 0 && sazim_rad < 0)
        azim_rad = M_PI - asin(sazim_rad);		// azim is in the range 180 - 270 degrees
    else if (cazim_rad > 0 && sazim_rad < 0)
        azim_rad = TWO_PI + asin(sazim_rad);	// azim is in the range 270 - 360 degrees
    else
        azim_rad = 0.0; // ERROR...

    elev = elev_rad * RADTODEGREE;
    azim = azim_rad * RADTODEGREE;

    //printf("\n\nHORIZONTAL COORDINATES:\nAZIMUTH (from South to West)=%f\nELEVATION=%f\n\n",*pazim,*pelev);
}

/** Returns the hour angle given the longitude and equation of time. */
double Common::GetHourAngle(double hr, double lon, double EqT)
{
    double H = 15.0 * (hr + lon / 15 + EqT - 12);
    //    printf("HOUR ANGLE (from noon,increasing with time): %f\n",H);
    return(H);
}

