/**
  Common.h is a general header file for storing definitions and constants
      which are used throughout the whole program.
*/

#ifndef COMMON_H
#define COMMON_H

#include <PPPLib/PPPLib.h>
#include <PPPLib/Logging.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/MFC/CList.h>
#include <PPPLib/Measurement.h>
#include <PPPLib/Configuration/InstrumentType.h>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

#include <cmath>

namespace Geometry
{
class CPlumeHeight;
}
namespace Meteorology
{
class CWindField;
}

// ---------------------------------------------------------------
// --------------------------- LOGGING ---------------------------
// ---------------------------------------------------------------

class PocoLogger : public ILogger
{
public:
    /** Logs an informative message to the log */
    virtual void Information(const std::string& message) override;

    /** Logs an error message to the log */
    virtual void Error(const std::string& message) override;
};


/** Update the top line of list box */
void UpdateMessage(const novac::CString& message);

// ---------------------------------------------------------------
// ---------------- DEFINED CONSTANTS ----------------------------
// ---------------------------------------------------------------


// The maximum number of references that can be fitted to a single spectrum
#define MAX_N_REFERENCES 10

//size of buffer to receive data
#define SIZEOFBUFFER 2048

/** Not initialized values are set to 'NOT_A_NUMBER' to indicate
    that a value is missing.
    This can e.g. be the case if only the wind-direction (and not the wind speed)
        is known by an element in the database (which can be the case if a the
        wind direction was calculated by combining two scans).
    */
#define NOT_A_NUMBER -9999.0


    // ----------------------------------------------------------------
    // ---------------- MATHEMATICAL CONSTANTS ------------------------
    // ----------------------------------------------------------------

    // converts degrees to radians
#define DEGREETORAD 0.017453 

// converts radians to degrees
#define RADTODEGREE 57.295791

// a quite familiar constant
#define TWO_PI 6.28318
#define HALF_PI 1.5708
#ifndef M_PI
#define M_PI 3.141592
#endif

// -----------------------------------------------------------------
// -------------------------- MESSAGES -----------------------------
// -----------------------------------------------------------------

// signals to the 'NovacMasterProgramView' -class that a scan has been sucessfully evaluated
//  also signals to the reevaluation dialog that a spectrum has been sucessfully evaluated
#define WM_EVAL_SUCCESS		        WM_USER + 12


// -------------------------------------------------------
// ---------------- CLASS COMMON.H -----------------------
// -------------------------------------------------------

/** The class <b>Common</b> contains misc. functions that we need but
    do not really fit anywhere else.
*/

class Common
{

public:
    Common();

    // --------------------------------------------------------------------
    // ------------------------- FILE -------------------------------------
    // --------------------------------------------------------------------

    /** Get file size in bytes.
        @param - The file name (including path)
        @return - The file size (in bytes)
    */
    static long RetrieveFileSize(novac::CString& fileName);

    /** Compares two files to see if their contents are the same */
    static bool AreIdenticalFiles(const novac::CString& fileName1, const novac::CString& fileName2);

    /** If there's a file with the given input name, then it will be renamed to
        PATH\\FILENAME_creationDate_creationTime.FILEENDING */
    static bool ArchiveFile(const novac::CString& fileName);

    // --------------------------------------------------------------------
    // ------------------------ SYSTEM FUNCTIONS  -------------------------
    // --------------------------------------------------------------------

    // --------------------------------------------------------------------
    // ------------------------- PATH -------------------------------------
    // --------------------------------------------------------------------

    /** Take out the file name from a long path containing both a directory name and a file name.
        The path separator can be either '/' or '\'.
        @param fileName path of the file, will be set to only contain the filename, without the path . */
    static void GetFileName(novac::CString& fileName);

    /** Take out the directory from a long path name.
        @param fileName - the complete path of the file */
    static void GetDirectory(novac::CString& fileName);

    /** Copies the file to the new file location */
    static void CopyFile(const novac::CString& oldName, const novac::CString& newName);

    /** m_exePath will be set to the path where the application resides,
        This is set once at application startup and should never be changed. */
    const novac::CString m_exePath;

    /** m_exeFileName will be set to the filename of the program.
        This is set once at application startup and should never be changed. */
    const novac::CString m_exeFileName;

    // --------------------------------------------------------------------
    // -------------------- SUN - FUNCTIONS -------------------------------
    // --------------------------------------------------------------------

    /** Retrieves the solar zenith angle (SZA) and the solar azimuth angle (SAZ)
            for the site specified by (lat, lon) and for the time given in gmtTime.
            Note that the returned angles are in degrees and that the specified
            time _must_ be GMT-time. */
    static RETURN_CODE GetSunPosition(const novac::CDateTime& gmtTime, double lat, double lon, double& SZA, double& SAZ);

    /** ????
        D  is the JulianDay counted from the 1st of January 2000 @ midnight
        Ra is the  Right Ascencion
        dec	is the declination
        EQT is the equation of time (in hours)
    */
    static void EquatorialCoordinates(double D, double& RA, double& dec, double& EQT);

    /** ???
            lat - latitude in degrees
            H		- The hour angle, in degrees
            dec - The declination, degres
            elev - the returned elevation of the sun above the horizon
            azim - the returned azimuth angle of the sun, counted from south to west			*/
    static void HorizontalCoordinates(double lat, double H, double dec, double& elev, double& azim);

    /** Returns the hour angle given the longitude and equation of time. */
    static double GetHourAngle(double hr, double lon, double EqT);

    // --------------------------------------------------------------------
    // -------------------- GPS FUNCTIONS ---------------------------------
    // --------------------------------------------------------------------

    /* This function returns the distance in <b>meters</b> between the two points defined
        by (lat1,lon1) and (lat2, lon2). <b>All angles must be in degrees</b> */
    static double GPSDistance(double lat1, double lon1, double lat2, double lon2);

    /* This function returns the initial bearing (<b>degrees</b>) when travelling from
      the point defined by (lat1, lon1) to the point (lat2, lon2). <b>All angles must be in degrees</b> */
    static double GPSBearing(double lat1, double lon1, double lat2, double lon2);

    /** This function calculates the latitude and longitude for point
            which is the distance 'dist' m and bearing 'az' degrees from
            the point defied by 'lat1' and 'lon1' */
    void CalculateDestination(double lat1, double lon1, double dist, double az, double& lat2, double& lon2);

    // --------------------------------------------------------------------
    // -------------------- CALCULATING FLUX ------------------------------
    // --------------------------------------------------------------------

    /** Calculates the flux using the supplied data. Automatically decides which
            algorithm to use based on the given cone angle. */
    static double CalculateFlux(const double* scanAngle, const double* scanAngle2, const double* column, double offset, int nDataPoints, const Meteorology::CWindField& wind, const Geometry::CPlumeHeight& relativePlumeHeight, double compass, INSTRUMENT_TYPE type, double coneAngle = 90.0, double tilt = 0.0);

    // ---------------------------- MISC ----------------------------------

    /** Guesses the name of the specie from the name of the reference file. */
    static void GuessSpecieName(const novac::CString& string, novac::CString& specie);

};// end of class Common

// --------------------------------------------------------------------
// --------------- TEMPLATE ARRAY FUNCTIONS ---------------------------
// --------------------------------------------------------------------


/** Searches for the minimum element in the array.
        @param pBuffer - The array in which to search for an element.
        @param bufLen - The length of the array.
        @return - The minimum value in the array */
template <class T> T Min(T* pBuffer, long bufLen)
{
    T minValue = pBuffer[0];
    for (long i = 1; i < bufLen; i++)
    {
        if (pBuffer[i] < minValue)
            minValue = pBuffer[i];
    }
    return minValue;
}

/** This function finds the 'N' lowest values in the supplied array.
        On successfull return the array 'output' will be filled with the N lowest
        values in the input array, sorted in ascending order.
        @param array - the array to look into
        @param nElements - the number of elements in the supplied array
        @param output - the output array, must be at least 'N'- elements long
        @param N - the number of values to take out. 	*/
template <class T> bool FindNLowest(const T array[], long nElements, T output[], int N, int* indices = NULL)
{
    for (int i = 0; i < N; ++i)
        output[i] = 1e16; // to get some initial value

    // loop through all elements in the array
    for (int i = 0; i < nElements; ++i)
    {

        // compare this element with all elements in the output array.
        for (int j = 0; j < N; ++j)
        {
            if (array[i] < output[j])
            {
                // If we found a higher value, shift all other values down one step...
                for (int k = N - 1; k > j; --k)
                {
                    output[k] = output[k - 1];
                    if (indices)							indices[k] = indices[k - 1];
                }
                output[j] = array[i];
                if (indices)								indices[j] = i;
                break;
            }
        }
    }
    return true;
}

/** This function calculates the average value of all the elements
        in the supplied array.
        @param array - the array of which to calculate the average value
        @param nElements - the length of the array. */
template <class T> double Average(T array[], long nElements)
{
    if (nElements <= 0)
        return 0.0;

    double sum = 0;
    for (int k = 0; k < nElements; ++k)
    {
        sum += array[k];
    }
    return (sum / nElements);
}

/** This function calculates the variance of all the elements
        in the supplied array.
        @param array - the array of which to calculate the average value
        @param nElements - the length of the array. */
template <class T> double Variance(T array[], long nElements)
{
    if (nElements <= 0)
        return 0.0;

    // First get the mean-value
    T mean = Average(array, nElements);

    double sum = 0;
    for (int k = 0; k < nElements; ++k)
    {
        sum += (array[k] - mean) * (array[k] - mean);
    }
    sum = sum / nElements;
    return sum;
}

/** This function calculates the standard deviation of all the elements
        in the supplied array.
        @param array - the array of which to calculate the average value
        @param nElements - the length of the array. */
template <class T> double Std(T array[], long nElements)
{
    if (nElements <= 0)
        return 0.0;

    return std::sqrt(Variance(array, nElements));
}

#endif