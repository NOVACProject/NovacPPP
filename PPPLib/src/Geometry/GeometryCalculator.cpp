#include <PPPLib/Geometry/GeometryCalculator.h>
#include <PPPLib/VolcanoInfo.h>
#include <PPPLib/File/EvaluationLogFileHandler.h>

#include <SpectralEvaluation/GPSData.h>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

#include <Poco/Path.h>
#include <algorithm>
#include <cmath>

#undef min
#undef max

using namespace Geometry;
using namespace novac;

extern novac::CVolcanoInfo					g_volcanoes;   // <-- A list of all known volcanoes
extern Configuration::CUserConfiguration	g_userSettings;// <-- The settings of the user

CGeometryCalculator::CGeometryCalculationInfo::CGeometryCalculationInfo()
{
    Clear();
}
CGeometryCalculator::CGeometryCalculationInfo::~CGeometryCalculationInfo()
{
}
void CGeometryCalculator::CGeometryCalculationInfo::Clear()
{
    for (int k = 0; k < 2; ++k)
    {
        scanner[k].m_altitude = 0;
        scanner[k].m_latitude = 0.0;
        scanner[k].m_longitude = 0.0;
        plumeCentre[k] = 0.0;
    }
}
Geometry::CGeometryCalculator::CGeometryCalculationInfo& CGeometryCalculator::CGeometryCalculationInfo::operator =(const Geometry::CGeometryCalculator::CGeometryCalculationInfo& info2)
{
    for (int k = 0; k < 2; ++k)
    {
        scanner[k] = info2.scanner[k];
        plumeCentre[k] = info2.plumeCentre[k];
    }
    return *this;
}

/** Rotates the given vector the given angle [degrees] around the given axis
        @param vec - the coordiates of the vector
        @param angle - the angle to rotate, in degrees
        @param axis - the axis to rotate around (1,2 or 3) */
void CGeometryCalculator::Rotate(double vec[3], double angle, int axis)
{
    double COS = cos(angle * DEGREETORAD);
    double SIN = sin(angle * DEGREETORAD);
    double a = vec[0], b = vec[1], c = vec[2];

    if (axis == 1)
    {
        /** Rotation around X - axis*/
        a = vec[0];
        b = COS * vec[1] + SIN * vec[2];
        c = -SIN * vec[1] + COS * vec[2];
    }
    else if (axis == 2)
    {
        /** Rotation around Y - axis*/
        a = COS * vec[0] - SIN * vec[2];
        b = vec[1];
        c = SIN * vec[0] + COS * vec[2];
    }
    else if (axis == 3)
    {
        /** Rotation around Z - axis*/
        a = COS * vec[0] + SIN * vec[1];
        b = -SIN * vec[0] + COS * vec[1];
        c = vec[2];
    }

    vec[0] = a;
    vec[1] = b;
    vec[2] = c;
}

bool CGeometryCalculator::Intersection(const double o1[3], const double d1[3], const double o2[3], const double d2[3], double& t1, double& t2)
{
    double eps = 1e-19;
    double d1_cross_d2[3], point1[3], point2[3];
    double o2_minus_o1[3];

    // calculate the cross-product (d1 x d2)
    Cross(d1, d2, d1_cross_d2);

    // calculate the squared norm: ||d1 x d2||^2
    double N2 = Norm2(d1_cross_d2);

    if (fabs(N2) < eps)
    {
        /** The lines are parallel */
        t1 = 0;		t2 = 0;
        return false;
    }

    // calculate the distance between the origins
    o2_minus_o1[0] = o2[0] - o1[0];
    o2_minus_o1[1] = o2[1] - o1[1];
    o2_minus_o1[2] = o2[2] - o1[2];

    // Calculate the first determinant
    double det1 = Det(o2_minus_o1, d2, d1_cross_d2);

    // Calculate the second determinant
    double det2 = Det(o2_minus_o1, d1, d1_cross_d2);

    // The result...
    t1 = det1 / N2;
    t2 = det2 / N2;

    // See if the lines do intersect or not
    PointOnRay(o1, d1, t1, point1);
    PointOnRay(o2, d2, t2, point2);

    if (fabs(point1[0] - point2[0]) > eps && fabs(point1[1] - point2[1]) > eps && fabs(point1[2] - point2[2]) > eps)
        return false;

    return true;
}

void CGeometryCalculator::PointOnRay(const double origin[3], const double direction[3], double t, double point[3])
{
    for (int k = 0; k < 3; ++k)
        point[k] = origin[k] + t * direction[k];
}

bool CGeometryCalculator::GetPlumeHeight_Exact(const Configuration::CInstrumentLocation locations[2], const double plumeCentre[2], double& plumeHeight)
{
    CGPSData gps[2] = { CGPSData(locations[0].m_latitude, locations[0].m_longitude, locations[0].m_altitude),
                        CGPSData(locations[1].m_latitude, locations[1].m_longitude, locations[1].m_altitude) };
    double compass[2] = { locations[0].m_compass,		locations[1].m_compass };
    double coneAngle[2] = { locations[0].m_coneangle,	locations[1].m_coneangle };
    double tilt[2] = { locations[0].m_tilt,			locations[1].m_tilt };

    return GetPlumeHeight_Exact(gps, compass, plumeCentre, coneAngle, tilt, plumeHeight);
}

bool CGeometryCalculator::GetPlumeHeight_Exact(const CGPSData gps[2], const double compass[2], const double plumeCentre[2], const double coneAngle[2], const double tilt[2], double& plumeHeight)
{
    double distance, bearing;
    double posLower[3] = { 0, 0, 0 }; // <-- the position of the lower scanner in our changed coordinate system
    double posUpper[3];						// <-- the position of the higher scanner in our changed coordinate system

    // 1. To make the calculations easier, we put a changed coordinate system
    //		on the lowest of the two scanners and calculate the position of the 
    //		other scanner in this coordinate system.
    int lowerScanner = (gps[0].m_altitude < gps[1].m_altitude) ? 0 : 1;
    int upperScanner = 1 - lowerScanner;

    // 2. The distance between the two systems
    distance = GpsMath::Distance(gps[lowerScanner].m_latitude, gps[lowerScanner].m_longitude,
        gps[upperScanner].m_latitude, gps[upperScanner].m_longitude);

    // 3. The bearing from the lower to the higher system (degrees from north, counted clock-wise)
    bearing = GpsMath::Bearing(gps[lowerScanner].m_latitude, gps[lowerScanner].m_longitude,
        gps[upperScanner].m_latitude, gps[upperScanner].m_longitude);

    // 4. The position of the upper scanner 
    posUpper[0] = distance * cos(bearing * DEGREETORAD);
    posUpper[1] = distance * sin(-bearing * DEGREETORAD);
    posUpper[2] = gps[upperScanner].m_altitude - gps[lowerScanner].m_altitude;

    // 5. The directions of the two plume-center rays (defined in the coordinate systems of each scanner)
    double dirLower[3], dirUpper[3]; // <-- the directions
    GetDirection(dirLower, plumeCentre[lowerScanner], coneAngle[lowerScanner], tilt[lowerScanner]);
    GetDirection(dirUpper, plumeCentre[upperScanner], coneAngle[upperScanner], tilt[upperScanner]);

    // 6. Find the direction of the plume-center ray of the upper scanner
    //		in the coordinate system of the lower scanner
    Rotate(dirUpper, compass[upperScanner] - compass[lowerScanner], 3);

    // 7. Find the position of the upper scanner in the coordinate-system
    //		of the lower scanner.
    CGeometryCalculator::Rotate(posUpper, -compass[lowerScanner], 3);

    // 8. Calculate the intersection point of the two rays
    double t1, t2;
    Normalize(dirLower);
    Normalize(dirUpper);
    bool hit = Intersection(posLower, dirLower, posUpper, dirUpper, t1, t2);


    // 9. The plume-height (above the lower scanner) is the z-component of 
    //		the intersection-point
    if (hit)
    {
        double intersectionPoint[3];
        PointOnRay(posLower, dirLower, t1, intersectionPoint); // <-- calculate the intersection point
        plumeHeight = intersectionPoint[2];
    }
    else
    {
        // if the rays don't actually hit each other, calculate the distance between
        //	them. If this is small enough let's consider them as a hit.
        double point1[3], point2[3];
        PointOnRay(posLower, dirLower, t1, point1); // <-- calculate the intersection point
        PointOnRay(posUpper, dirUpper, t2, point2); // <-- calculate the intersection point
        double distance2 = pow(point1[0] - point2[0], 2) + pow(point1[1] - point2[1], 2) + pow(point1[2] - point2[2], 2);
        if (distance2 > 1600)
            return false; // the distance between the intersection points is > 400 m!!

        // take the plume-height as the average of the heights of the two intersection-points
        plumeHeight = (point1[2] + point2[2]) * 0.5;
    }

    return true;
}

bool CGeometryCalculator::GetPlumeHeight_Fuzzy(const CGPSData source, const Configuration::CInstrumentLocation locations[2], const double plumeCentre[2], double& plumeHeight, double& windDirection)
{
    CGPSData gps[2] = { CGPSData(locations[0].m_latitude, locations[0].m_longitude, locations[0].m_altitude),
                        CGPSData(locations[1].m_latitude, locations[1].m_longitude, locations[1].m_altitude) };
    double compass[2] = { locations[0].m_compass,		locations[1].m_compass };
    double coneAngle[2] = { locations[0].m_coneangle,	locations[1].m_coneangle };
    double tilt[2] = { locations[0].m_tilt,			locations[1].m_tilt };

    return GetPlumeHeight_Fuzzy(source, gps, compass, plumeCentre, coneAngle, tilt, plumeHeight, windDirection);
}

bool CGeometryCalculator::GetPlumeHeight_Fuzzy(const CGPSData source, const CGPSData gps[2], const double compass[2], const double plumeCentre[2], const double coneAngle[2], const double tilt[2], double& plumeHeight, double& windDirection)
{
    // 1. To make the calculations easier, we put a changed coordinate system
    //		on the lowest of the two scanners and calculate the position of the 
    //		other scanner in this coordinate system.
    int lowerScanner = (gps[0].m_altitude < gps[1].m_altitude) ? 0 : 1;
    int upperScanner = 1 - lowerScanner;
    double heightDifference = gps[upperScanner].m_altitude - gps[lowerScanner].m_altitude;

    // 2. Find the plume height that gives the same wind-direction for the two instruments
    double guess = 1000;	// the current guess for the plume height
    double h = 10.0;	// the step we use when searching for the plume height
    double maxDiff = 1.0;	// the maximum allowed difference in wind-direction, the convergence criterion

    // 2a. Make an initial guess of the plume height...
    if (gps[lowerScanner].m_altitude > 0 && source.m_altitude > 0)
    {
        guess = std::min(5000.0, std::max(0.0, source.m_altitude - gps[lowerScanner].m_altitude));
    }

    // ------------------------ HERE FOLLOW THE NEW ITERATION ALGORITHM -------------------
    double f = 1e9, f_plus = 1e9;
    double f1, f2;
    int nIterations = 0;
    while (1)
    {
        // Calculate the wind-direction for the current guess of the plume height
        f1 = GetWindDirection(source, guess, gps[lowerScanner], compass[lowerScanner], plumeCentre[lowerScanner], coneAngle[lowerScanner], tilt[lowerScanner]);
        f2 = GetWindDirection(source, guess - heightDifference, gps[upperScanner], compass[upperScanner], plumeCentre[upperScanner], coneAngle[upperScanner], tilt[upperScanner]);
        f = std::max(f1, f2) - std::min(f1, f2);
        if (f > 180.0)
            f = 360.0 - f;

        // Calculate the wind-direction for a plume height a little bit higher than the current guess of the plume height
        f1 = GetWindDirection(source, guess + h, gps[lowerScanner], compass[lowerScanner], plumeCentre[lowerScanner], coneAngle[lowerScanner], tilt[lowerScanner]);
        f2 = GetWindDirection(source, guess + h - heightDifference, gps[upperScanner], compass[upperScanner], plumeCentre[upperScanner], coneAngle[upperScanner], tilt[upperScanner]);
        f_plus = std::max(f1, f2) - std::min(f1, f2);
        if (f_plus > 180.0)
            f_plus = 360.0 - f_plus;

        // Check if we have a good enough result already
        if (f < maxDiff)
        {
            plumeHeight = guess;
            windDirection = (f1 + f2) / 2.0; // take the average wind-direction
            if (plumeHeight < 0 || plumeHeight > 10000)
                return false;
            else
                return true;
        }
        else if (f_plus < maxDiff)
        {
            plumeHeight = guess + h;
            windDirection = (f1 + f2) / 2.0; // take the average wind-direction
            if (plumeHeight < 0 || plumeHeight > 10000)
                return false;
            else
                return true;
        }

        // the local derivative
        double dfdx = (f_plus - f) / h;

        // one step using the Newton method, make a line-search
        //	of the step-size to guarantee that we do decrease
        //	the difference at each step
        double alpha = 0.5;
        double newGuess = guess - alpha * f / dfdx;
        f1 = GetWindDirection(source, newGuess, gps[lowerScanner], compass[lowerScanner], plumeCentre[lowerScanner], coneAngle[lowerScanner], tilt[lowerScanner]);
        f2 = GetWindDirection(source, newGuess - heightDifference, gps[upperScanner], compass[upperScanner], plumeCentre[upperScanner], coneAngle[upperScanner], tilt[upperScanner]);
        double f_new = fabs(f1 - f2);
        if (f_new > 180.0)
        {
            f_new = 360.0 - f_new;
        }

        int nIterations2 = 0;
        while (f_new > f)
        {
            alpha = alpha / 2;
            newGuess = guess - alpha * f / dfdx;
            f1 = GetWindDirection(source, newGuess, gps[lowerScanner], compass[lowerScanner], plumeCentre[lowerScanner], coneAngle[lowerScanner], tilt[lowerScanner]);
            f2 = GetWindDirection(source, newGuess - heightDifference, gps[upperScanner], compass[upperScanner], plumeCentre[upperScanner], coneAngle[upperScanner], tilt[upperScanner]);
            f_new = fabs(f1 - f2);
            if (f_new > 180.0)
            {
                f_new = 360.0 - f_new;
            }

            if (nIterations2++ > 1000)
            {
                return false;
            }
        }
        if (f_new < maxDiff)
        {
            plumeHeight = newGuess;
            windDirection = (f1 + f2) / 2.0; // take the average wind-direction
            return true;
        }
        guess = newGuess;


        // Increase and check the number of iterations
        if (nIterations++ > 100)
        {
            return false;
        }
    }

    // We should never end up here...
    return false;
}

void CGeometryCalculator::GetDirection(double direction[3], double scanAngle, double coneAngle, double tilt)
{
    double tan_coneAngle = tan(coneAngle * DEGREETORAD);
    double cos_tilt = cos(tilt * DEGREETORAD);
    double sin_tilt = sin(tilt * DEGREETORAD);
    double cos_alpha = cos(scanAngle * DEGREETORAD);
    double sin_alpha = sin(scanAngle * DEGREETORAD);
    double divisor = (cos_alpha * cos_tilt + sin_tilt / tan_coneAngle);

    direction[0] = (cos_tilt / tan_coneAngle - cos_alpha * sin_tilt) / divisor;
    direction[1] = sin_alpha / divisor;
    direction[2] = 1;
}

bool CGeometryCalculator::CalculateGeometry(const novac::CString& evalLog1, const novac::CString& evalLog2, const Configuration::CInstrumentLocation locations[2], Geometry::CGeometryResult& result)
{
    return CGeometryCalculator::CalculateGeometry(evalLog1, 0, evalLog2, 0, locations, result);
}

bool CGeometryCalculator::CalculateGeometry(const novac::CString& evalLog1, int scanIndex1, const novac::CString& evalLog2, int scanIndex2, const Configuration::CInstrumentLocation locations[2], Geometry::CGeometryResult& result)
{
    FileHandler::CEvaluationLogFileHandler reader[2];
    CGPSData source;
    CPlumeInScanProperty plume[2];
    CDateTime startTime[2];
    int k; // iterator

    // 1. Read the evaluation-logs
    reader[0].m_evaluationLog.Format("%s", (const char*)evalLog1);
    reader[1].m_evaluationLog.Format("%s", (const char*)evalLog2);
    if (RETURN_CODE::SUCCESS != reader[0].ReadEvaluationLog())
        return false;
    if (RETURN_CODE::SUCCESS != reader[1].ReadEvaluationLog())
        return false;

    // 2. Get the 'CPlumeInScanProperty' for the two scans and the start-times
    int index[2] = { scanIndex1, scanIndex2 };
    for (k = 0; k < 2; ++k)
    {
        reader[k].m_scan[index[k]].GetStartTime(0, startTime[k]);

        if (false == reader[k].m_scan[index[k]].CalculatePlumeCentre(CMolecule(g_userSettings.m_molecule), plume[k]))
        {
            return false; // <-- cannot see the plume
        }
        if (plume[k].completeness < g_userSettings.m_calcGeometry_CompletenessLimit + 0.01)
        {
            return false; // <-- cannot see enough of the plume
        }
    }

    // 3. Calculate the geometry
    return CalculateGeometry(plume[0], startTime[0], plume[1], startTime[1], locations, result);
}

bool CGeometryCalculator::CalculateGeometry(const CPlumeInScanProperty& plume1, const CDateTime& startTime1, const CPlumeInScanProperty& plume2, const CDateTime& startTime2, const Configuration::CInstrumentLocation locations[2], Geometry::CGeometryResult& result)
{
    CGPSData source;
    CDateTime startTime[2];
    double plumeCentre_perturbated[2];
    int k; // iterator

    // 2. Get the nearest volcanoes, if these are different then quit the calculations
    int volcanoIndex1 = g_volcanoes.GetVolcanoIndex(locations[0].m_volcano);
    int volcanoIndex2 = g_volcanoes.GetVolcanoIndex(locations[1].m_volcano);
    if ((volcanoIndex1 == -1) || (volcanoIndex1 != volcanoIndex2))
        return false; // if we couldn't find any volcano or we found two different volcanoes...

    source.m_latitude = g_volcanoes.GetPeakLatitude(volcanoIndex1);
    source.m_longitude = g_volcanoes.GetPeakLongitude(volcanoIndex1);
    source.m_altitude = (long)g_volcanoes.GetPeakAltitude(volcanoIndex1);

    // 4. Get the scan-angles around which the plumes are centred and the start-times of the scans
    double plumeCentre[2] = { plume1.plumeCenter, plume2.plumeCenter2 };

    // 5. Calculate the plume-height
    if (false == CGeometryCalculator::GetPlumeHeight_Fuzzy(source, locations, plumeCentre, result.m_plumeAltitude, result.m_windDirection))
    {
        return false; // <-- could not calculate plume-height
    }
    if (result.m_plumeAltitude < 0)
    {
        return false; // we failed to calculate anything reasonable
    }

    // 7. We also need an estimate of the errors in plume height and wind direction

    // 7a. The error in plume height and wind-direction due to uncertainty in finding the centre of the plume
    double ph_perp[4], wd_perp[4];
    for (k = 0; k < 4; ++k)
    {
        // make a small perturbation to the plume centre angles
        plumeCentre_perturbated[0] = plumeCentre[0] + plume1.plumeCenterError * ((k % 2 == 0) ? -1.0 : +1.0);
        plumeCentre_perturbated[1] = plumeCentre[1] + plume2.plumeCenterError * ((k < 2) ? -1.0 : +1.0);

        // make sure that the perturbation is not too large...
        if ((fabs(plumeCentre_perturbated[0]) > 89.0) || (fabs(plumeCentre_perturbated[1]) > 89.0))
        {
            ph_perp[k] = 1e99;
            continue;
        }

        // try to calculate the  plume height with the perturbed plume centre angles
        if (false == CGeometryCalculator::GetPlumeHeight_Fuzzy(source, locations, plumeCentre_perturbated, ph_perp[k], wd_perp[k]))
        {
            ph_perp[k] = 1e99; // <-- could not calculate plume-height
        }
    }
    result.m_plumeAltitudeError = (fabs(ph_perp[0] - result.m_plumeAltitude) +
        fabs(ph_perp[1] - result.m_plumeAltitude) +
        fabs(ph_perp[2] - result.m_plumeAltitude) +
        fabs(ph_perp[3] - result.m_plumeAltitude)) / 4;
    result.m_windDirectionError = (fabs(wd_perp[0] - result.m_windDirection) +
        fabs(wd_perp[1] - result.m_windDirection) +
        fabs(wd_perp[2] - result.m_windDirection) +
        fabs(wd_perp[3] - result.m_windDirection)) / 4;

    // 7b. Also scale the altitude error with the time difference between the two scans
    double timeDifference_Minutes = fabs(CDateTime::Difference(startTime[0], startTime[1])) / 60.0;
    result.m_plumeAltitudeError *= pow(2.0, timeDifference_Minutes / 30.0);

    // 7c. Remember to add the altitude of the lowest scanner to the plume height to get the total plume altitude
    result.m_plumeAltitude += std::min(locations[0].m_altitude, locations[1].m_altitude);
    // double plumeAltitudeRelativeToScanner0	= result.m_plumeAltitude - locations[0].m_altitude;

    // 8. Also store the date the measurements were made and the average-time
    double timeDifference = CDateTime::Difference(startTime1, startTime2);
    if (timeDifference < 0)
    {
        result.m_averageStartTime = startTime1;
        result.m_averageStartTime.Increment((int)fabs(timeDifference) / 2);
    }
    else
    {
        result.m_averageStartTime = startTime2;
        result.m_averageStartTime.Increment((int)(timeDifference / 2));
    }
    result.m_startTimeDifference = (int)fabs(timeDifference);

    // 9. The parameters about the scans that were combined
    result.m_plumeCentre1 = (float)plume1.plumeCenter;
    result.m_plumeCentreError1 = (float)plume1.plumeCenterError;
    result.m_plumeCentre2 = (float)plume2.plumeCenter2;       // changed 2019-02-20: Was plumeCenter
    result.m_plumeCentreError2 = (float)plume2.plumeCenterError2;  // changed 2019-02-20: Was plumeCenterError

    return true;
}

/** Calculates the cross product of the supplied vectors */
void CGeometryCalculator::Cross(const double u[3], const double v[3], double result[3])
{
    result[0] = u[1] * v[2] - u[2] * v[1];
    result[1] = u[2] * v[0] - u[0] * v[2];
    result[2] = u[0] * v[1] - u[1] * v[0];
}

/** Calculates the squared norm of the supplied vector */
double CGeometryCalculator::Norm2(const double v[3])
{
    return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/** Calculates the determinant of a matrix whose columns are defined
        by the three supplied vectors */
double CGeometryCalculator::Det(const double c1[3], const double c2[3], const double c3[3])
{
    double ret = c1[0] * c2[1] * c3[2] + c2[0] * c3[1] * c1[2] + c3[0] * c1[1] * c2[2];
    ret = ret - c1[0] * c3[1] * c2[2] - c2[0] * c1[1] * c3[2] - c3[0] * c2[1] * c1[2];

    return ret;
}

/** Normalizes the supplied vector */
void CGeometryCalculator::Normalize(double v[3])
{
    double norm_inv = 1 / sqrt(Norm2(v));
    v[0] *= norm_inv;
    v[1] *= norm_inv;
    v[2] *= norm_inv;
}

double CGeometryCalculator::GetWindDirection(const CGPSData source, double plumeHeight, const Configuration::CInstrumentLocation scannerLocation, double plumeCentre)
{
    CGPSData gps = CGPSData(scannerLocation.m_latitude, scannerLocation.m_longitude, scannerLocation.m_altitude);

    return GetWindDirection(source, plumeHeight, gps, scannerLocation.m_compass, plumeCentre, scannerLocation.m_coneangle, scannerLocation.m_tilt);
}

double CGeometryCalculator::GetWindDirection(const CGPSData source, double plumeHeight, const CGPSData scannerPos, double compass, double plumeCentre, double coneAngle, double tilt)
{
    if (plumeCentre == NOT_A_NUMBER)
        return NOT_A_NUMBER;

    // 1. Calculate the intersection-point
    double intersectionDistance, angle;
    if (fabs(coneAngle - 90.0) > 1)
    {
        // ------------ CONE SCANNERS -----------
        // 1a. the distance from the system to the intersection-point
        double x, y;
        double cos_tilt = cos(DEGREETORAD * tilt);
        double sin_tilt = sin(DEGREETORAD * tilt);
        double tan_coneAngle = tan(DEGREETORAD * coneAngle);
        double cos_alpha = cos(DEGREETORAD * plumeCentre);
        double sin_alpha = sin(DEGREETORAD * plumeCentre);

        // Calculate the projections of the intersection points in the ground-plane
        double commonDenominator = cos_alpha * cos_tilt + sin_tilt / tan_coneAngle;
        x = (cos_tilt / tan_coneAngle - cos_alpha * sin_tilt) / commonDenominator;
        y = (sin_alpha) / commonDenominator;

        intersectionDistance = plumeHeight * sqrt(pow(x, 2) + pow(y, 2));

        // 1b. the direction from the system to the intersection-point
        angle = atan2(y, x) / DEGREETORAD + compass;
    }
    else
    {
        // ------------- FLAT SCANNERS ---------------
        // 1a. the distance from the system to the intersection-point
        intersectionDistance = plumeHeight * tan(DEGREETORAD * plumeCentre);

        // 1b. the direction from the system to the intersection-point
        if (plumeCentre == 0)
            angle = 0;
        else if (plumeCentre < 0)
            angle = (compass + 90);
        else
            angle = (compass - 90);
    }

    // 1c. the intersection-point
    double lat2, lon2;
    GpsMath::CalculateDestination(scannerPos.m_latitude, scannerPos.m_longitude, intersectionDistance, angle, lat2, lon2);

    // 2. the wind-direction
    double windDirection = GpsMath::Bearing(lat2, lon2, source.m_latitude, source.m_longitude);

    return windDirection;
}

bool CGeometryCalculator::CalculatePlumeHeight(const novac::CString& evalLog, int scanIndex, Meteorology::CWindField& windField, Configuration::CInstrumentLocation location, Geometry::CGeometryResult& result)
{
    FileHandler::CEvaluationLogFileHandler reader;
    CPlumeInScanProperty plume;
    CGPSData source, scannerPos;

    // extract the location of the instrument
    scannerPos = CGPSData(location.m_latitude, location.m_longitude, location.m_altitude);

    // Extract the location of the source
    int volcanoIndex1 = g_volcanoes.GetVolcanoIndex(location.m_volcano);
    if (volcanoIndex1 == -1)
        return false; // if we couldn't find any volcano 
    source.m_latitude = g_volcanoes.GetPeakLatitude(volcanoIndex1);
    source.m_longitude = g_volcanoes.GetPeakLongitude(volcanoIndex1);
    source.m_altitude = (long)g_volcanoes.GetPeakAltitude(volcanoIndex1);

    // 3. Read the evaluation-log
    reader.m_evaluationLog.Format("%s", (const char*)evalLog);
    if (RETURN_CODE::SUCCESS != reader.ReadEvaluationLog())
        return false;

    // 4. Get the scan-angles around which the plumes are centred
    if (false == reader.m_scan[scanIndex].CalculatePlumeCentre(CMolecule(g_userSettings.m_molecule), plume))
    {
        return false; // <-- cannot see the plume
    }
    if (plume.completeness < g_userSettings.m_calcGeometry_CompletenessLimit + 0.01)
    {
        return false; // <-- cannot see enough of the plume
    }

    // calculate the plume height
    double plumeHeight = CGeometryCalculator::GetPlumeHeight(source, windField.GetWindDirection(), scannerPos, location.m_compass, plume.plumeCenter, location.m_coneangle, location.m_tilt);

    // Check that the plume height is reasonable
    if (plumeHeight < 0)
        return false;

    // Also try to estimate the error in the plume height measurement
    double windDirection_plus = windField.GetWindDirection() + std::max(windField.GetWindDirectionError(), 5.0);
    double windDirection_minus = windField.GetWindDirection() - std::max(windField.GetWindDirectionError(), 5.0);

    // the error in plume height due to the uncertainty in wind direction
    double plumeHeight_plus_wd = CGeometryCalculator::GetPlumeHeight(source, windDirection_plus, scannerPos, location.m_compass, plume.plumeCenter, location.m_coneangle, location.m_tilt);
    double plumeHeight_minus_wd = CGeometryCalculator::GetPlumeHeight(source, windDirection_minus, scannerPos, location.m_compass, plume.plumeCenter, location.m_coneangle, location.m_tilt);

    // the error in plume height due to the uncertainty in the plume centre position
    double plumeHeight_plus_pc = CGeometryCalculator::GetPlumeHeight(source, windField.GetWindDirection(), scannerPos, location.m_compass, plume.plumeCenter + plume.plumeCenterError, location.m_coneangle, location.m_tilt);
    double plumeHeight_minus_pc = CGeometryCalculator::GetPlumeHeight(source, windField.GetWindDirection(), scannerPos, location.m_compass, plume.plumeCenter - plume.plumeCenterError, location.m_coneangle, location.m_tilt);

    // the total error in plume height
    double plumeHeightErr = sqrt(pow(plumeHeight_plus_wd - plumeHeight_minus_wd, 2.0) + pow(plumeHeight_plus_pc - plumeHeight_minus_pc, 2.0));

#ifdef _DEBUG
    novac::CString fileName;
    fileName.Format("%s%cdebugGeometrySingleInstr.txt", (const char*)g_userSettings.m_outputDirectory, Poco::Path::separator());
    FILE* f = fopen(fileName, "a");
    if (f > 0)
    {
        fprintf(f, "%.1lf\t%.2lf\t%.2lf\n", plume.plumeCenter, plumeHeight, plumeHeightErr);
        fclose(f);
    }
#endif

    if (plumeHeightErr > g_userSettings.m_calcGeometry_MaxPlumeAltError)
        return false;

    reader.m_scan[scanIndex].GetStartTime(0, result.m_averageStartTime);
    result.m_plumeAltitude = plumeHeight + location.m_altitude;
    result.m_plumeAltitudeError = plumeHeightErr;
    result.m_windDirection = NOT_A_NUMBER;
    result.m_windDirectionError = 0.0;
    result.m_plumeCentre1 = (float)plume.plumeCenter;
    result.m_plumeCentreError1 = (float)plume.plumeCenterError;
    result.m_calculationType = Meteorology::MET_GEOMETRY_CALCULATION_SINGLE_INSTR;

    return true;
}

bool CGeometryCalculator::CalculateWindDirection(const novac::CString& evalLog, int scanIndex, Geometry::CPlumeHeight& absolutePlumeHeight, Configuration::CInstrumentLocation location, Geometry::CGeometryResult& result)
{
    FileHandler::CEvaluationLogFileHandler reader;
    CPlumeInScanProperty plume;
    CGPSData source, scannerPos;

    // extract the location of the instrument
    scannerPos = CGPSData(location.m_latitude, location.m_longitude, location.m_altitude);

    // Extract the location of the source
    int volcanoIndex1 = g_volcanoes.GetVolcanoIndex(location.m_volcano);
    if (volcanoIndex1 == -1)
        return false; // if we couldn't find any volcano 
    source.m_latitude = g_volcanoes.GetPeakLatitude(volcanoIndex1);
    source.m_longitude = g_volcanoes.GetPeakLongitude(volcanoIndex1);
    source.m_altitude = (long)g_volcanoes.GetPeakAltitude(volcanoIndex1);

    // 3. Read the evaluation-log
    reader.m_evaluationLog.Format("%s", (const char*)evalLog);
    if (RETURN_CODE::SUCCESS != reader.ReadEvaluationLog())
        return false;

    // 4. Get the scan-angles around which the plumes are centred
    if (false == reader.m_scan[scanIndex].CalculatePlumeCentre(CMolecule(g_userSettings.m_molecule), plume))
    {
        return false; // <-- cannot see the plume
    }
    if (plume.completeness < g_userSettings.m_calcGeometry_CompletenessLimit + 0.01)
    {
        return false; // <-- cannot see enough of the plume
    }

    // the relative plume height
    double plumeHeight = absolutePlumeHeight.m_plumeAltitude - scannerPos.m_altitude;
    if (plumeHeight <= 0.0)
        return false; // failure!

    // calculate the wind direction
    double windDirection = CGeometryCalculator::GetWindDirection(source, plumeHeight, scannerPos, location.m_compass, plume.plumeCenter, location.m_coneangle, location.m_tilt);

    // Check that the wind direction is reasonable
    if (windDirection <= NOT_A_NUMBER)
        return false;

    // the error in wind direction due to the uncertainty in plume height
    double windDirection_plus_ph = CGeometryCalculator::GetWindDirection(source, plumeHeight + absolutePlumeHeight.m_plumeAltitudeError, scannerPos, location.m_compass, plume.plumeCenter, location.m_coneangle, location.m_tilt);
    double windDirection_minus_ph = CGeometryCalculator::GetWindDirection(source, plumeHeight - absolutePlumeHeight.m_plumeAltitudeError, scannerPos, location.m_compass, plume.plumeCenter, location.m_coneangle, location.m_tilt);

    // the error in wind direction due to the uncertainty in the plume centre position
    double windDirection_plus_pc = CGeometryCalculator::GetWindDirection(source, plumeHeight, scannerPos, location.m_compass, plume.plumeCenter + plume.plumeCenterError, location.m_coneangle, location.m_tilt);
    double windDirection_minus_pc = CGeometryCalculator::GetWindDirection(source, plumeHeight, scannerPos, location.m_compass, plume.plumeCenter - plume.plumeCenterError, location.m_coneangle, location.m_tilt);

    // the total error in wind direction
    double windDirectionErr = sqrt(pow(windDirection_plus_ph - windDirection_minus_ph, 2.0) + pow(windDirection_plus_pc - windDirection_minus_pc, 2.0));


#ifdef _DEBUG
    novac::CString fileName;
    fileName.Format("%s%cdebugGeometrySingleInstr.txt", (const char*)g_userSettings.m_outputDirectory, Poco::Path::separator());
    FILE* f = fopen(fileName, "a");
    if (f > 0)
    {
        fprintf(f, "wd\t%.1lf\t%.2lf\t%.2lf\n", plume.plumeCenter, windDirection, windDirectionErr);
        fclose(f);
    }
#endif

    if (windDirectionErr > g_userSettings.m_calcGeometry_MaxWindDirectionError)
        return false;

    reader.m_scan[scanIndex].GetStartTime(0, result.m_averageStartTime);
    result.m_plumeAltitude = NOT_A_NUMBER;
    result.m_plumeAltitudeError = 0.0;
    result.m_windDirection = windDirection;
    result.m_windDirectionError = windDirectionErr;
    result.m_plumeCentre1 = (float)plume.plumeCenter;
    result.m_plumeCentreError1 = (float)plume.plumeCenterError;
    result.m_calculationType = Meteorology::MET_GEOMETRY_CALCULATION_SINGLE_INSTR;

    return true;
}

double CGeometryCalculator::GetPlumeHeight(const CGPSData source, double windDirection, const CGPSData scannerPos, double compass, double plumeCentre, double coneAngle, double tilt)
{
    if (plumeCentre == NOT_A_NUMBER)
        return NOT_A_NUMBER;

    // 1. prepare by calculating the sine and cosine fo the wind direction
    double sin_wd = sin(DEGREETORAD * (windDirection - compass));
    double cos_wd = cos(DEGREETORAD * (windDirection - compass));

    // 2. Calculate the location of the source in the coordinate system that has its
    //		origin in the scanner
    double distanceToSource = GpsMath::Distance(source.m_latitude, source.m_longitude, scannerPos.m_latitude, scannerPos.m_longitude);
    double directionToSource = GpsMath::Bearing(scannerPos.m_latitude, scannerPos.m_longitude, source.m_latitude, source.m_longitude);
    double xs = distanceToSource * cos(DEGREETORAD * (compass - directionToSource));
    double ys = distanceToSource * sin(DEGREETORAD * (compass - directionToSource));

    // 1. Formulate the line emerging from the scanner
    double dx, dy;
    if (fabs(coneAngle - 90.0) > 1)
    {
        // ------------ CONE SCANNERS -----------
        double cos_tilt = cos(DEGREETORAD * tilt);
        double sin_tilt = sin(DEGREETORAD * tilt);
        double tan_coneAngle = tan(DEGREETORAD * coneAngle);
        double cos_alpha = cos(DEGREETORAD * plumeCentre);
        double sin_alpha = sin(DEGREETORAD * plumeCentre);

        double commonDenominator = cos_alpha * cos_tilt + sin_tilt / tan_coneAngle;
        dx = (cos_tilt / tan_coneAngle - cos_alpha * sin_tilt) / commonDenominator;
        dy = (sin_alpha) / commonDenominator;
    }
    else
    {
        // ------------- FLAT SCANNERS ---------------
        dx = 0;
        dy = tan(DEGREETORAD * plumeCentre);
    }

    // 2. Calculate the intersection point between the line emerging from the scanner
    //		and the plane which contains the source and the wind-direction vector
    double denominator = dx * sin_wd - dy * cos_wd;

    if (fabs(denominator) < 0.001)
    {
        // the line does not intersect the plane
        return NOT_A_NUMBER;
    }

    double plumeHeight = (xs * sin_wd - ys * cos_wd) / denominator;

    return plumeHeight;
}

double CGeometryCalculator::GetWindDirection(const CGPSData source, const CGPSData scannerPos, double plumeHeight, double alpha_center_of_mass, double phi_center_of_mass)
{
    //longitudinal distance between instrument and source:
    double x_source = GpsMath::Distance(scannerPos.m_latitude, source.m_longitude, scannerPos.m_latitude, scannerPos.m_longitude);
    if (source.m_longitude < scannerPos.m_longitude)
        x_source = -fabs(x_source);
    else
        x_source = fabs(x_source);

    //latitudinal distance between instrument and source:
    double y_source = GpsMath::Distance(source.m_latitude, scannerPos.m_longitude, scannerPos.m_latitude, scannerPos.m_longitude);
    if (source.m_latitude < scannerPos.m_latitude)
        y_source = -fabs(y_source);
    else
        y_source = fabs(y_source);

    //the two angles for the measured center of mass of the plume converted to rad:
    double alpha_cm_rad = DEGREETORAD * alpha_center_of_mass;
    double phi_cm_rad = DEGREETORAD * phi_center_of_mass;

    double wd = atan2((x_source - plumeHeight * tan(alpha_cm_rad) * sin(phi_cm_rad)), (y_source - plumeHeight * tan(alpha_cm_rad) * cos(phi_cm_rad))) / DEGREETORAD;
    if (wd < 0)
        wd += 360;		//because atan2 returns values between -pi...+pi

    return wd;
}

double CGeometryCalculator::GetPlumeHeight_OneInstrument(const CGPSData source, const CGPSData gps, double WindDirection, double alpha_center_of_mass, double phi_center_of_mass)
{
    //horizontal distance between instrument and source:
    double distance_to_source = GpsMath::Distance(gps.m_latitude, gps.m_longitude, source.m_latitude, source.m_longitude);

    //angle (in rad) pointing from instrument to source (with respect to north, clockwise):
    double angle_to_source_rad = DEGREETORAD * GpsMath::Bearing(gps.m_latitude, gps.m_longitude, source.m_latitude, source.m_longitude);

    //the two angles for the measured center of mass of the plume converted to rad:
    double alpha_cm_rad = DEGREETORAD * alpha_center_of_mass;
    double phi_cm_rad = DEGREETORAD * phi_center_of_mass;

    double WindDirection_rad = DEGREETORAD * WindDirection;

    return 1 / tan(alpha_cm_rad) * sin(angle_to_source_rad - WindDirection_rad) / sin(phi_cm_rad - WindDirection_rad) * distance_to_source;

}
