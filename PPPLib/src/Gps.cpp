#include <PPPLib/Gps.h>
#include <SpectralEvaluation/Definitions.h>
#include <cmath>
#include <algorithm>

namespace Gps
{

double GpsMath::GPSDistance(double lat1, double lon1, double lat2, double lon2)
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

}