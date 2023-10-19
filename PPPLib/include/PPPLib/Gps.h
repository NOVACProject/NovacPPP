#pragma once

namespace Gps
{
    class GpsMath
    {
        public:

            /* This function returns the distance in <b>meters</b> between the two points defined
            by (lat1,lon1) and (lat2, lon2). <b>All angles must be in degrees</b> */
            static double GPSDistance(double lat1, double lon1, double lat2, double lon2);
    };
}