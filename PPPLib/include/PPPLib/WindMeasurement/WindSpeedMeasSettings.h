#pragma once
#include <SpectralEvaluation/Definitions.h>

namespace WindSpeedMeasurement
{

/** The <b>CWindSpeedMeasSettings</b> contains the parameters
        necessary for calculating the correlation between data series
        in order to calculate wind speeds. 	*/
class CWindSpeedMeasSettings
{
public:
    /** The number of pixels to average in the low pass filtering.
            If lowPassFilterAverage is 0, then no filtering will be done */
    unsigned int lowPassFilterAverage = 20;

    /** The maximum number of seconds to shift */
    unsigned int shiftMax = 90;

    /** The length of the test-region, in seconds */
    unsigned testLength = 120;

    /** The minimum column value that will be taken into account (normally not used) */
    double columnMin = -1e99;

    /** The minimum sigma - level (???) */
    double sigmaMin = 0.1;

    /** The plume height */
    double plumeHeight = 1000.0;

    /** The angle separation in the instrument. In degrees */
    double angleSeparation = 0.08 / DEGREETORAD;
};
}
