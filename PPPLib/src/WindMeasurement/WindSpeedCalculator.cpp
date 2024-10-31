#include <PPPLib/WindMeasurement/WindSpeedCalculator.h>
#include <PPPLib/File/EvaluationLogFileHandler.h>
#include <PPPLib/Meteorology/WindField.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/Logging.h>
#include <cstring>
#include <cmath>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

using namespace WindSpeedMeasurement;
using namespace novac;


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


CWindSpeedCalculator::CMeasurementSeries::CMeasurementSeries()
{
    column = NULL;
    time = NULL;
    length = 0;
}

CWindSpeedCalculator::CMeasurementSeries::CMeasurementSeries(size_t len)
{
    column = new double[len];
    time = new double[len];
    length = len;
}

RETURN_CODE CWindSpeedCalculator::CMeasurementSeries::SetLength(size_t len)
{
    if (column != NULL)
    {
        delete[] column;
        delete[] time;
    }
    column = new double[len];
    if (column == NULL)
        return RETURN_CODE::FAIL;
    time = new double[len];
    if (time == NULL)
        return RETURN_CODE::FAIL;
    length = len;

    return RETURN_CODE::SUCCESS;
}

double CWindSpeedCalculator::CMeasurementSeries::AverageColumn(int from, int to) const
{
    double sum = 0.0;

    // check input
    if (from > to || from < 0 || to < 0 || from >= length || to >= length)
        return 0.0;

    for (int k = from; k <= to; ++k)
    {
        sum += column[k];
    }
    return sum / (to - from + 1);
}

// calculates and returns the average time between two measurements
double CWindSpeedCalculator::CMeasurementSeries::SampleInterval()
{
    if (length <= 0)
        return 0.0;

    double totalTime = time[length - 1] - time[0];
    return totalTime / (length - 1);
}

CWindSpeedCalculator::CMeasurementSeries::~CMeasurementSeries()
{
    if (length != 0)
    {
        delete[] column;
        delete[] time;
    }
}


CWindSpeedCalculator::CWindSpeedCalculator(novac::ILogger& log, const Configuration::CUserConfiguration& userSettings)
    : m_log(log), m_userSettings(userSettings)
{
    shift = NULL;
    corr = NULL;
    used = NULL;
    delays = NULL;
    m_length = 0;
    m_firstDataPoint = -1;
}

CWindSpeedCalculator::~CWindSpeedCalculator(void)
{
    delete[] shift;
    delete[] corr;
    delete[] used;
}

RETURN_CODE CWindSpeedCalculator::CalculateDelay(
    double& /*delay*/,
    const CMeasurementSeries* upWindSerie,
    const CMeasurementSeries* downWindSerie,
    const CWindSpeedMeasSettings& settings)
{

    int		comparisonLength;

    CMeasurementSeries modifiedUpWind;
    CMeasurementSeries modifiedDownWind;

    // 0. Error checking of the input
    if (upWindSerie == NULL || upWindSerie->length == 0)
        return RETURN_CODE::FAIL;
    if (downWindSerie == NULL || downWindSerie->length == 0)
        return RETURN_CODE::FAIL;

    // 1. Start with the low pass filtering
    if (RETURN_CODE::SUCCESS != LowPassFilter(upWindSerie, &modifiedUpWind, settings.lowPassFilterAverage))
        return RETURN_CODE::FAIL;
    if (RETURN_CODE::SUCCESS != LowPassFilter(downWindSerie, &modifiedDownWind, settings.lowPassFilterAverage))
        return RETURN_CODE::FAIL;

    // 1b. Get the sample time
    double sampleInterval = modifiedDownWind.SampleInterval();
    if (fabs(modifiedUpWind.SampleInterval() - sampleInterval) > 0.5)
    {
        return RETURN_CODE::FAIL; // <-- we cannot have different sample intervals of the two time series
    }

    // 1c. Calculate the how many pixels that we should shift maximum
    int maximumShift = (int)round(settings.shiftMax / sampleInterval);

    // 1d. Calculate the length of the comparison-interval
    //		in data-points instead of in seconds
    if (m_userSettings.m_fUseMaxTestLength_DualBeam)
    {
        comparisonLength = static_cast<int>(modifiedDownWind.length) - 2 * maximumShift - 10;
    }
    else
    {
        comparisonLength = (int)round(settings.testLength / sampleInterval);
    }

    // 1e. check that the resulting series is long enough
    if (modifiedDownWind.length - maximumShift - comparisonLength < maximumShift + 1)
        return RETURN_CODE::FAIL; // <-- data series to short to use current settings of test length and shiftmax

    // 2. Allocate some assistance arrays 
    //		(Note that it is the down wind data series that is shifted)
    m_length = static_cast<int>(modifiedDownWind.length);
    InitializeArrays();

    // The number of datapoints skipped because we cannot see the plume.
    int skipped = 0;

    m_arrayLength = m_length - maximumShift - comparisonLength - 1;

    // 3. Iterate over the set of sub-arrays in the down-wind data series
    //		Offset is the starting-point in this sub-array whos length is 'comparisonLength'
    for (int offset = 0; offset < m_length - maximumShift - comparisonLength; ++offset)
    {
        double highestCorr;
        int bestShift;

        // 3a. Pick out the sub-vectors
        double* series1 = modifiedUpWind.column + offset;
        double* series2 = modifiedDownWind.column + offset;
        size_t series1Length = modifiedUpWind.length - offset;
        size_t series2Length = static_cast<size_t>(comparisonLength);

        // 3b. The midpoint in the subvector
        int midPoint = (int)round(offset + comparisonLength / 2);

        //	3c. Check if we see the plume...
        if (upWindSerie->AverageColumn(offset, offset + comparisonLength) < settings.columnMin)
        {
            // we don't see the plume. 
            skipped = skipped + 1;
            continue;
        }

        // 3d. Do a shifting...
        FindBestCorrelation(series1, series1Length, series2, series2Length, maximumShift, highestCorr, bestShift);

        // 3e. Calculate the time-shift
        delays[midPoint] = bestShift * sampleInterval;
        corr[midPoint] = highestCorr;
        shift[midPoint] = bestShift - 1;
        used[midPoint] = 1;

        if (m_firstDataPoint == -1)
            m_firstDataPoint = midPoint;

    }

    return RETURN_CODE::SUCCESS;
}

/** Performs a low pass filtering on the supplied measurement series.
    The number of iterations in the filtering is given by 'nIterations'
    if nIterations is zero, nothing will be done. */
RETURN_CODE CWindSpeedCalculator::LowPassFilter(const CMeasurementSeries* series, CMeasurementSeries* result, unsigned int nIterations)
{

    // 1. Check for errors in input
    if (series == NULL || series->length == 0 || result == NULL)
        return RETURN_CODE::FAIL;

    // 2. Calculate the old and the new data series lengths
    int length = static_cast<int>(series->length);							// <-- the length of the old data series
    int newLength = length - nIterations - 1;		// <-- the length of the new data series

    if (newLength <= 0)
        return RETURN_CODE::FAIL;

    // 3. If no iterations should be done, the output should be a copy of the input...
    if (nIterations == 0)
    {
        if (RETURN_CODE::SUCCESS != result->SetLength(length))
            return RETURN_CODE::FAIL;

        for (int i = 0; i < length; ++i)
        {
            result->column[i] = series->column[i];
            result->time[i] = series->time[i];
        }
        result->length = series->length;
        return RETURN_CODE::SUCCESS;
    }

    // 4. Change the length of the resulting series.
    if (RETURN_CODE::SUCCESS != result->SetLength(newLength))
        return RETURN_CODE::FAIL;

    // 5. Calculate the factorials
    double* factorial = new double[nIterations];
    factorial[0] = 1;
    if (nIterations > 1)
        factorial[1] = 1;
    for (unsigned int k = 2; k < nIterations; ++k)
        factorial[k] = factorial[k - 1] * (double)k;

    double coeffSum = 0; // <-- the sum of all coefficients, to make sure that we dont increase the values...

    // 6. Allocate temporary arrays for the time-stamps and the column values
    double* tmpCol = new double[newLength];
    double* tmpTime = new double[newLength];
    memset(tmpCol, 0, newLength * sizeof(double));
    memset(tmpTime, 0, newLength * sizeof(double));

    // 7. Do the filtering...
    for (int k = 1; k < (int)nIterations + 1; ++k)
    {
        // 7a. The coefficient in the binomial - expansion
        double coefficient = factorial[nIterations - 1] / (factorial[nIterations - k] * factorial[k - 1]);
        coeffSum += coefficient;

        // 7b. Do the filtering for all data points in the new time series
        for (int i = 0; i < newLength; ++i)
        {
            tmpCol[i] += coefficient * series->column[k - 1 + i];
            tmpTime[i] += coefficient * series->time[k - 1 + i];
        }
    }

    // 8. Divide by the coeffsum to preserve the energy...
    for (int i = 0; i < newLength; ++i)
    {
        result->time[i] = tmpTime[i] / coeffSum;
        result->column[i] = tmpCol[i] / coeffSum;
    }
    result->length = newLength;

    delete[] factorial;
    delete[] tmpCol;
    delete[] tmpTime;

    return RETURN_CODE::SUCCESS;
}

/** Shifts the vector 'shortVector' against the vector 'longVector' and returns the
            shift for which the correlation between the two is highest. */
RETURN_CODE CWindSpeedCalculator::FindBestCorrelation(
    const double* longVector, size_t longLength,
    const double* shortVector, size_t shortLength,
    unsigned int maximumShift,
    double& highestCorr, int& bestShift)
{

    // 0. Check for errors in the input
    if (longLength == 0 || shortLength == 0)
        return RETURN_CODE::FAIL;
    if (longVector == NULL || shortVector == NULL)
        return RETURN_CODE::FAIL;
    if (longLength <= shortLength)
        return RETURN_CODE::FAIL;

    // Reset
    highestCorr = 0;
    bestShift = 0;

    // To calculate the correlation, we need to pick out a subvector (with length 'shortLength)
    //	from the longVector and compare this with 'shortVector' and calculate the correlation.
    // left is the startingpoint of this subvector
    int left = 0;

    // 1. Start shifting
    while ((int)(left + shortLength) < (int)longLength && left < (int)maximumShift)
    {
        double C = correlation(shortVector, longVector + left, shortLength);
        if (C > highestCorr)
        {
            highestCorr = C;
            bestShift = left;
        }
        ++left;
    }

    return RETURN_CODE::SUCCESS;
}

/** Calculates the correlation between the two vectors 'x' and 'y', both of length 'length'
        @return - the correlation between the two vectors. */
double CWindSpeedCalculator::correlation(const double* x, const double* y, size_t length)
{
    double s_xy = 0; // <-- the dot-product X*Y
    double s_x2 = 0; // <-- the dot-product X*X
    double s_x = 0; // <-- sum of all elements in X
    double s_y = 0; // <-- sum of all elements in Y
    double s_y2 = 0; // <-- the dot-product Y*Y
    double c = 0; // <-- the final correlation
    double eps = 1e-5;

    if (length == 0)
        return 0;

    for (size_t k = 0; k < length; ++k)
    {
        s_xy += x[k] * y[k];
        s_x2 += x[k] * x[k];
        s_x += x[k];
        s_y += y[k];
        s_y2 += y[k] * y[k];
    }

    double nom = (length * s_xy - s_x * s_y);
    double denom = sqrt(((length * s_x2 - s_x * s_x) * (length * s_y2 - s_y * s_y)));

    if ((fabs(nom - denom) < eps) && (fabs(denom) < eps))
        c = 1.0;
    else
        c = nom / denom;

    return c;
}

void CWindSpeedCalculator::InitializeArrays()
{
    delete[] shift;
    delete[] corr;
    delete[] used;
    delete[] delays;
    shift = new double[m_length];
    corr = new double[m_length];
    used = new double[m_length];
    delays = new double[m_length]; // <-- the delays
    memset(corr, 0, m_length * sizeof(double));
    memset(shift, 0, m_length * sizeof(double));
    memset(used, 0, m_length * sizeof(double));
    memset(delays, 0, m_length * sizeof(double));
}

/** Calculates the wind speed from the two time series found in the given evaluation-log files.
    If the instrument is a Heidelberg instrument then the second parameter will be ignored.
    @param evalLog1 - the full path and file name of the evaluation-log file
        containing the first of the two time series to correlate
    @param evalLog2 - the full path and file name of the evaluation-log file
        containing the first of the two time series to correlate.
        This will be ignored if the instrument is a Heidelberg type.
    @param plumeHeight - the altitude of the plume at the time of the wind speed
        measurement. This is the full altitude of the plume, in meters above sea level.
    @windField - will on successful return be filled with the derived wind-speed.
        NOTE THAT THE WIND-DIRECTION WILL NOT BE CALCULATED AND WILL THUS NOT MAKE
        ANY SENSE...
*/
int CWindSpeedCalculator::CalculateWindSpeed(const novac::CString& evalLog1, const novac::CString& evalLog2, const Configuration::CInstrumentLocation& location, const Geometry::PlumeHeight& plumeHeight, Meteorology::WindField& windField)
{
    double distance = 0; // the distance between the two viewing directions at the altitude of the plume.

    // Extract the relative plume height
    Geometry::PlumeHeight relativePlumeHeight = plumeHeight;
    relativePlumeHeight.m_plumeAltitude -= location.m_altitude;

    // Calculate the correlation between the time series
    if (location.m_instrumentType == NovacInstrumentType::Gothenburg)
    {
        if (RETURN_CODE::SUCCESS != CalculateCorrelation(evalLog1, evalLog2))
            return 1;
    }
    else if (location.m_instrumentType == NovacInstrumentType::Heidelberg)
    {
        if (RETURN_CODE::SUCCESS != CalculateCorrelation_Heidelberg(evalLog1))
            return 1;
    }

    // we've successfully calculated a correlation between the two time series... 
    //	Now calculate the speed of the plume...

    // 1. Get the distance between the two viewing directions...
    double scanAngle = 0.0;
    if (location.m_instrumentType == NovacInstrumentType::Gothenburg)
    {
        if (fabs(location.m_coneangle - 90.0) < 1)
        {
            // Flat scanner
            distance = relativePlumeHeight.m_plumeAltitude * (1.0 / cos(DEGREETORAD * scanAngle)) * tan(DEGREETORAD * m_settings.angleSeparation);
        }
        else
        {
            // Cone scanner
            double angle = DEGREETORAD * (90.0 - (location.m_coneangle - fabs(location.m_tilt)));
            distance = relativePlumeHeight.m_plumeAltitude * fabs(tan(angle) - tan(angle - DEGREETORAD * m_settings.angleSeparation));
        }
        // If the scanners are not looking straight up then the distance between
        //		the two directions gets decreased with the scanAngle
        distance *= cos(scanAngle * DEGREETORAD);
    }
    else
    {
        //for Heidelberg instrument: recover the distance along the plume
        //	from the plume height and the two measurement directions used
        //	(alpha=angle from zenith; beta=azimuth angle)
        //int halfLength= scan.GetEvaluatedNum() / 2;
        //double alpha1 = DEGREETORAD * scan.GetScanAngle(halfLength);
        //double alpha2 = DEGREETORAD * scan.GetScanAngle(halfLength+1);
        //double beta1  = DEGREETORAD * scan.GetScanAngle2(halfLength);
        //double beta2  = DEGREETORAD * scan.GetScanAngle2(halfLength+1);
        //distance		= m_settings.plumeHeight *
        //							sqrt ( (tan(alpha1)*sin(beta1)-tan(alpha2)*sin(beta2))*
        //											(tan(alpha1)*sin(beta1)-tan(alpha2)*sin(beta2)) + 
        //											(tan(alpha1)*cos(beta1)-tan(alpha2)*cos(beta2))*
        //											(tan(alpha1)*cos(beta1)-tan(alpha2)*cos(beta2)) );
    }

    // 2. Get the average delay, for all datapoints where the correlation is > 0.9
    double* goodDelays = new double[m_length];
    int nGoodPoints = 0;
    for (int k = m_firstDataPoint; k < m_firstDataPoint + m_arrayLength; ++k)
    {
        if (corr[k] > 0.9)
        {
            goodDelays[nGoodPoints++] = delays[k];
        }
    }
    if (nGoodPoints < 50)
    {
        delete[] goodDelays;
        return 1; // if the measurement does not contain enough many good points
    }

    // 3. Calculate the wind-speed as the distance divided by the average delay
    double averageDelay = Average(goodDelays, nGoodPoints);
    double stdDelay = Std(goodDelays, nGoodPoints);

    windField.SetWindSpeed(distance / averageDelay, Meteorology::MeteorologySource::DualBeamMeasurement);

    // Estimate the error in the wind-speed
    windField.SetWindSpeedError(fabs(relativePlumeHeight.m_plumeAltitudeError / averageDelay) + fabs(relativePlumeHeight.m_plumeAltitude * stdDelay / (averageDelay * averageDelay)));

    // Also set the time for which the measurement is valid.
    CDateTime validFrom = m_startTime;
    CDateTime validTo = m_stopTime;
    validFrom.Decrement(m_userSettings.m_dualBeam_ValidTime / 2);
    validTo.Increment(m_userSettings.m_dualBeam_ValidTime / 2);
    windField.SetValidTimeFrame(validFrom, validTo);

    return 0;
}

RETURN_CODE CWindSpeedCalculator::CalculateCorrelation(const novac::CString& evalLog1, const novac::CString& evalLog2)
{
    CDateTime time;
    Meteorology::WindField wf;
    novac::CString errorMessage;
    WindSpeedMeasurement::CWindSpeedCalculator::CMeasurementSeries* series[2];
    size_t scanIndex[2];
    double delay;

    // 1. Read the evaluation-logs
    std::vector<FileHandler::CEvaluationLogFileHandler> reader;

    FileHandler::CEvaluationLogFileHandler reader1(m_log, evalLog1.std_str(), m_userSettings.m_molecule);
    if (RETURN_CODE::SUCCESS != reader1.ReadEvaluationLog())
        return RETURN_CODE::FAIL;
    reader.push_back(reader1);

    FileHandler::CEvaluationLogFileHandler reader2(m_log, evalLog2.std_str(), m_userSettings.m_molecule);
    if (RETURN_CODE::SUCCESS != reader2.ReadEvaluationLog())
        return RETURN_CODE::FAIL;
    reader.push_back(reader2);

    // 2. Find the wind-speed measurement series in the log-files
    for (int k = 0; k < 2; ++k)
    {
        for (scanIndex[k] = 0; scanIndex[k] < reader[k].m_scan.size(); ++scanIndex[k])
        {
            if (novac::IsWindMeasurement(reader[k].m_scan[scanIndex[k]]))
            {
                break;
            }
        }
        if (scanIndex[k] == reader[k].m_scan.size())
        {
            return RETURN_CODE::FAIL; // <-- no wind-speed measurement found
        }
    }

    // 2a. Check the measurement series to make sure that they are the same 
    //		length and that they can be combined
    if (reader[0].m_scan[scanIndex[0]].GetEvaluatedNum() != reader[1].m_scan[scanIndex[1]].GetEvaluatedNum())
    {
        errorMessage.Format("Cannot correlate series %s and %s. The number of spectra is not same", (const char*)evalLog1, (const char*)evalLog2);
        ShowMessage(errorMessage);
        return RETURN_CODE::FAIL;
    }

    // 2b. Find the start and stop-time of the measurement
    Evaluation::CScanResult& scan = reader[0].m_scan[scanIndex[0]];
    scan.GetStartTime(0, m_startTime);
    scan.GetStopTime(scan.GetEvaluatedNum() - 1, m_stopTime);

    // 3. Create the wind-speed measurement series
    for (int k = 0; k < 2; ++k)
    {
        // 3a. The scan we're looking at
        Evaluation::CScanResult& secondScan = reader[k].m_scan[scanIndex[k]];

        // 3c. The length of the measurement
        int	length = secondScan.GetEvaluatedNum();

        // 3d. Allocate the memory for the series
        series[k] = new CMeasurementSeries(length);

        // 3e. Copy the relevant data in the scan
        for (int i = 0; i < length; ++i)
        {
            secondScan.GetStartTime(i, time);

            // get the column value
            series[k]->column[i] = secondScan.GetColumn(i, 0);

            // calculate the time-difference between the start of the
            //	time-series and this measurement
            series[k]->time[i] = 3600.0 * (time.hour - m_startTime.hour) +
                60.0 * (time.minute - m_startTime.minute) +
                1.0 * (time.second - m_startTime.second);
        }
    }

    // 4. Perform the correlation calculations...

    // 4a. Calculate the correlation, assuming that series[0] is the upwind series
    if (RETURN_CODE::SUCCESS != CalculateDelay(delay, series[0], series[1], m_settings))
    {
        ShowMessage("Failed to correlate time-series, no windspeed could be derived");

        // Tell the world that we've tried to make a correlation calculation but failed
        // Evaluation::CScanResult &scan = reader[0].m_scan[scanIndex[0]];
        //scan.GetDate(0, date);

        // Clean up and return
        delete series[0];		delete series[1];
        return RETURN_CODE::FAIL;
    }

    // 4b. Calculate the average correlation
    double avgCorr1 = Average(corr + m_firstDataPoint, m_arrayLength);

    // 4c. Calculate the correlation, assuming that series[1] is the upwind series
    if (RETURN_CODE::SUCCESS != CalculateDelay(delay, series[1], series[0], m_settings))
    {
        ShowMessage("Failed to correlate time-series, no windspeed could be derived");

        // Tell the world that we've tried to make a correlation calculation but failed
        // Evaluation::CScanResult &scan = reader[0].m_scan[scanIndex[0]];
        // scan.GetDate(0, date);
        // PostWindMeasurementResult(0, 0, 0, m_startTime, stopTime, scan.GetSerial());

        // Clean up and return
        delete series[0];		delete series[1];
        return RETURN_CODE::FAIL;
    }

    // 4d. Calculate the average correlation
    double avgCorr2 = Average(corr + m_firstDataPoint, m_arrayLength);

    // 4e. Use the result which gave the higest correlation
    if (avgCorr1 > avgCorr2)
    {
        if (RETURN_CODE::SUCCESS != CalculateDelay(delay, series[0], series[1], m_settings))
        {
            // this should never happen
            ShowMessage("Failed to correlate time-series, no windspeed could be derived");
            delete series[0];		delete series[1];
            return RETURN_CODE::FAIL;
        }
    }

    // 5. Return the results of the calculation

//	WriteWindMeasurementLog(calc, evalLog1, reader[0].m_scan[scanIndex[0]], INSTR_GOTHENBURG);

    // 6. Clean up a little bit.
    delete series[0];
    delete series[1];

    return RETURN_CODE::SUCCESS;
}

/** Calculate the correlation between the two time-series found in the
        given evaluation-file. */
RETURN_CODE CWindSpeedCalculator::CalculateCorrelation_Heidelberg(const novac::CString& evalLog)
{
    WindSpeedMeasurement::CWindSpeedCalculator::CMeasurementSeries* series[2];
    Meteorology::WindField wf;
    CDateTime time;
    size_t scanIndex;
    double delay;

    // 1. Read the evaluation-log
    FileHandler::CEvaluationLogFileHandler reader(m_log, evalLog.std_str(), m_userSettings.m_molecule);
    if (RETURN_CODE::SUCCESS != reader.ReadEvaluationLog())
    {
        return RETURN_CODE::FAIL;
    }

    // 2. Find the wind-speed measurement series in the log-files
    for (scanIndex = 0; scanIndex < reader.m_scan.size(); ++scanIndex)
    {
        if (novac::IsWindMeasurement_Heidelberg(reader.m_scan[scanIndex]))
        {
            break;
        }
    }
    if (scanIndex == reader.m_scan.size())
    {
        return RETURN_CODE::FAIL; // <-- no wind-speed measurement found
    }

    // 3. Create the wind-speed measurement series

    // 3a. The scan we're looking at
    Evaluation::CScanResult& scan = reader.m_scan[scanIndex];

    // 3b. The start-time of the whole measurement
    scan.GetStartTime(0, m_startTime);
    scan.GetStopTime(scan.GetEvaluatedNum() - 1, m_stopTime);

    // 3c. The length of the measurement
    int	length = scan.GetEvaluatedNum();

    // 3d. Allocate the memory for the series
    series[0] = new CWindSpeedCalculator::CMeasurementSeries(length / 2);
    series[1] = new CWindSpeedCalculator::CMeasurementSeries(length / 2);

    // 3e. Copy the relevant data in the scan
    for (int k = 0; k < length; ++k)
    {
        scan.GetStartTime(k, time);

        series[k % 2]->column[k / 2] = scan.GetColumn(k, 0);

        // Save the time difference
        series[k % 2]->time[k / 2] =
            3600 * (time.hour - m_startTime.hour) +
            60 * (time.minute - m_startTime.minute) +
            (time.second - m_startTime.second);
    }

    // 3f. Adjust the settings to have the correct angle
    int midpoint = (int)(length / 2);
    double d1 = scan.GetScanAngle(midpoint) - scan.GetScanAngle(midpoint + 1);
    double d2 = scan.GetScanAngle2(midpoint) - scan.GetScanAngle2(midpoint + 1);
    m_settings.angleSeparation = sqrt(d1 * d1 + d2 * d2);

    // 4. Perform the correlation calculations...

    // 4a. Calculate the correlation, assuming that series[0] is the upwind series
    if (RETURN_CODE::SUCCESS != CalculateDelay(delay, series[0], series[1], m_settings))
    {
        ShowMessage("Failed to correlate time-series, no windspeed could be derived");
        delete series[0];		delete series[1];
        return RETURN_CODE::FAIL;
    }

    // 4b. Calculate the average correlation
    double avgCorr1 = Average(corr + m_firstDataPoint, m_length);

    // 4c. Calculate the correlation, assuming that series[1] is the upwind series
    if (RETURN_CODE::SUCCESS != CalculateDelay(delay, series[1], series[0], m_settings))
    {
        ShowMessage("Failed to correlate time-series, no windspeed could be derived");
        delete series[0];		delete series[1];
        return RETURN_CODE::FAIL;
    }

    // 4d. Calculate the average correlation
    double avgCorr2 = Average(corr + m_firstDataPoint, m_length);

    // 4e. Use the result which gave the higest correlation
    if (avgCorr1 > avgCorr2)
    {
        if (RETURN_CODE::SUCCESS != CalculateDelay(delay, series[0], series[1], m_settings))
        {
            // this should never happen
            ShowMessage("Failed to correlate time-series, no windspeed could be derived");
            delete series[0];		delete series[1];
            return RETURN_CODE::FAIL;
        }
    }

    // 5. Write the results of our calculations to file
//	WriteWindMeasurementLog(calc, evalLog, reader.m_scan[scanIndex], INSTR_HEIDELBERG);

    // 6. Clean up a little bit.
    delete series[0];
    delete series[1];

    return RETURN_CODE::SUCCESS;
}

/** Writes the header of a dual-beam wind speed log file to the given
    file. */
void CWindSpeedCalculator::WriteWindSpeedLogHeader(const novac::CString& fileName)
{
    CDateTime now;

    if (Filesystem::IsExistingFile(fileName))
        return; // don't write the header if the file already exists..

    // get the current time
    now.SetToNow();

    FILE* f = fopen(fileName, "w");
    if (f == NULL)
        return;

    // Write the header and the starting comments
    fprintf(f, "# This is result of the dual-beam calculations of the NOVAC Post Processing Program \n");
    fprintf(f, "#   File generated on %04d.%02d.%02d at %02d:%02d:%02d \n\n", now.year, now.month, now.day, now.hour, now.minute, now.second);

    fprintf(f, "#StartTimeOfMeas\tWsValidFrom\tWsValidTo\tCalculatedWindSpeed\tWindSpeedError\tPlumeHeightUsed\tPlumeHeightError\n");

    fclose(f);
    return;
}

/** Appends a dual-beam wind speed result to the given file */
void CWindSpeedCalculator::AppendResultToFile(const novac::CString& fileName, const CDateTime& startTime,
    const Configuration::CInstrumentLocation& /*location*/,
    const Geometry::PlumeHeight& plumeHeight,
    Meteorology::WindField& windField)
{
    CDateTime validFrom, validTo;

    // try to open the file
    FILE* f = fopen(fileName, "a");
    if (f == NULL)
        return;

    // Write the result
    fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t",
        startTime.year, startTime.month, startTime.day,
        startTime.hour, startTime.minute, startTime.second);

    // valid from and valid to
    windField.GetValidTimeFrame(validFrom, validTo);

    fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t",
        validFrom.year, validFrom.month, validFrom.day,
        validFrom.hour, validFrom.minute, validFrom.second);

    fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t",
        validTo.year, validTo.month, validTo.day,
        validTo.hour, validTo.minute, validTo.second);

    // the calculated wind speed and it's error
    fprintf(f, "%.2lf\t", windField.GetWindSpeed());
    fprintf(f, "%.2lf\t", windField.GetWindSpeedError());

    // the plume altitude used, and it's error
    fprintf(f, "%.2lf\t", plumeHeight.m_plumeAltitude);
    fprintf(f, "%.2lf\n", plumeHeight.m_plumeAltitudeError);



    // remember to close the file before we return
    fclose(f);
    return;


}