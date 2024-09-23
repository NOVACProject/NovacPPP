#pragma once

#include "WindSpeedMeasSettings.h"
#include <PPPLib/PPPLib.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Geometry/PlumeHeight.h>
#include <PPPLib/Meteorology/WindField.h>
#include <PPPLib/MFC/CString.h>

namespace WindSpeedMeasurement
{

/** The <b>CWindSpeedCalculator</b> class contains the basic
        algorithms for calculating wind speeds from measured data series
        of column variations using the correlation between the two
        data series. The idea being that the two time series have been
        measured in such a way that one of them measured on a more up-wind
        position along the plume than the other one. Assuming a certain
        altitude of the plume, the speed of the plume can be calculated
        by calculating the temporal delay between the two time series. */
class CWindSpeedCalculator
{
public:

    class CMeasurementSeries
    {
    public:
        CMeasurementSeries();			// <-- Creates an empty measurement series
        CMeasurementSeries(int len);	// <-- Creates a measurement series of length 'len'
        ~CMeasurementSeries();
        RETURN_CODE SetLength(int len); // <-- changes the length of the measurment series to 'len'
        double	AverageColumn(int from, int to) const; // <-- calculated the average column value between 'from' and 'to'
        double	SampleInterval();		// <-- calculates and returns the average time between two measurements
        double* column;
        double* time;
        long		length;
    };

    CWindSpeedCalculator(void);
    ~CWindSpeedCalculator(void);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The settings for how to perform the correlation calculations */
    CWindSpeedMeasSettings	m_settings;

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------


    /** Calculates the wind speed from the two time series found in the given evaluation-log files.
        If the instrument is a Heidelberg instrument then the second parameter will be ignored.
        @param evalLog1 - the full path and file name of the evaluation-log file
            containing the first of the two time series to correlate
        @param evalLog2 - the full path and file name of the evaluation-log file
            containing the first of the two time series to correlate.
            This will be ignored if the instrument is a Heidelberg type.
        @param location - the location of the instrument at the time of the measurement.
        @param plumeHeight - the altitude of the plume at the time of the wind speed
            measurement. This is the full altitude of the plume, in meters above sea level.
        @windField - will on successful return be filled with the derived wind-speed.
            The time-period that the measurement is valid for will be filled in.
            NOTE THAT THE WIND-DIRECTION WILL NOT BE CALCULATED AND WILL THUS NOT MAKE
            ANY SENSE...
        @return 0 on success, else non-zero.
    */
    int CalculateWindSpeed(const novac::CString& evalLog1, const novac::CString& evalLog2,
        const Configuration::CInstrumentLocation& location,
        const Geometry::CPlumeHeight& plumeHeight,
        Meteorology::CWindField& windField);

    /** Writes the header of a dual-beam wind speed log file to the given
        file. */
    void WriteWindSpeedLogHeader(const novac::CString& fileName);

    /** Appends a dual-beam wind speed result to the given file */
    void AppendResultToFile(const novac::CString& fileName, const novac::CDateTime& startTime,
        const Configuration::CInstrumentLocation& location,
        const Geometry::CPlumeHeight& plumeHeight,
        Meteorology::CWindField& windField);

private:
    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------
    /** The calculated values. These will be filled in after a call to 'CalculateDelay'
            Before that they are null and cannot be used. The length of these arrays are 'm_length' */
    double* shift, * corr, * used, * delays;
    int			m_arrayLength;
    int			m_length;
    int			m_firstDataPoint;

    /** This is the start-time and the stop time of the measurement.
        This is filled in after a call to 'CalculateCorrelation' */
    novac::CDateTime	m_startTime;
    novac::CDateTime	m_stopTime;

    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------

    /** Calculate the correlation between the two time-series found in the
            given evaluation-files.
        @param evalLog1 - the full path and file name of the evaluation-log file
            containing the first of the two time series to correlate
        @param evalLog2 - the full path and file name of the evaluation-log file
            containing the first of the two time series to correlate
        The results of the calculations will be filled into the buffers 'shift',
        'corr', 'used' and 'delays'
         */
    RETURN_CODE CalculateCorrelation(const novac::CString& evalLog1, const novac::CString& evalLog2);

    /** Calculate the correlation between the two time-series found in the
            given evaluation-file.
        @param evalLog - the full path and file name of the evaluation-log file
            containing the two time series to correlate
        The results of the calculations will be filled into the buffers 'shift',
        'corr', 'used' and 'delays'
        */
    RETURN_CODE CalculateCorrelation_Heidelberg(const novac::CString& evalLog);


    /** Calculate the time delay between the two provided time series
        @param delay - The return value of the function. Will be set to the time
            delay between the two data series, in seconds. Will be positive if the
            'upWindColumns' comes temporally before the 'downWindColumns', otherwise negative.
        @param upWindSerie - the measurement for the more upwind time series
        @param downWindSerie - the measurement for the more downwind time series
        @param settings - The settings for how the calculation should be done
    */
    RETURN_CODE CalculateDelay(double& delay,
        const CMeasurementSeries* upWindSerie,
        const CMeasurementSeries* downWindSerie,
        const CWindSpeedMeasSettings& settings);

    /** Intializes the arrays 'shift', 'corr', 'used' and 'delays' before they are used*/
    void InitializeArrays();

    /** Performs a low pass filtering on the supplied measurement series.
            The number of iterations in the filtering is given by 'nIterations'
            if nIterations is zero, nothing will be done. */
    static RETURN_CODE LowPassFilter(const CMeasurementSeries* series, CMeasurementSeries* result, unsigned int nIterations);

    /** Shifts the vector 'shortVector' against the vector 'longVector' and returns the
                shift for which the correlation between the two is highest.
                The length of the longVector must be larger than the length of the short vector! */
    static RETURN_CODE FindBestCorrelation(
        const double* longVector, unsigned long longLength,
        const double* shortVector, unsigned long shortLength,
        unsigned int maximumShift,
        double& highestCorr, int& bestShift);


    /** Calculates the correlation between the two vectors 'x' and 'y', both of length 'length'
            @return - the correlation between the two vectors. */
    static double	correlation(const double* x, const double* y, long length);
};
}