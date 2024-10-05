#pragma once

#include <PPPLib/Evaluation/ExtendedScanResult.h>
#include <SpectralEvaluation/GPSData.h>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Log.h>

#include <list>

namespace Configuration
{
    class CNovacPPPConfiguration;
    class CUserConfiguration;
}

namespace Stratosphere
{
class CStratosphereCalculator
{
public:
    CStratosphereCalculator(
        novac::ILogger& log,
        const Configuration::CNovacPPPConfiguration& setup,
        const Configuration::CUserConfiguration& userSettings);
    ~CStratosphereCalculator(void);

    // -----------------------------------------------------------
    // ---------------------- PUBLIC DATA ------------------------
    // -----------------------------------------------------------


    // -----------------------------------------------------------
    // --------------------- PUBLIC METHODS ----------------------
    // -----------------------------------------------------------

    /** Calculates the stratospheric columns based on the results int the
        supplied list of evaluation logs.
    */
    void CalculateVCDs(const std::list <Evaluation::CExtendedScanResult>& results);

private:
    // -----------------------------------------------------------
    // ---------------------- PRIVATE DATA -----------------------
    // -----------------------------------------------------------
    class CZenithMeasurement
    {
    public:
        CZenithMeasurement();
        ~CZenithMeasurement();
        novac::CDateTime time;	// the time of the measurement
        double column; // Sn
        double AMF; // the calculated air mass factor at the time of measurement
        CZenithMeasurement& operator=(const CZenithMeasurement& z);
    };

    class CMeasurementDay
    {
    public:
        CMeasurementDay();
        ~CMeasurementDay();
        novac::CDateTime day; // the date of the measurement
        std::list <CZenithMeasurement> measList;
        CMeasurementDay& operator=(const CMeasurementDay& m);
    };

    /** The set of measurements that we have */
    std::list <CMeasurementDay> m_measurementDays;

    novac::ILogger& m_log;

    const Configuration::CNovacPPPConfiguration& m_setup;

    const Configuration::CUserConfiguration& m_userSettings;

    // -----------------------------------------------------------
    // --------------------- PRIVATE METHODS ---------------------
    // -----------------------------------------------------------

    /** Builds the list of measurements. I.e. takes the supplied
        result files and builds the list 'm_measurementDays' */
    void BuildMeasurementList(const std::list <Evaluation::CExtendedScanResult>& results);

    /** Inserts the results in the supplied measurement day structure
        into the appropriate place in the list of measurements */
    void InsertIntoMeasurementList(const CMeasurementDay& mday);

    /** Retrieves the Air Mass Factor for a zenith measuremnt performed
        at the given location and at the given time of day (UTC). */
    double GetAMF_ZenithMeasurement(const novac::CGPSData& location, const novac::CDateTime& gmtTime);

    /** Calculates the VCD for a given measurement day
        @param VCD - will on return be the calculated VCD
        @param VCDErr - will on return be the estimated error in the calculated VCD
        @param S0 - will on return be the calculated slant column in the sky-spectrum
        */
    void CalculateVCD(CMeasurementDay& measDay, double& VCD, double& VCDErr, double& S0);

    /** Writes the given results to the output file */
    void WriteResultToFile(CMeasurementDay& measDay, double VCD, double VCDErr, double S0);

};
}
