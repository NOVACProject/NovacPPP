#include "StratosphereCalculator.h"

#include "../Evaluation/ScanResult.h"
#include <PPPLib/Molecule.h>
#include "../Common/Common.h"
#include "../Common/EvaluationLogFileHandler.h"
#include <PPPLib/File/Filesystem.h>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

// This is the configuration of the network
#include <PPPLib/Configuration/NovacPPPConfiguration.h>

#include <Poco/Path.h>


extern Configuration::CNovacPPPConfiguration        g_setup;	   // <-- The settings
extern Configuration::CUserConfiguration			g_userSettings;// <-- The settings of the user

using namespace Stratosphere;
using namespace novac;

// ------------- The private class Zenith Measurement -------------
CStratosphereCalculator::CZenithMeasurement::CZenithMeasurement()
{
    this->column = 0;
    this->AMF = 0;
    this->time = CDateTime(0, 0, 0, 0, 0, 0);
}

CStratosphereCalculator::CZenithMeasurement::~CZenithMeasurement()
{

}

CStratosphereCalculator::CZenithMeasurement& CStratosphereCalculator::CZenithMeasurement::operator =(const CStratosphereCalculator::CZenithMeasurement& z)
{
    this->time = z.time;
    this->column = z.column;
    this->AMF = z.AMF;
    return *this;
}

// ------------- The private class Measurement Day -------------

CStratosphereCalculator::CMeasurementDay::CMeasurementDay()
{

}

CStratosphereCalculator::CMeasurementDay::~CMeasurementDay()
{

}

CStratosphereCalculator::CMeasurementDay& CStratosphereCalculator::CMeasurementDay::operator =(const CStratosphereCalculator::CMeasurementDay& m)
{
    this->day = m.day;
    this->measList.clear();

    std::list<CZenithMeasurement>::const_iterator p = m.measList.begin();
    while (p != m.measList.end())
    {
        CZenithMeasurement& zm = (CZenithMeasurement&)*(p++);
        this->measList.push_back(zm);
    }

    return *this;
}


CStratosphereCalculator::CStratosphereCalculator(void)
{
}

CStratosphereCalculator::~CStratosphereCalculator(void)
{
}


/** Calculates the stratospheric columns based on the results int the
    supplied list of evaluation logs.
*/
void CStratosphereCalculator::CalculateVCDs(const std::list <Evaluation::CExtendedScanResult>& results)
{
    double vcd, vcdErr, s0;
    novac::CString fileName;

    // remove the old stratospheric output-file
    // TODO: ImplementMe
    //fileName.Format("%s%cStratosphere.txt", (const char*)g_userSettings.m_outputDirectory,  Poco::Path::separator());
    //DeleteFile(fileName);


    // make the list of measurements
    this->BuildMeasurementList(results);

    // write the results to files
    std::list<CMeasurementDay>::const_iterator p = m_measurementDays.begin();
    while (p != m_measurementDays.end())
    {
        CMeasurementDay& m = (CMeasurementDay&)*(p++);

        CalculateVCD(m, vcd, vcdErr, s0);

        WriteResultToFile(m, vcd, vcdErr, s0);
    }
}

/** Builds the list of measurements. I.e. takes the supplied
    result files and builds the list 'm_measurementDays' */
void CStratosphereCalculator::BuildMeasurementList(const std::list <Evaluation::CExtendedScanResult>& results)
{
    std::list<Evaluation::CExtendedScanResult>::const_iterator p = results.begin();
    novac::CString mainFitWindowName = novac::CString(g_userSettings.m_fitWindowsToUse[g_userSettings.m_mainFitWindow]);
    novac::CString evalLogfileToRead;
    CMolecule specie = CMolecule(g_userSettings.m_molecule);
    Configuration::CInstrumentLocation instrLocation;

    // loop through all the evaluation log files
    while (p != results.end())
    {
        const Evaluation::CExtendedScanResult& result = (Evaluation::CExtendedScanResult&)*(p++);

        evalLogfileToRead.Format("");
        for (int k = 0; k < g_userSettings.m_nFitWindowsToUse; ++k)
        {
            if (Equals(result.m_fitWindowName[k], mainFitWindowName))
            {
                evalLogfileToRead.Format(result.m_evalLogFile[k]);
                k = 1000;
                continue;
            }
        }
        if (evalLogfileToRead.GetLength() < 3)
            continue; // file not found...

        // REad the evaluation-log file
        FileHandler::CEvaluationLogFileHandler reader;
        reader.m_evaluationLog.Format(evalLogfileToRead);
        if (RETURN_CODE::SUCCESS != reader.ReadEvaluationLog())
            continue;

        // Get a handle to the result
        Evaluation::CScanResult& scanResult = reader.m_scan[0];

        // loop through the measurements in this scan and insert them into the correct measurement day
        CMeasurementDay measurementDay;
        scanResult.GetSkyStartTime(measurementDay.day);

        for (int k = 0; k < scanResult.GetEvaluatedNum(); ++k)
        {
            // make sure that we only insert zenith measurements
            if (fabs(scanResult.GetScanAngle(k)) > 1.0)
                continue;

            CZenithMeasurement meas;
            scanResult.GetStartTime(k, meas.time);
            meas.column = scanResult.GetColumn(k, specie);

            // find the location of this instrument
            if (g_setup.GetInstrumentLocation(scanResult.GetSerial(), meas.time, instrLocation))
                continue;
            CGPSData location = CGPSData(instrLocation.m_latitude, instrLocation.m_longitude, instrLocation.m_altitude);

            // calculate the AMF
            meas.AMF = GetAMF_ZenithMeasurement(location, meas.time);

            measurementDay.measList.push_back(meas);
        }
        InsertIntoMeasurementList(measurementDay);
    }

}

/** Inserts the results in the supplied measurement day structure
    into the appropriate place in the list of measurements */
void CStratosphereCalculator::InsertIntoMeasurementList(const CMeasurementDay& mday)
{
    CDateTime measDate = CDateTime(mday.day.year, mday.day.month, mday.day.day, 0, 0, 0);
    CMeasurementDay r;
    r = mday; // make a local copy of the results


    std::list<CMeasurementDay>::iterator p = m_measurementDays.begin();
    while (p != m_measurementDays.end())
    {
        CMeasurementDay& d = (CMeasurementDay&)*p;
        if (d.day == measDate)
        {
            // insert the data into this day
            d.measList.insert(d.measList.end(), r.measList.begin(), r.measList.end());
            return;

        }
        else if (measDate < d.day)
        {
            // insert the result at the position before this day
            p = m_measurementDays.insert(p, r);
            return;
        }

        // move on to the next measurement day
        ++p;
    }

    // we've passed the whole list without finding anything that's larger than this day
    //	insert the result as a new measurement day in the end of the list
    r.day = measDate;
    m_measurementDays.push_back(r);

    return;

}

/** Retrieves the Air Mass Factor for a zenith measuremnt performed
    at the given location and at the given time of day (UTC). */
double CStratosphereCalculator::GetAMF_ZenithMeasurement(const CGPSData& location, const CDateTime& gmtTime)
{
    double SZA, SAZ;

    if (RETURN_CODE::SUCCESS != Common::GetSunPosition(gmtTime, location.m_latitude, location.m_longitude, SZA, SAZ))
        return 1.0;

    return 1.0 / cos(DEGREETORAD * SZA);
}


/** Calculates the VCD for a given measurement day
    @param VCD - will on return be the calculated VCD
    @param VCDErr - will on return be the estimated error in the calculated VCD
    @param S0 - will on return be the calculated slant column in the sky-spectrum
    */
void CStratosphereCalculator::CalculateVCD(CMeasurementDay& /*measDay*/, double& /*VCD*/, double& /*VCDErr*/, double& /*S0*/)
{


}

/** Writes the given results to the output file */
void CStratosphereCalculator::WriteResultToFile(CMeasurementDay& measDay, double /*VCD*/, double /*VCDErr*/, double /*S0*/)
{
    bool writeHeader = false;
    novac::CString fileName;

    // the name of the output file
    fileName.Format("%s%cStratosphere.txt", (const char*)g_userSettings.m_outputDirectory, Poco::Path::separator());

    if (!Filesystem::IsExistingFile(fileName))
        writeHeader = true;

    FILE* f = fopen(fileName, "a");
    if (f == NULL)
        return;

    if (writeHeader)
    {
        fprintf(f, "Start\tColumn\tAMF");
        fprintf(f, "\n");
    }

    std::list<CZenithMeasurement>::const_iterator p = measDay.measList.begin();
    while (p != measDay.measList.end())
    {
        CZenithMeasurement& m = (CZenithMeasurement&)*(p++);

        fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t", m.time.year, m.time.month, m.time.day, m.time.hour, m.time.minute, m.time.second);
        fprintf(f, "%.5e\t", m.column);
        fprintf(f, "%.2lf\t", m.AMF);
        fprintf(f, "\n");
    }

    fclose(f);
    return;
}

