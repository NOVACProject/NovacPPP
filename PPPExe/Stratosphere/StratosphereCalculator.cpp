#include "StratosphereCalculator.h"
#include <SpectralEvaluation/Geometry.h>

#include <PPPLib/Logging.h>
#include <PPPLib/Evaluation/ScanResult.h>
#include <PPPLib/File/EvaluationLogFileHandler.h>
#include <PPPLib/File/Filesystem.h>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

// This is the configuration of the network
#include <PPPLib/Configuration/NovacPPPConfiguration.h>

#include <Poco/Path.h>
#include <cmath>

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


CStratosphereCalculator::CStratosphereCalculator(
    novac::ILogger& log,
    const Configuration::CNovacPPPConfiguration& setup,
    const Configuration::CUserConfiguration& userSettings)
    : m_log(log), m_setup(setup), m_userSettings(userSettings)
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
    //fileName.Format("%s%cStratosphere.txt", (const char*)m_userSettings.m_outputDirectory,  Poco::Path::separator());
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
    novac::CString mainFitWindowName = novac::CString(m_userSettings.m_fitWindowsToUse[m_userSettings.m_mainFitWindow]);
    Molecule specie = Molecule(m_userSettings.m_molecule);

    // loop through all the evaluation log files
    while (p != results.end())
    {
        const Evaluation::CExtendedScanResult& result = (Evaluation::CExtendedScanResult&)*(p++);

        std::string evalLogfileToRead;
        for (size_t k = 0; k < m_userSettings.m_nFitWindowsToUse; ++k)
        {
            if (Equals(result.m_fitWindowName[k], mainFitWindowName))
            {
                evalLogfileToRead = result.m_evalLogFile[k];
                break;
            }
        }
        if (evalLogfileToRead.size() < 3)
        {
            continue; // file not found...
        }

        // Read the evaluation-log file
        FileHandler::CEvaluationLogFileHandler reader(m_log, evalLogfileToRead, m_userSettings.m_molecule);
        if (RETURN_CODE::SUCCESS != reader.ReadEvaluationLog())
        {
            continue;
        }

        // Get a handle to the result
        Evaluation::CScanResult& scanResult = reader.m_scan[0];

        // loop through the measurements in this scan and insert them into the correct measurement day
        CMeasurementDay measurementDay;
        scanResult.GetSkyStartTime(measurementDay.day);

        for (size_t k = 0; k < scanResult.GetEvaluatedNum(); ++k)
        {
            try
            {
                // make sure that we only insert zenith measurements
                if (fabs(scanResult.GetScanAngle(k)) > 1.0)
                    continue;

                CZenithMeasurement meas;
                scanResult.GetStartTime(k, meas.time);
                meas.column = scanResult.GetColumn(k, specie);

                // find the location of this instrument
                auto instrLocation = m_setup.GetInstrumentLocation(scanResult.GetSerial(), meas.time);
                CGPSData location = instrLocation.GpsData();

                // calculate the AMF
                meas.AMF = GetAMF_ZenithMeasurement(location, meas.time);

                measurementDay.measList.push_back(meas);
            }
            catch (novac::NotFoundException& ex)
            {
                ShowMessage(ex.message);
            }
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

double CStratosphereCalculator::GetAMF_ZenithMeasurement(const CGPSData& location, const CDateTime& gmtTime)
{
    novac::SolarPosition sun = novac::GetSunPosition(gmtTime, location);

    return 1.0 / std::cos(sun.zenithAngle.InRadians());
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
    fileName.Format("%s%cStratosphere.txt", (const char*)m_userSettings.m_outputDirectory, Poco::Path::separator());

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

