#include <PPPLib/Geometry/PlumeDataBase.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <SpectralEvaluation/VectorUtils.h>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

#undef min
#undef max

using namespace Geometry;

extern Configuration::CUserConfiguration g_userSettings;// <-- The settings of the user

// ----------- THE SUB-CLASS CPlumeData --------------
CPlumeDataBase::CPlumeData::CPlumeData()
{
    this->altitude = 1000.0;
    this->altitudeError = 1000.0;
    this->altitudeSource = Meteorology::MeteorologySource::Default;

    this->validFrom = novac::CDateTime(0, 0, 0, 0, 0, 0);
    this->validTo = novac::CDateTime(9999, 12, 31, 23, 59, 59);
}

CPlumeDataBase::CPlumeData::CPlumeData(const CPlumeDataBase::CPlumeData& p)
{
    this->altitude = p.altitude;
    this->altitudeError = p.altitudeError;
    this->altitudeSource = p.altitudeSource;

    this->validFrom = p.validFrom;
    this->validTo = p.validTo;
}

CPlumeDataBase::CPlumeData::~CPlumeData()
{
}

CPlumeDataBase::CPlumeData& CPlumeDataBase::CPlumeData::operator =(const CPlumeDataBase::CPlumeData& p)
{
    this->altitude = p.altitude;
    this->altitudeError = p.altitudeError;
    this->altitudeSource = p.altitudeSource;

    this->validFrom = p.validFrom;
    this->validTo = p.validTo;

    return *this;
}

// --------- THE CLASS CPlumeDataBase ----------

CPlumeDataBase::CPlumeDataBase(void)
{
}

CPlumeDataBase::~CPlumeDataBase(void)
{
    m_dataBase.clear();
}

/** Retrieves the plume height at a given time.
    @param time - the time for which the wind field should be retrieved
    @param plumeHeight - will on successful return be filled with the information
        about the plume height at the requested time.
    @return true if the wind field could be retrieved, otherwise false.
    */
bool CPlumeDataBase::GetPlumeHeight(const novac::CDateTime& time, CPlumeHeight& plumeHeight) const
{
    std::list <CPlumeData> validData;

    // There can be more than one piece of wind-information valid for this given moment
    //	extract the ones which are valid and put them into the list 'validData'
    std::list <CPlumeData>::const_iterator pos = m_dataBase.begin();
    while (pos != m_dataBase.end())
    {
        const CPlumeData& data = (CPlumeData&)*pos;

        if ((data.validFrom < time || data.validFrom == time) && (time < data.validTo || time == data.validTo))
        {
            validData.push_back(CPlumeData(data));
        }

        ++pos; // go to the next element in the list
    }
    if (validData.size() == 0)
        return false; // there's no datapoints which are valid for the given time...

    // If there's only one time, then return that one
    if (validData.size() == 1)
    {
        CPlumeData& data = (CPlumeData&)*validData.begin();
        plumeHeight.m_plumeAltitude = data.altitude;
        plumeHeight.m_plumeAltitudeError = data.altitudeError;
        plumeHeight.m_plumeAltitudeSource = data.altitudeSource;
        plumeHeight.m_validFrom = data.validFrom;
        plumeHeight.m_validTo = data.validTo;
        return true;
    }

    // If there are several, then the priority is to take the one which is 
    //	calculated (GeometryCalculationTwoInstruments) over the others. If there are 
    //	several calculated then use their average value
    std::list <CPlumeData> calculatedData_2instr;
    std::list <CPlumeData> calculatedData_1instr;
    pos = validData.begin();
    while (pos != validData.end())
    {
        CPlumeData& data = (CPlumeData&)*pos;

        if (Meteorology::MeteorologySource::GeometryCalculationTwoInstruments == data.altitudeSource)
        {
            calculatedData_2instr.push_back(CPlumeData(data));
        }
        else if (Meteorology::MeteorologySource::GeometryCalculationSingleInstrument == data.altitudeSource)
        {
            calculatedData_1instr.push_back(CPlumeData(data));
        }

        ++pos; // go to the next element in the list
    }
    if (calculatedData_2instr.size() == 1)
    {
        // There's only one geometry calculation. Return that one...
        CPlumeData& data = (CPlumeData&)*calculatedData_2instr.begin();
        plumeHeight.m_plumeAltitude = data.altitude;
        plumeHeight.m_plumeAltitudeError = data.altitudeError;
        plumeHeight.m_plumeAltitudeSource = Meteorology::MeteorologySource::GeometryCalculationTwoInstruments;
        plumeHeight.m_validFrom = data.validFrom;
        plumeHeight.m_validTo = data.validTo;
        return true;
    }
    else if (calculatedData_2instr.size() > 1)
    {
        double avgAltitude, altitudeError;
        CalculateAverageHeight(calculatedData_2instr, avgAltitude, altitudeError);

        plumeHeight.m_plumeAltitude = avgAltitude;
        plumeHeight.m_plumeAltitudeError = altitudeError;
        plumeHeight.m_plumeAltitudeSource = Meteorology::MeteorologySource::GeometryCalculationTwoInstruments;
        plumeHeight.m_validFrom = time;
        plumeHeight.m_validTo = time;
        return true;
    }
    else if (calculatedData_1instr.size() == 1)
    {
        // There's only one geometry calculation from a single instrument. Return that one...
        CPlumeData& data = (CPlumeData&)*calculatedData_1instr.begin();
        plumeHeight.m_plumeAltitude = data.altitude;
        plumeHeight.m_plumeAltitudeError = data.altitudeError;
        plumeHeight.m_plumeAltitudeSource = Meteorology::MeteorologySource::GeometryCalculationSingleInstrument;
        plumeHeight.m_validFrom = data.validFrom;
        plumeHeight.m_validTo = data.validTo;
        return true;
    }
    else if (calculatedData_1instr.size() > 1)
    {
        double avgAltitude, altitudeError;
        CalculateAverageHeight(calculatedData_1instr, avgAltitude, altitudeError);

        plumeHeight.m_plumeAltitude = avgAltitude;
        plumeHeight.m_plumeAltitudeError = altitudeError;
        plumeHeight.m_plumeAltitudeSource = Meteorology::MeteorologySource::GeometryCalculationSingleInstrument;
        plumeHeight.m_validFrom = time;
        plumeHeight.m_validTo = time;
        return true;
    }

    // If we get this far, then there's no calculated data...
    //	use the average of the available data..
    double avgAltitude, altitudeError;
    CalculateAverageHeight(calculatedData_2instr, avgAltitude, altitudeError);

    plumeHeight.m_plumeAltitude = avgAltitude;
    plumeHeight.m_plumeAltitudeError = altitudeError;
    plumeHeight.m_plumeAltitudeSource = Meteorology::MeteorologySource::GeometryCalculationTwoInstruments;
    plumeHeight.m_validFrom = time;
    plumeHeight.m_validTo = time;

    return false;
}

/** Inserts a plume height into the database */
void CPlumeDataBase::InsertPlumeHeight(const CPlumeHeight& plumeHeight)
{
    CPlumeData data;

    // generate a copy of the CPlumeHeight
    data.altitude = (float)plumeHeight.m_plumeAltitude;
    data.altitudeError = (float)plumeHeight.m_plumeAltitudeError;
    data.altitudeSource = plumeHeight.m_plumeAltitudeSource;
    data.validFrom = plumeHeight.m_validFrom;
    data.validTo = plumeHeight.m_validTo;

    // insert the copy into the database
    m_dataBase.push_back(data);
}

/** Inserts a calculated plume height into the database */
void CPlumeDataBase::InsertPlumeHeight(const CGeometryResult& geomResult)
{
    CPlumeData data;
    novac::CDateTime validFrom = geomResult.m_averageStartTime;
    novac::CDateTime validTo = geomResult.m_averageStartTime;

    // make the result of the measurement valid within +- 'g_userSettings.m_calcGeometryValidTime'/2 minutes of the measurement occasion
    validFrom.Decrement(g_userSettings.m_calcGeometryValidTime / 2);
    validTo.Increment(g_userSettings.m_calcGeometryValidTime / 2);

    // generate a copy of the CGeometryResult
    data.altitude = (float)geomResult.m_plumeAltitude;
    data.altitudeError = (float)geomResult.m_plumeAltitudeError;
    data.altitudeSource = geomResult.m_calculationType;
    data.validFrom = validFrom;
    data.validTo = validTo;

    // insert the copy into the database
    m_dataBase.push_back(data);
}

/** Writes the contents of this database to file.
    @return 0 on success. */
int CPlumeDataBase::WriteToFile(const novac::CString& /*fileName*/) const
{

    return 1;
}

void CPlumeDataBase::CalculateAverageHeight(const std::list <CPlumeData>& plumeList, double& averageAltitude, double& altitudeError) const
{
    std::vector<double> plumeAltitudes;
    std::vector<double> plumeAltitudeErrors;

    // loop through the altitudes to extract the average and the errors
    std::list <CPlumeData>::const_iterator pos = plumeList.begin();
    while (pos != plumeList.end())
    {
        const CPlumeData& data = (CPlumeData&)*pos;
        plumeAltitudes.push_back(data.altitude);
        plumeAltitudeErrors.push_back(data.altitudeError);

        ++pos; // go to the next element in the list
    }

    // Calculate the error
    averageAltitude = Average(plumeAltitudes);
    altitudeError = std::max(Stdev(plumeAltitudes), Average(plumeAltitudeErrors));
}
