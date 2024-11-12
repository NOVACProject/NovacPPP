#include <PPPLib/Geometry/PlumeDataBase.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <SpectralEvaluation/VectorUtils.h>

#undef min
#undef max

using namespace Geometry;

// ----------- THE SUB-CLASS PlumeData --------------
CPlumeDataBase::PlumeData::PlumeData()
{
    this->altitude = 1000.0;
    this->altitudeError = 1000.0;
    this->altitudeSource = Meteorology::MeteorologySource::Default;

    this->validFrom = novac::CDateTime(0, 0, 0, 0, 0, 0);
    this->validTo = novac::CDateTime(9999, 12, 31, 23, 59, 59);
}

CPlumeDataBase::PlumeData::PlumeData(const CPlumeDataBase::PlumeData& p)
{
    this->altitude = p.altitude;
    this->altitudeError = p.altitudeError;
    this->altitudeSource = p.altitudeSource;

    this->validFrom = p.validFrom;
    this->validTo = p.validTo;
}

CPlumeDataBase::PlumeData::~PlumeData()
{
}

CPlumeDataBase::PlumeData& CPlumeDataBase::PlumeData::operator =(const CPlumeDataBase::PlumeData& p)
{
    this->altitude = p.altitude;
    this->altitudeError = p.altitudeError;
    this->altitudeSource = p.altitudeSource;

    this->validFrom = p.validFrom;
    this->validTo = p.validTo;

    return *this;
}

// --------- THE CLASS CPlumeDataBase ----------

CPlumeDataBase::CPlumeDataBase(const Configuration::CUserConfiguration& userSettings)
    : m_userSettings(userSettings)
{
}

bool CPlumeDataBase::GetPlumeHeight(const novac::CDateTime& time, PlumeHeight& plumeHeight) const
{
    std::list <PlumeData> validData;

    // There can be more than one piece of wind-information valid for this given moment
    //	extract the ones which are valid and put them into the list 'validData'
    std::list <PlumeData>::const_iterator pos = m_dataBase.begin();
    while (pos != m_dataBase.end())
    {
        const PlumeData& data = (PlumeData&)*pos;

        if ((data.validFrom < time || data.validFrom == time) && (time < data.validTo || time == data.validTo))
        {
            validData.push_back(PlumeData(data));
        }

        ++pos; // go to the next element in the list
    }
    if (validData.size() == 0)
        return false; // there's no datapoints which are valid for the given time...

    // If there's only one time, then return that one
    if (validData.size() == 1)
    {
        PlumeData& data = (PlumeData&)*validData.begin();
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
    std::list <PlumeData> calculatedData_2instr;
    std::list <PlumeData> calculatedData_1instr;
    pos = validData.begin();
    while (pos != validData.end())
    {
        PlumeData& data = (PlumeData&)*pos;

        if (Meteorology::MeteorologySource::GeometryCalculationTwoInstruments == data.altitudeSource)
        {
            calculatedData_2instr.push_back(PlumeData(data));
        }
        else if (Meteorology::MeteorologySource::GeometryCalculationSingleInstrument == data.altitudeSource)
        {
            calculatedData_1instr.push_back(PlumeData(data));
        }

        ++pos; // go to the next element in the list
    }
    if (calculatedData_2instr.size() == 1)
    {
        // There's only one geometry calculation. Return that one...
        PlumeData& data = (PlumeData&)*calculatedData_2instr.begin();
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
        PlumeData& data = (PlumeData&)*calculatedData_1instr.begin();
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

void CPlumeDataBase::InsertPlumeHeight(const PlumeHeight& plumeHeight)
{
    PlumeData data;

    // generate a copy of the PlumeHeight
    data.altitude = (float)plumeHeight.m_plumeAltitude;
    data.altitudeError = (float)plumeHeight.m_plumeAltitudeError;
    data.altitudeSource = plumeHeight.m_plumeAltitudeSource;
    data.validFrom = plumeHeight.m_validFrom;
    data.validTo = plumeHeight.m_validTo;

    // insert the copy into the database
    m_dataBase.push_back(data);
}

void CPlumeDataBase::InsertPlumeHeight(const CGeometryResult& geomResult)
{
    PlumeData data;
    novac::CDateTime validFrom = geomResult.m_averageStartTime;
    novac::CDateTime validTo = geomResult.m_averageStartTime;

    // make the result of the measurement valid within +- 'm_userSettings.m_calcGeometryValidTime'/2 minutes of the measurement occasion
    validFrom.Decrement(m_userSettings.m_calcGeometryValidTime / 2);
    validTo.Increment(m_userSettings.m_calcGeometryValidTime / 2);

    // generate a copy of the CGeometryResult
    data.altitude = geomResult.m_plumeAltitude.Value();
    data.altitudeError = geomResult.m_plumeAltitudeError.Value();
    data.altitudeSource = geomResult.m_calculationType;
    data.validFrom = validFrom;
    data.validTo = validTo;

    // insert the copy into the database
    m_dataBase.push_back(data);
}

int CPlumeDataBase::WriteToFile(const novac::CString& /*fileName*/) const
{

    return 1;
}

void CPlumeDataBase::CalculateAverageHeight(const std::list <PlumeData>& plumeList, double& averageAltitude, double& altitudeError) const
{
    std::vector<double> plumeAltitudes;
    std::vector<double> plumeAltitudeErrors;

    // loop through the altitudes to extract the average and the errors
    std::list <PlumeData>::const_iterator pos = plumeList.begin();
    while (pos != plumeList.end())
    {
        const PlumeData& data = (PlumeData&)*pos;
        plumeAltitudes.push_back(data.altitude);
        plumeAltitudeErrors.push_back(data.altitudeError);

        ++pos; // go to the next element in the list
    }

    // Calculate the error
    averageAltitude = Average(plumeAltitudes);
    altitudeError = std::max(Stdev(plumeAltitudes), Average(plumeAltitudeErrors));
}
