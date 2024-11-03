#include <PPPLib/Evaluation/ScanResult.h>
#include <PPPLib/VolcanoInfo.h>
#include <PPPLib/Geometry/GeometryCalculator.h>
#include <PPPLib/Logging.h>
#include <PPPLib/Flux/FluxCalculator.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/Geometry.h>
#include <algorithm>
#include <cmath>

using namespace Evaluation;
using namespace novac;

extern novac::CVolcanoInfo g_volcanoes; // <-- A list of all known volcanoes

CScanResult::CScanResult(const CScanResult& s2) :
    m_flux(s2.m_flux)
{
    this->m_specNum = s2.m_specNum;
    this->m_spec = s2.m_spec;
    this->m_specInfo = s2.m_specInfo;
    this->m_corruptedSpectra = s2.m_corruptedSpectra;

    this->m_plumeProperties = s2.m_plumeProperties;

    this->m_skySpecInfo = s2.m_skySpecInfo;
    this->m_darkSpecInfo = s2.m_darkSpecInfo;

    this->m_plumeProperties = s2.m_plumeProperties;
    this->m_instrumentType = s2.m_instrumentType;
    this->m_measurementMode = s2.m_measurementMode;
}

CScanResult& CScanResult::operator=(const CScanResult& s2)
{
    // The calculated flux and offset
    this->m_flux = s2.m_flux;

    this->m_plumeProperties = s2.m_plumeProperties;

    this->m_plumeProperties = s2.m_plumeProperties;

    this->m_spec = s2.m_spec;
    this->m_specInfo = s2.m_specInfo;
    this->m_corruptedSpectra = s2.m_corruptedSpectra;
    this->m_specNum = s2.m_specNum;

    this->m_skySpecInfo = s2.m_skySpecInfo;
    this->m_darkSpecInfo = s2.m_darkSpecInfo;

    this->m_measurementMode = s2.m_measurementMode;
    this->m_instrumentType = s2.m_instrumentType;

    return *this;
}

void CScanResult::InitializeArrays(size_t specNum)
{
    if (specNum > 1024)
    {
        return;
    }

    m_spec.reserve(specNum);
    m_specInfo.reserve(specNum);
}

void CScanResult::AppendResult(const CEvaluationResult& evalRes, const CSpectrumInfo& specInfo)
{
    m_spec.push_back(CEvaluationResult(evalRes));
    m_specInfo.push_back(CSpectrumInfo(specInfo));

    // Increase the number of spectra in this result-set.
    ++m_specNum;
}

void CScanResult::MarkAsCorrupted(size_t specIndex)
{
    m_corruptedSpectra.push_back(specIndex);
}

const CEvaluationResult* CScanResult::GetResult(size_t specIndex) const
{
    if (specIndex >= m_specNum)
        return nullptr; // not a valid index

    return &m_spec.at(specIndex);
}

int CScanResult::RemoveResult(size_t specIndex)
{
    if (specIndex >= m_specNum)
        return 1; // not a valid index

    // Remove the desired value
    auto it = m_specInfo.begin() + static_cast<std::int64_t>(specIndex);
    m_specInfo.erase(it);

    auto it2 = m_spec.begin() + static_cast<std::int64_t>(specIndex);
    m_spec.erase(it2);

    // Decrease the number of values in the list
    m_specNum -= 1;

    return 0;
}

/** Stores the information about the sky-spectrum used */
void CScanResult::SetSkySpecInfo(const CSpectrumInfo& skySpecInfo)
{
    this->m_skySpecInfo = skySpecInfo;
}

/** Stores the information about the dark-spectrum used */
void CScanResult::SetDarkSpecInfo(const CSpectrumInfo& darkSpecInfo)
{
    this->m_darkSpecInfo = darkSpecInfo;
}

/** Stores the information about the offset-spectrum used */
void CScanResult::SetOffsetSpecInfo(const CSpectrumInfo& offsetSpecInfo)
{
    this->m_offsetSpecInfo = offsetSpecInfo;
}

/** Stores the information about the dark-current-spectrum used */
void CScanResult::SetDarkCurrentSpecInfo(const CSpectrumInfo& darkCurSpecInfo)
{
    this->m_darkCurSpecInfo = darkCurSpecInfo;
}

bool CScanResult::CheckGoodnessOfFit(const CSpectrumInfo& info, const SpectrometerModel* spectrometer, float chi2Limit, float upperLimit, float lowerLimit)
{
    return CheckGoodnessOfFit(info, m_specNum - 1, spectrometer, chi2Limit, upperLimit, lowerLimit);
}

bool CScanResult::CheckGoodnessOfFit(const CSpectrumInfo& info, size_t index, const SpectrometerModel* spectrometer, float chi2Limit, float upperLimit, float lowerLimit)
{
    if (index >= m_specNum)
    {
        return false;
    }

    // remember the electronic offset (NB. this is not same as the scan-offset)
    //  m_specInfo[index].m_offset    = (float)offsetLevel;

    return m_spec[index].CheckGoodnessOfFit(info, spectrometer, chi2Limit, upperLimit, lowerLimit);
}

double CScanResult::GetColumn(size_t spectrumNum, size_t specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, COLUMN);
}

double CScanResult::GetColumn(size_t spectrumNum, Molecule& molec) const
{
    const int index = this->GetSpecieIndex(molec.name);
    if (index == -1)
    {
        return NOT_A_NUMBER;
    }

    return GetFitParameter(spectrumNum, static_cast<size_t>(index), COLUMN);
}

double CScanResult::GetColumnError(size_t spectrumNum, size_t specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, COLUMN_ERROR);
}

double CScanResult::GetShift(size_t spectrumNum, size_t specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, SHIFT);
}

double CScanResult::GetShiftError(size_t spectrumNum, size_t specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, SHIFT_ERROR);
}

double CScanResult::GetSqueeze(size_t spectrumNum, size_t specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, SQUEEZE);
}

double CScanResult::GetSqueezeError(size_t spectrumNum, size_t specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, SQUEEZE_ERROR);
}

double CScanResult::GetDelta(size_t spectrumNum) const
{
    return this->m_spec[spectrumNum].m_delta;
}

double CScanResult::GetChiSquare(size_t spectrumNum) const
{
    return this->m_spec[spectrumNum].m_chiSquare;
}

/** Returns the desired fit parameter */
double CScanResult::GetFitParameter(size_t specIndex, size_t specieIndex, FIT_PARAMETER parameter) const
{
    if (specIndex >= m_specNum)
        return 0.0f;

    if (specieIndex >= this->m_spec[specIndex].m_referenceResult.size())
        return 0.0f;

    switch (parameter)
    {
    case COLUMN:        return this->m_spec[specIndex].m_referenceResult[specieIndex].m_column;
    case COLUMN_ERROR:  return this->m_spec[specIndex].m_referenceResult[specieIndex].m_columnError;
    case SHIFT:         return this->m_spec[specIndex].m_referenceResult[specieIndex].m_shift;
    case SHIFT_ERROR:   return this->m_spec[specIndex].m_referenceResult[specieIndex].m_shiftError;
    case SQUEEZE:       return this->m_spec[specIndex].m_referenceResult[specieIndex].m_squeeze;
    case SQUEEZE_ERROR: return this->m_spec[specIndex].m_referenceResult[specieIndex].m_squeezeError;
    case DELTA:         return this->m_spec[specIndex].m_delta;
    default:            return 0.0f;
    }
}

const CSpectrumInfo& CScanResult::GetSpectrumInfo(size_t index) const
{
    return m_specInfo[index];
}

/** Returns a reference to the spectrum info-structure of the sky-spectrum used */
const CSpectrumInfo& CScanResult::GetSkySpectrumInfo() const
{
    return m_skySpecInfo;
}

/** Returns a reference to the spectrum info-structure of the dark-spectrum used */
const CSpectrumInfo& CScanResult::GetDarkSpectrumInfo() const
{
    return m_darkSpecInfo;
}

bool CScanResult::MarkAs(size_t index, int MARK_FLAG)
{
    if (!IsValidSpectrumIndex(index))
        return false;

    return m_spec[index].MarkAs(MARK_FLAG);
}

bool CScanResult::RemoveMark(size_t index, int MARK_FLAG)
{
    if (!IsValidSpectrumIndex(index))
        return false;

    return m_spec[index].RemoveMark(MARK_FLAG);
}

/** Returns the latitude of the system */
double CScanResult::GetLatitude() const
{
    for (unsigned int k = 0; k < m_specNum; ++k)
    {
        const CSpectrumInfo& info = m_specInfo[k];
        if (fabs(info.m_gps.m_latitude) > 1e-2)
            return info.m_gps.m_latitude;
    }
    return 0.0;
}

/** Returns the longitude of the system */
double CScanResult::GetLongitude() const
{
    for (unsigned int k = 0; k < m_specNum; ++k)
    {
        const CSpectrumInfo& info = m_specInfo[k];
        if (fabs(info.m_gps.m_longitude) > 1e-2)
            return info.m_gps.m_longitude;
    }
    return 0.0;
}
/** Returns the altitude of the system */
double CScanResult::GetAltitude() const
{
    for (unsigned int k = 0; k < m_specNum; ++k)
    {
        const CSpectrumInfo& info = m_specInfo[k];
        if (fabs(info.m_gps.m_altitude) > 1e-2)
            return info.m_gps.m_altitude;
    }
    return 0.0;
}

/** Returns the compass-direction of the system */
double CScanResult::GetCompass() const
{
    if (m_specNum == 0)
        return 0.0;

    const CSpectrumInfo& info = m_specInfo.front();
    return info.m_compass;
}

/** Returns the battery-voltage of the sky spectrum */
float CScanResult::GetBatteryVoltage() const
{
    if (m_specNum == 0)
        return 0.0;

    const CSpectrumInfo& info = m_specInfo.front();
    return info.m_batteryVoltage;
}

/** Returns the cone angle of the scanning instrument */
double CScanResult::GetConeAngle() const
{
    if (m_specNum == 0)
        return 90.0;

    const CSpectrumInfo& info = m_specInfo.front();
    return info.m_coneAngle;
}

/** Returns the pitch of the scanning instrument */
double CScanResult::GetPitch() const
{
    if (m_specNum == 0)
        return 90.0;

    const CSpectrumInfo& info = m_specInfo.front();
    return info.m_pitch;
}

/** Returns the roll of the scanning instrument */
double CScanResult::GetRoll() const
{
    if (m_specNum == 0)
        return 90.0;

    const CSpectrumInfo& info = m_specInfo.front();
    return info.m_roll;
}

std::string CScanResult::GetSerial() const
{
    for (unsigned int k = 0; k < m_specNum; ++k)
    {
        const CSpectrumInfo& info = m_specInfo[k];
        if (info.m_device.size() > 0)
        {
            return info.m_device;
        }
    }
    return "";
}


/** returns the time and date (UMT) when evaluated spectrum number 'index' was started.
    @param index - the zero based index into the list of evaluated spectra */
RETURN_CODE CScanResult::GetStartTime(size_t index, CDateTime& t) const
{
    if (!IsValidSpectrumIndex(index))
        return RETURN_CODE::FAIL;

    // The start-time
    t = m_specInfo[index].m_startTime;

    return RETURN_CODE::SUCCESS;
}

void CScanResult::GetSkyStartTime(CDateTime& t) const
{
    t = m_skySpecInfo.m_startTime;
    return;
}
CDateTime CScanResult::GetSkyStartTime() const
{
    return m_skySpecInfo.m_startTime;
}

RETURN_CODE CScanResult::GetStopTime(size_t index, CDateTime& t) const
{
    if (!IsValidSpectrumIndex(index))
        return RETURN_CODE::FAIL;

    t = m_specInfo[index].m_stopTime;

    return RETURN_CODE::SUCCESS;
}


/** Sets the type of the instrument used */
void CScanResult::SetInstrumentType(NovacInstrumentType type)
{
    this->m_instrumentType = type;
}

/** Sets the type of the instrument used */
NovacInstrumentType CScanResult::GetInstrumentType() const
{
    return this->m_instrumentType;
}
