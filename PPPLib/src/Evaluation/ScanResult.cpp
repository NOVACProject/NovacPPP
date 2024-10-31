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


void CScanResult::InitializeArrays(long specNum)
{
    if (specNum < 0 || specNum > 1024)
    {
        return;
    }

    m_spec.reserve(specNum);
    m_specInfo.reserve(specNum);
}

/** Appends the result to the list of calculated results */
int CScanResult::AppendResult(const CEvaluationResult& evalRes, const CSpectrumInfo& specInfo)
{

    // Append the evaluationresult to the end of the 'm_spec'-vector
    m_spec.push_back(CEvaluationResult(evalRes));

    // Append the spectral information to the end of the 'm_specInfo'-vector
    m_specInfo.push_back(CSpectrumInfo(specInfo));

    // Increase the numbers of spectra in this result-set.
    ++m_specNum;
    return 0;
}

const CEvaluationResult* CScanResult::GetResult(unsigned int specIndex) const
{
    if (specIndex >= m_specNum)
        return NULL; // not a valid index

    return &m_spec.at(specIndex);
}

void CScanResult::MarkAsCorrupted(unsigned int specIndex)
{
    m_corruptedSpectra.push_back(specIndex);
}

int CScanResult::GetCorruptedNum() const
{
    return (int)m_corruptedSpectra.size();
}

int CScanResult::RemoveResult(unsigned int specIndex)
{
    if (specIndex >= m_specNum)
        return 1; // not a valid index

    // Remove the desired value
    auto it = m_specInfo.begin() + specIndex;
    m_specInfo.erase(it);
    // m_specInfo.RemoveAt(specIndex, 1);

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

bool CScanResult::CheckGoodnessOfFit(const CSpectrumInfo& info, int index, const SpectrometerModel* spectrometer, float chi2Limit, float upperLimit, float lowerLimit)
{
    if (index < 0 || (unsigned int)index >= m_specNum)
    {
        return false;
    }

    // remember the electronic offset (NB. this is not same as the scan-offset)
    //  m_specInfo[index].m_offset    = (float)offsetLevel;

    return m_spec[index].CheckGoodnessOfFit(info, spectrometer, chi2Limit, upperLimit, lowerLimit);
}

double CScanResult::GetColumn(unsigned long spectrumNum, unsigned long specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, COLUMN);
}

double CScanResult::GetColumn(unsigned long spectrumNum, Molecule& molec) const
{
    int index = this->GetSpecieIndex(molec.name);
    if (index == -1)
        return NOT_A_NUMBER;
    else
        return GetFitParameter(spectrumNum, index, COLUMN);
}

double CScanResult::GetColumnError(unsigned long spectrumNum, unsigned long specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, COLUMN_ERROR);
}

double CScanResult::GetShift(unsigned long spectrumNum, unsigned long specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, SHIFT);
}

double CScanResult::GetShiftError(unsigned long spectrumNum, unsigned long specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, SHIFT_ERROR);
}

double CScanResult::GetSqueeze(unsigned long spectrumNum, unsigned long specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, SQUEEZE);
}

double CScanResult::GetSqueezeError(unsigned long spectrumNum, unsigned long specieNum) const
{
    return this->GetFitParameter(spectrumNum, specieNum, SQUEEZE_ERROR);
}

/** @return the delta of the fit for spectrum number 'spectrumNum'
      @param spectrumNum - the spectrum number (zero-based) for which the delta value is desired         */
double CScanResult::GetDelta(unsigned long spectrumNum) const
{
    return this->m_spec[spectrumNum].m_delta;
}

/** @return the chi-square of the fit for spectrum number 'spectrumNum'
      @param spectrumNum - the spectrum number (zero-based) for which the delta value is desired         */
double CScanResult::GetChiSquare(unsigned long spectrumNum) const
{
    return this->m_spec[spectrumNum].m_chiSquare;
}

/** Returns the desired fit parameter */
double CScanResult::GetFitParameter(unsigned long specIndex, unsigned long specieIndex, FIT_PARAMETER parameter) const
{
    if (specIndex < 0 || specIndex > m_specNum)
        return 0.0f;

    if (specieIndex < 0 || specieIndex > this->m_spec[specIndex].m_referenceResult.size())
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

const CSpectrumInfo& CScanResult::GetSpectrumInfo(unsigned long index) const
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

bool CScanResult::MarkAs(unsigned long index, int MARK_FLAG)
{
    if (!IsValidSpectrumIndex(index))
        return false;

    return m_spec[index].MarkAs(MARK_FLAG);
}

bool CScanResult::RemoveMark(unsigned long index, int MARK_FLAG)
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

/** Returns the name of the requested spectrum */
novac::CString CScanResult::GetName(int index) const
{
    if (!IsValidSpectrumIndex(index))
        return novac::CString("");

    const CSpectrumInfo& info = m_specInfo[index];
    return info.m_name;
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
RETURN_CODE CScanResult::GetStartTime(unsigned long index, CDateTime& t) const
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

/** returns the time and date (UMT) when evaluated spectrum number 'index' was stopped.
    @param index - the zero based index into the list of evaluated spectra.
        @return SUCCESS if the index is valid */
RETURN_CODE CScanResult::GetStopTime(unsigned long index, CDateTime& t) const
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
