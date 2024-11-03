#pragma once

#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Evaluation/FitParameter.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>

#include <PPPLib/PPPLib.h>
#include <PPPLib/Flux/FluxResult.h>
#include <PPPLib/MFC/CString.h>

namespace novac
{
struct SpectrometerModel;
}

namespace Evaluation
{
/** <b>CScanResult</b> is a class designed to handle the results
    from evaluating a set of spectra (e.g. a scan).

    It contains a set of CEvaluationResult's which describes the
    result of evaluating each spectrum.

    The CScanResult also handles information about the whole set
    of evaluation results such as the offset or the calculated flux of the scan,
    or a judgement wheather each evaluated spectrum is judged to be an ok
    spectrum or not. */
class CScanResult : public novac::BasicScanEvaluationResult
{
public:
    CScanResult() = default;

    CScanResult(const CScanResult&);
    CScanResult& operator=(const CScanResult& s2);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The calculated flux and the parameters used to calculate the flux */
    Flux::FluxResult m_flux;

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Appends the result to the list of calculated results */
    void AppendResult(const novac::CEvaluationResult& evalRes, const novac::CSpectrumInfo& specInfo);

    /** Removes the spectrum number 'specIndex' from the list of calculated results */
    int RemoveResult(size_t specIndex);

    /** Intializes the memory arrays to have, initially, space for 'specNum' spectra. */
    void InitializeArrays(size_t specNum);

    /** Retrieves the evaluation result for spectrum number
        'specIndex' from the list of calculated results.
            @return - NULL if specIndex is out of bounds... */
    const novac::CEvaluationResult* GetResult(size_t specIndex) const;

    /** Adds spectrum number 'specIndex' into the list of spectra in the .pak -file
            which are corrupted and could not be evaluated */
    void MarkAsCorrupted(size_t specIndex);

    /** Stores the information about the sky-spectrum used */
    void SetSkySpecInfo(const novac::CSpectrumInfo& skySpecInfo);

    /** Stores the information about the dark-spectrum used */
    void SetDarkSpecInfo(const novac::CSpectrumInfo& darkSpecInfo);

    /** Stores the information about the offset-spectrum used */
    void SetOffsetSpecInfo(const novac::CSpectrumInfo& offsetSpecInfo);

    /** Stores the information about the dark-current-spectrum used */
    void SetDarkCurrentSpecInfo(const novac::CSpectrumInfo& darkCurSpecInfo);

    /** Check the last spectrum point for goodness of fit.
        The parameters 'deltaLimit', 'upperLimit' and 'lowerLimit' are for
        development purposes only. */
    bool CheckGoodnessOfFit(const novac::CSpectrumInfo& info, const novac::SpectrometerModel* spectrometer, float chi2Limit = -1, float upperLimit = -1, float lowerLimit = -1);

    /** Check spectrum number 'index' for goodness of fit.
        The parameters 'deltaLimit', 'upperLimit' and 'lowerLimit' are for
        development purposes only. */
    bool CheckGoodnessOfFit(const novac::CSpectrumInfo& info, size_t index, const novac::SpectrometerModel* spectrometer, float chi2Limit = -1, float upperLimit = -1, float lowerLimit = -1);

    /** Returns the calculated flux */
    double GetFlux() const { return m_flux.m_flux; }

    const Flux::FluxResult& GetFluxResult() const { return m_flux; }

    /** Returns true if the automatic judgment considers this flux
        measurement to be a good measurement */
    bool IsFluxOk() const { return (m_flux.m_fluxQualityFlag != FluxQuality::Red); }

    /** Set the flux to the given value. ONLY USED FOR READING EVALUATION-LOGS */
    void SetFlux(double flux) { this->m_flux.m_flux = flux; }

    /** returns the temperature of the system at the time of the measurement */
    double GetTemperature() const { return m_skySpecInfo.m_temperature; }

    /** returns the number of spectra evaluated */
    size_t GetEvaluatedNum() const { return m_specNum; }

    /** Returns the latitude of the system */
    double GetLatitude() const;

    /** Returns the longitude of the system */
    double GetLongitude() const;

    /** Returns the altitude of the system */
    double GetAltitude() const;

    /** Returns the compass-direction of the system */
    double GetCompass() const;

    /** Returns the cone angle of the scanning instrument */
    double GetConeAngle() const;

    /** Returns the pitch (tilt) of the scanning instrument */
    double GetPitch() const;

    /** Returns the roll (scan-angle offset) of the scanning instrument */
    double GetRoll() const;

    /** Returns the battery-voltage of the sky spectrum */
    float GetBatteryVoltage() const;

    /** Returns the serial-number of the spectrometer that collected this scan */
    std::string GetSerial() const;

    /** returns the goodness of fit for the fitting of the evaluated
            spectrum number 'index'.
        This function is the complementary of IsBad(size_t index)*/
    int  IsOk(size_t index) const { return m_spec[index].IsOK(); }

    /** returns the goodness of fit for the fitting of the evaluated
            spectrum number 'index'.
        This function is the complementary of IsOK(size_t index). */
    int  IsBad(size_t index) const { return m_spec[index].IsBad(); }

    /** returns true if the evaluated spectrum number 'index' is marked
            as deleted */
    int  IsDeleted(size_t index) const { return m_spec[index].IsDeleted(); }

    /** Marks the desired spectrum with the supplied mark_flag.
        Mark flag must be MARK_BAD_EVALUATION, or MARK_DELETED
        @return true on success. */
    bool MarkAs(size_t index, int MARK_FLAG);

    /** Removes the desired mark from the desired spectrum
        Mark flag must be MARK_BAD_EVALUATION, or MARK_DELETED
        @return true on success. */
    bool RemoveMark(size_t index, int MARK_FLAG);

    /** Returns a reference to the desired spectrum info-structure */
    const novac::CSpectrumInfo& GetSpectrumInfo(size_t index) const;

    /** Returns a reference to the spectrum info-structure of the sky-spectrum used */
    const novac::CSpectrumInfo& GetSkySpectrumInfo() const;

    /** Returns a reference to the spectrum info-structure of the dark-spectrum used */
    const novac::CSpectrumInfo& GetDarkSpectrumInfo() const;

    /** returns the scan angle of evaluated spectrum number 'index'.
        @param index - the zero based index into the list of  evaluated spectra */
    float GetScanAngle(size_t index) const { return (IsValidSpectrumIndex(index)) ? m_specInfo[index].m_scanAngle : 0; }

    /** returns the azimuth angle (the angle of the second motor) of
            evaluated spectrum number 'index'.
        @param index - the zero based index into the list of  evaluated spectra */
    float GetScanAngle2(size_t index) const { return (IsValidSpectrumIndex(index)) ? m_specInfo[index].m_scanAngle2 : 0; }

    /** returns the time and date (UMT) when evaluated spectrum number
            'index' was started.
        @param index - the zero based index into the list of evaluated spectra.
        @return SUCCESS if the index is valid */
    RETURN_CODE GetStartTime(size_t index, novac::CDateTime& time) const;

    /** returns the time and date (UMT) when the sky-spectrum was started. */
    void GetSkyStartTime(novac::CDateTime& t) const;
    novac::CDateTime GetSkyStartTime() const;

    /** returns the time and date (UMT) when evaluated spectrum number 'index'
            was stopped.
        @param index - the zero based index into the list of evaluated spectra.
        @return SUCCESS if the index is valid */
    RETURN_CODE GetStopTime(size_t index, novac::CDateTime& time) const;

    /** returns the evaluated column for specie number 'specieNum' and
            spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of evaluated
            spectra.*/
    double GetColumn(size_t spectrumNum, size_t specieNum) const;
    double GetColumn(size_t spectrumNum, novac::Molecule& mol) const;

    /** returns the error for the evaluated column for specie number
            'specieNum' and spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of evaluated
            spectra.*/
    double GetColumnError(size_t spectrumNum, size_t specieNum) const;

    /** returns the SHIFT parameter for specie number 'specieNum' and
            spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of
            evaluated spectra.*/
    double GetShift(size_t spectrumNum, size_t specieNum) const;

    /** returns the error for the SHIFT parameter for specie number
            'specieNum' and spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of
            evaluated spectra.*/
    double GetShiftError(size_t spectrumNum, size_t specieNum) const;

    /** returns the SQUEEZE parameter for specie number 'specieNum' and
            spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of
            evaluated spectra.*/
    double GetSqueeze(size_t spectrumNum, size_t specieNum) const;

    /** returns the error for the SQUEEZE parameter for specie number
            'specieNum' and spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of
            evaluated spectra.*/
    double GetSqueezeError(size_t spectrumNum, size_t specieNum) const;

    /** @return the delta of the fit for spectrum number 'spectrumNum'
        @param spectrumNum - the spectrum number (zero-based) for
            which the delta value is desired */
    double GetDelta(size_t spectrumNum) const;

    /** @return the chi-square of the fit for spectrum number 'spectrumNum'
        @param spectrumNum - the spectrum number (zero-based) for
            which the delta value is desired */
    double GetChiSquare(size_t spectrumNum) const;

    /** returns the number of spectra averaged to get evaluated
            spectrum number 'spectrumNum'
        @return the number of spectra averaged. */
    long GetSpecNum(size_t spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_numSpec : 0; }

    /** returns the expsure time of evaluated spectrum number 'spectrumNum'
            in ms  */
    long GetExposureTime(size_t spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_exposureTime : 0; }

    /** returns the peak intensity of evaluated spectrum number 'spectrumNum'
        (the maximum intensity of the whole spectrum). */
    float GetPeakIntensity(size_t spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_peakIntensity : 0; }

    /** returns the fit intensity of evaluated spectrum number 'spectrumNum'
        (the maximum intensity int the fit region of the spectrum). */
    float GetFitIntensity(size_t spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_fitIntensity : 0; }

    /** returns the electronic offset in spectrum number 'spectrumNum' */
    float GetElectronicOffset(size_t spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_offset : 0; }

    /** returns the number of species that were used in the evaluation of a
        given spectrum */
    int GetSpecieNum(size_t spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? (int)m_spec[spectrumNum].m_referenceResult.size() : 0; }

    /** returns the specie name */
    const std::string GetSpecieName(size_t spectrumNum, size_t specieNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_spec[spectrumNum].m_referenceResult[specieNum].m_specieName : "NA"; }

    /** Sets the type of the instrument used */
    void SetInstrumentType(novac::NovacInstrumentType type);

    /** Returns the type of the instrument used */
    novac::NovacInstrumentType GetInstrumentType() const;

private:

    // ----------------------------------------------------------------------
    // -------------------- PRIVATE METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** makes a sanity check of the parameters and returns fit parameter number 'index'.
        @param specIndex - the zero based into the list of evaluated spectra.
        @param specieIndex - the zero based into the list of species to evaluate for.
        @param fitParameter - a parameter to return.
        @return NaN if any parameter is wrong */
    double GetFitParameter(size_t spectrumNum, size_t specieIndex, FIT_PARAMETER parameter) const;

    /** returns true if the given index is a valid spectrum index */
    inline bool IsValidSpectrumIndex(size_t spectrumIndex) const { return spectrumIndex < m_specNum; }
};
}