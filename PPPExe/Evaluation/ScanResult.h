#pragma once

#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/Evaluation/FitParameter.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>

#include "../Flux/FluxResult.h"
#include <PPPLib/Molecule.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/MFC/CArray.h>

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
    CScanResult();

    CScanResult(const CScanResult&);
    CScanResult& operator=(const CScanResult& s2);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Appends the result to the list of calculated results */
    int AppendResult(const novac::CEvaluationResult& evalRes, const novac::CSpectrumInfo& specInfo);

    /** Removes the spectrum number 'specIndex' from the list of calcualted results */
    int RemoveResult(unsigned int specIndex);

    /** Intializes the memory arrays to have, initially, space for
        'specNum' spectra. */
    void InitializeArrays(long specNum);

    /** Retrieves the evaluation result for spectrum number
        'specIndex' from the list of calculated results.
            @return - NULL if specIndex is out of bounds... */
    const novac::CEvaluationResult* GetResult(unsigned int specIndex) const;

    /** Adds spectrum number 'specIndex' into the list of spectra in the .pak -file
            which are corrupted and could not be evaluated */
    void MarkAsCorrupted(unsigned int specIndex);

    /** Retrieves how many spectra are corrupted in the scan */
    int GetCorruptedNum() const;

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
    bool CheckGoodnessOfFit(const novac::CSpectrumInfo& info, int index, const novac::SpectrometerModel* spectrometer, float chi2Limit = -1, float upperLimit = -1, float lowerLimit = -1);

    /** Gets the offset of the scan. The offset is calculated as the average of the
      three lowest columns values (bad values are skipped). After this function has
      been called the actual offset can be retrieved by a call to 'GetOffset'.
      @param specie - The name of the specie for which the offset should be found.
      @return 0 on success. @return 1 - if any error occurs. */
    int CalculateOffset(const CMolecule& specie);

    /** Calculate the flux in this scan, using the supplied compass direction
            and coneAngle.
        The result is saved in the private parameter 'm_flux', whose vale can be
        retrieved by a call to 'GetFlux()'
        @param specie - the name of the spece to calculate the flux for
        @param wind - the wind speed and wind direction to use in the calculation
        @param relativePlumeHeight - the height of the plume, in meters above the instrument
            that should be used to calculate this flux.
        @param compass - the compass direction of the instrument
        @param coneAngle - the cone-angle of the instrument
        @param tilt - the tilt of the instrument.
        @return 0 if all is ok. @return 1 if any error occurs. */
    int CalculateFlux(const CMolecule& specie, const Meteorology::CWindField& wind, const Geometry::CPlumeHeight& relativePlumeHeight, double compass, double coneAngle = 90.0, double tilt = 0.0);

    /** Tries to find a plume in the last scan result. If the plume is found
            this function returns true, and the centre of the plume (in scanAngles)
            is given in 'plumeCentre', the width of the plume (in scanAngles)
            is given in 'plumeWidth' and the estimated completeness of the plume
            is given in 'plumeCompleteness' (ranging from 0.0 to 1.0)
            */
    bool CalculatePlumeCentre(const CMolecule& specie, novac::CPlumeInScanProperty& plumeProperties);

    /** Tries to find a plume in the last scan result. If the plume is found
        this function returns true. The result of the calculations is stored in
        the member-variables 'm_plumeCentre' and 'm_plumeCompleteness' */
    bool CalculatePlumeCentre(const CMolecule& specie);

    /** Checks the kind of measurement that we have here and sets the
        flag 'm_measurementMode' to the appropriate value... */
    MEASUREMENT_MODE CheckMeasurementMode();

    /** Checks the kind of measurement that we have here without
        setting the flag 'm_measurementMode' */
    MEASUREMENT_MODE GetMeasurementMode() const;

    /** Returns true if this is a flux measurement */
    bool IsFluxMeasurement();

    /** Returns true if this is a wind-speed measurement */
    bool IsWindMeasurement() const;

    /** Returns true if this is a wind-speed measurement of Gothenburg type */
    bool IsWindMeasurement_Gothenburg() const;

    /** Returns true if this is a wind-speed measurement of Heidelberg type */
    bool IsWindMeasurement_Heidelberg() const;

    /** Returns true if this is a stratospheric mode measurement */
    bool IsStratosphereMeasurement() const;

    /** Returns true if this is a direct-sun mode measurement */
    bool IsDirectSunMeasurement() const;

    /** Returns true if this is a lunar mode measurement */
    bool IsLunarMeasurement() const;

    /** Returns true if this is a composition mode measurement */
    bool IsCompositionMeasurement() const;

    /** Calculates the maximum good column value in the scan,
        corrected for the offset.
        NB!! The function 'CalculateOffset' must have been called
        before this function is called. */
    double GetMaxColumn(const novac::CString& specie) const;

    /** Returns the calculated flux */
    double GetFlux() const { return m_flux.m_flux; }

    const Flux::CFluxResult& GetFluxResult() const { return m_flux; }

    /** Returns true if the automatic judgment considers this flux
        measurement to be a good measurement */
    bool IsFluxOk() const { return (m_flux.m_fluxQualityFlag != FLUX_QUALITY_RED); }

    /** Set the flux to the given value. ONLY USED FOR READING EVALUATION-LOGS */
    void SetFlux(double flux) { this->m_flux.m_flux = flux; }

    /** returns the offset of the measurement */
    double GetOffset() const { return m_plumeProperties.offset; }

    /** returns the temperature of the system at the time of the measurement */
    double GetTemperature() const { return m_skySpecInfo.m_temperature; }

    /** Fills in the supplied CPlumeInScanProperty object with the calculated properties of this scan */
    void GetCalculatedPlumeProperties(novac::CPlumeInScanProperty& properties) const { properties = m_plumeProperties; }

    /** Returns the calculated plume-centre position */
    double GetCalculatedPlumeCentre(int motor = 0) const;

    /** Returns the calculated plume edges */
    void GetCalculatedPlumeEdges(double& lowEdge, double& highEdge) const;

    /** Returns the calculated plume-completeness */
    double GetCalculatedPlumeCompleteness() const { return m_plumeProperties.completeness; }

    /** Sets the offset of the measurement */
    void SetOffset(double offset) { m_plumeProperties.offset = offset; }

    /** returns the number of spectra evaluated */
    long  GetEvaluatedNum() const { return m_specNum; }

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

    /** Returns the name of the requested spectrum */
    novac::CString GetName(int index) const;

    /** Returns the serial-number of the spectrometer that collected this scan */
    novac::CString GetSerial() const;

    /** returns the goodness of fit for the fitting of the evaluated
            spectrum number 'index'.
        This function is the complementary of IsBad(unsigned long index)*/
    int  IsOk(unsigned long index) const { return m_spec[index].IsOK(); }

    /** returns the goodness of fit for the fitting of the evaluated
            spectrum number 'index'.
        This function is the complementary of IsOK(unsigned long index). */
    int  IsBad(unsigned long index) const { return m_spec[index].IsBad(); }

    /** returns true if the evaluated spectrum number 'index' is marked
            as deleted */
    int  IsDeleted(unsigned long index) const { return m_spec[index].IsDeleted(); }

    /** Marks the desired spectrum with the supplied mark_flag.
        Mark flag must be MARK_BAD_EVALUATION, or MARK_DELETED
        @return true on success. */
    bool MarkAs(unsigned long index, int MARK_FLAG);

    /** Removes the desired mark from the desired spectrum
        Mark flag must be MARK_BAD_EVALUATION, or MARK_DELETED
        @return true on success. */
    bool RemoveMark(unsigned long index, int MARK_FLAG);

    /** Returns a reference to the desired spectrum info-structure */
    const novac::CSpectrumInfo& GetSpectrumInfo(unsigned long index) const;

    /** Returns a reference to the spectrum info-structure of the sky-spectrum used */
    const novac::CSpectrumInfo& GetSkySpectrumInfo() const;

    /** Returns a reference to the spectrum info-structure of the dark-spectrum used */
    const novac::CSpectrumInfo& GetDarkSpectrumInfo() const;

    /** returns the scan angle of evaluated spectrum number 'index'.
        @param index - the zero based index into the list of  evaluated spectra */
    float GetScanAngle(unsigned long index) const { return (IsValidSpectrumIndex(index)) ? m_specInfo[index].m_scanAngle : 0; }

    /** returns the azimuth angle (the angle of the second motor) of
            evaluated spectrum number 'index'.
        @param index - the zero based index into the list of  evaluated spectra */
    float GetScanAngle2(unsigned long index) const { return (IsValidSpectrumIndex(index)) ? m_specInfo[index].m_scanAngle2 : 0; }

    /** returns the time and date (UMT) when evaluated spectrum number
            'index' was started.
        @param index - the zero based index into the list of evaluated spectra.
        @return SUCCESS if the index is valid */
    RETURN_CODE GetStartTime(unsigned long index, novac::CDateTime& time) const;

    /** returns the time and date (UMT) when the sky-spectrum was started. */
    void GetSkyStartTime(novac::CDateTime& t) const;
    novac::CDateTime GetSkyStartTime() const;

    /** returns the time and date (UMT) when evaluated spectrum number 'index'
            was stopped.
        @param index - the zero based index into the list of evaluated spectra.
        @return SUCCESS if the index is valid */
    RETURN_CODE GetStopTime(unsigned long index, novac::CDateTime& time) const;

    /** returns the evaluated column for specie number 'specieNum' and
            spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of evaluated
            spectra.*/
    double GetColumn(unsigned long spectrumNum, unsigned long specieNum) const;
    double GetColumn(unsigned long spectrumNum, CMolecule& mol) const;

    /** returns the error for the evaluated column for specie number
            'specieNum' and spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of evaluated
            spectra.*/
    double GetColumnError(unsigned long spectrumNum, unsigned long specieNum) const;

    /** returns the SHIFT parameter for specie number 'specieNum' and
            spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of
            evaluated spectra.*/
    double GetShift(unsigned long spectrumNum, unsigned long specieNum) const;

    /** returns the error for the SHIFT parameter for specie number
            'specieNum' and spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of
            evaluated spectra.*/
    double GetShiftError(unsigned long spectrumNum, unsigned long specieNum) const;

    /** returns the SQUEEZE parameter for specie number 'specieNum' and
            spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of
            evaluated spectra.*/
    double GetSqueeze(unsigned long spectrumNum, unsigned long specieNum) const;

    /** returns the error for the SQUEEZE parameter for specie number
            'specieNum' and spectrum number 'specNum'
        @param specieNum - the zero based index into the list of species
            to evaluate for
        @param spectrumNum - the zero based index into the list of
            evaluated spectra.*/
    double GetSqueezeError(unsigned long spectrumNum, unsigned long specieNum) const;

    /** @return the delta of the fit for spectrum number 'spectrumNum'
        @param spectrumNum - the spectrum number (zero-based) for
            which the delta value is desired */
    double GetDelta(unsigned long spectrumNum) const;

    /** @return the chi-square of the fit for spectrum number 'spectrumNum'
        @param spectrumNum - the spectrum number (zero-based) for
            which the delta value is desired */
    double GetChiSquare(unsigned long spectrumNum) const;

    /** returns the number of spectra averaged to get evaluated
            spectrum number 'spectrumNum'
        @return the number of spectra averaged. */
    long GetSpecNum(unsigned long spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_numSpec : 0; }

    /** returns the expsure time of evaluated spectrum number 'spectrumNum'
            in ms  */
    long GetExposureTime(unsigned long spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_exposureTime : 0; }

    /** returns the peak intensity of evaluated spectrum number 'spectrumNum'
        (the maximum intensity of the whole spectrum). */
    float GetPeakIntensity(unsigned long spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_peakIntensity : 0; }

    /** returns the fit intensity of evaluated spectrum number 'spectrumNum'
        (the maximum intensity int the fit region of the spectrum). */
    float GetFitIntensity(unsigned long spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_fitIntensity : 0; }

    /** returns the electronic offset in spectrum number 'spectrumNum' */
    float GetElectronicOffset(unsigned long spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_specInfo[spectrumNum].m_offset : 0; }

    /** returns true if the spectra have been evaluated for the supplied specie.
        @param specie - a string containing the name of the specie to
            search for, e.g. "SO2" (case insensitive)*/
    bool IsEvaluatedSpecie(const novac::CString& specie) const { return (-1 != GetSpecieIndex(specie)); }

    /** returns the number of species that were used in the evaluation of a
        given spectrum */
    int GetSpecieNum(unsigned long spectrumNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? (int)m_spec[spectrumNum].m_referenceResult.size() : 0; }

    /** returns the specie name */
    const novac::CString GetSpecieName(unsigned long spectrumNum, unsigned long specieNum) const { return (IsValidSpectrumIndex(spectrumNum)) ? m_spec[spectrumNum].m_referenceResult[specieNum].m_specieName : 0; }

    /** Sets the type of the instrument used */
    void SetInstrumentType(INSTRUMENT_TYPE type);

    /** Returns the type of the instrument used */
    INSTRUMENT_TYPE GetInstrumentType() const;

    /** Getting the estimated geometrical error */
    double GetGeometricalError() const;

    /** Getting the scattering Error */
    double GetScatteringError() const;

    /** Getting the spectroscopical error */
    double GetSpectroscopicalError() const;

    /** Getting the estimated geometrical error */
    void SetGeometricalError(double err);

    /** Getting the scattering Error */
    void SetScatteringError(double err);

    /** Getting the spectroscopical error */
    void SetSpectroscopicalError(double err);
private:

    // ----------------------------------------------------------------------
    // --------------------- PRIVATE DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The calculated flux and the parameters used to
         calculate the flux */
    Flux::CFluxResult m_flux;

    /** The estimated error (in percent) in the geometrical setup for
        this flux-calculation. */
    double m_geomError;

    /** The estimated error (in percent) due to scattering inside or below the
        plume for this flux-calculation */
    double m_scatteringError;

    /** The estimated error (in percent) of the spectral results due
        to incertainties in cross-sections, slit-functions, stray-light etc. */
    double m_spectroscopyError;

    /** This contains the parameters of the plume that is seen
        in this scan, such as the completeness or the centre angle of
        the plume. */
    novac::CPlumeInScanProperty m_plumeProperties;

    /** The number of evaluations */
    unsigned long m_specNum;

    /** The type of the instrument used for this scan */
    INSTRUMENT_TYPE m_instrumentType;

    /** Flag to signal if this is a wind measurement, a scan, or something else. */
    MEASUREMENT_MODE m_measurementMode;

    // ----------------------------------------------------------------------
    // -------------------- PRIVATE METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** makes a sanity check of the parameters and returns fit parameter number 'index'.
        @param specIndex - the zero based into the list of evaluated spectra.
        @param specieIndex - the zero based into the list of species to evaluate for.
        @param fitParameter - a parameter to return.
        @return NaN if any parameter is wrong */
    double GetFitParameter(unsigned long spectrumNum, unsigned long specieIndex, FIT_PARAMETER parameter) const;

    /** returns true if the given index is a valid spectrum index */
    inline bool IsValidSpectrumIndex(unsigned long spectrumNum) const { return (spectrumNum >= 0 && spectrumNum < m_specNum); }
};
}