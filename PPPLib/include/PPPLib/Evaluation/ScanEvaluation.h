#pragma once

#include <PPPLib/Evaluation/ScanResult.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/ScanEvaluationBase.h>
#include <SpectralEvaluation/Log.h>

namespace novac
{
class CFitWindow;
class CEvaluationBase;
}

namespace Evaluation
{

/**
    An object of the <b>CScanEvaluation</b>-class handles the evaluation of one
    scan. The actual evaluation is performed by calling the function
    <b>EvaluateScan</b> with the name of the .pak-file that we want to evaluate,
    an evaluator to use for the evaluation and parameters for the dark correction.
*/
class CScanEvaluation : public novac::ScanEvaluationBase
{

public:
    CScanEvaluation(const Configuration::CUserConfiguration& userSettings, novac::ILogger& log);

    ~CScanEvaluation();

    /** Called to evaluate the one scan.
    *   @param scan an opened scan file, containing one single scan.
    *   @param fitWindow The settings for the fit to use.
    *   @param spectrometerModel The model of the spectrometer, needed to verify intensities and device specific settings.
    *   @param darkSettings An optional set of settings for the dark correction. If null then the default settings will be used.
        @return the result of the evaluation. This is null if the evaluation failed.
        @throws std::exception (or subclass of this) if the given fit window is not ok */
    std::unique_ptr<CScanResult> EvaluateScan(
        novac::CScanFileHandler& scan,
        const novac::CFitWindow& fitWindow,
        const novac::SpectrometerModel& spectrometerModel,
        const Configuration::CDarkSettings* darkSettings = nullptr);

private:

    const Configuration::CUserConfiguration& m_userSettings;

    novac::ILogger& m_log;

    // ----------------------- PRIVATE METHODS ---------------------------

    /** Performs the evaluation using the supplied evaluator
        @return the result of the evaluation, or null if something goes wrong */
    std::unique_ptr<CScanResult> EvaluateOpenedScan(
        novac::LogContext logContext,
        novac::CScanFileHandler& scan,
        std::unique_ptr<novac::CEvaluationBase>& eval,
        const novac::SpectrometerModel& spectrometer,
        const Configuration::CDarkSettings* darkSettings = nullptr);

    /** This returns the sky spectrum that is to be used in the fitting.
        Which spectrum to be used is taken from the given settings.
        @return true on success. */
    bool GetSky(novac::CScanFileHandler& scan, const Configuration::CSkySettings& settings, novac::CSpectrum& sky);

    /** This returns the dark spectrum that is to be used in the fitting.
        @param scan - the scan-file handler from which to get the dark spectrum
        @param spec - the spectrum for which the dark spectrum should be retrieved
        @param dark - will on return be filled with the dark spectrum
        @param darkSettings - the settings for how to get the dark spectrum from this spectrometer.
        @return true on success. */
    bool GetDark(novac::CScanFileHandler& scan, const novac::CSpectrum& spec, novac::CSpectrum& dark, const Configuration::CDarkSettings* darkSettings = NULL);

    /** Finds the optimum shift and squeeze for an evaluated scan by looking at
            the spectrum with the highest absorption of the evaluated specie
            and evaluate it with shift and squeeze free
        @param fitWindow The old fit-window where we should try to improve the settings.
        @param index indexOfMostAbsorbingSpectrum The index of the spectrum which has the highest absorption,
            this is the spectrum to evaluate again.
        @param scan The scan to improve the evaluation for.
        @return a new evaluator with the fit-window set to the new optimum values.
        @return nullptr if the evaluation failed. */
    novac::CEvaluationBase* FindOptimumShiftAndSqueeze(
        novac::LogContext logContext,
        const novac::CFitWindow& fitWindow,
        int indexOfMostAbsorbingSpectrum,
        novac::CScanFileHandler& scan,
        const Configuration::CDarkSettings* darkSettings = nullptr);

    // ------------------------ THE PARAMETERS FOR THE EVALUATION ------------------

    /** Remember the index of the spectrum with the highest absorption, to be able to
        adjust the shift and squeeze with it later */
    int m_indexOfMostAbsorbingSpectrum = -1;

    /** Performs a basic validation on the setup of the given fit window.
        @throws std::exception (or subclass of this) if the window is not ok */
    static void ValidateSetup(const novac::CFitWindow& window);
};
}