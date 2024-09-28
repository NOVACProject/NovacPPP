#pragma once

#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

#include <PPPLib/Configuration/InstrumentConfiguration.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/PostProcessingStatistics.h>

namespace Evaluation
{

/** Checks the supplied scan if it's good enough to bother evaluating.
    @returns false if the scan is too bad and should be ignored. Else return true. */
bool IsGoodEnoughToEvaluate(
    const novac::CScanFileHandler& scan,
    const novac::CFitWindow& window,
    const novac::SpectrometerModel& model,
    const Configuration::CInstrumentLocation& instrLocation,
    const Configuration::CUserConfiguration& userSettings,
    ReasonForScanRejection& reason,
    std::string& reasonMessage);

}