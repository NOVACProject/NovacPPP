#pragma once

#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

#include <PPPLib/Configuration/InstrumentConfiguration.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/PostProcessingStatistics.h>
#include <PPPLib/Measurement.h>

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


/** Gets the filename under which the scan-file should be stored.
    @return true if a filename is found. */
bool GetArchivingfileName(
    novac::CString& pakFile,
    novac::CString& txtFile,
    const novac::CString& fitWindowName,
    const novac::CString& temporaryScanFile,
    std::string outputDirectory,
    MEASUREMENT_MODE mode);

}