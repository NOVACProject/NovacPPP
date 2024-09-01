#pragma once

#include <map>
#include <string>
#include <vector>
#include <PPPLib/SpectrometerId.h>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Evaluation/ReferenceFile.h>

namespace novac
{
class CPostCalibrationStatistics
{
public:

    void RememberCalibrationPerformed(
        const SpectrometerId& instrument,
        const novac::CDateTime& gpsTimeOfScan,
        const std::vector<novac::CReferenceFile>& referencesCreated);

    int GetNumberOfCalibrationsPerformedFor(const SpectrometerId& instrument) const;

    /** Retrieves the properties of one calibration performed for the provided spectrometer.
        The calibrations will be retrieved in order increasing time, such that validFrom and validTo will always increase
            with increasing index.
        @param validFrom The time stamp midway between the retrieved calibration and the current. For this first, this is all zeroes.
        @param validTo The timestamp midway between the retrieved calibration and the next. For the last calibration, this is the year 9999.
        @throws std::invalid_argument if index < 0 or index >= GetNumberOfCalibrationsPerformedFor for the given instrument.*/
    void GetCalibration(
        const SpectrometerId& instrument,
        int index,
        novac::CDateTime& validFrom,
        novac::CDateTime& validTo,
        std::vector<novac::CReferenceFile>& referencesCreated) const;

private:
    struct InstrumentCalibration
    {
        // The time stamp of the calibrated scans. This is the time the scan started (utc).
        novac::CDateTime calibrationTimeStamp;

        std::vector<novac::CReferenceFile> referenceFiles;
    };

    struct InstrumentCalibrationCollection
    {
        // The list of calibrations performed for this particular instrument.
        // This list shall be sorted in increasing time such that the last scan is last in the list.
        std::vector<InstrumentCalibration> calibrationsPerformed;
    };

    // Map of instrument serial&channel number vs calibration properties
    std::map<SpectrometerId, InstrumentCalibrationCollection> m_calibrations;

};
}