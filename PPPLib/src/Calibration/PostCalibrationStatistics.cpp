#include <PPPLib/Calibration/PostCalibrationStatistics.h>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <assert.h>
#include <cmath>

using namespace novac;

void CPostCalibrationStatistics::RememberCalibrationPerformed(
    const SpectrometerId& instrument,
    const novac::CDateTime& gpsTimeOfScan,
    const std::vector<novac::CReferenceFile>& referencesCreated)
{
    InstrumentCalibration newCalibration;
    newCalibration.calibrationTimeStamp = gpsTimeOfScan;
    newCalibration.referenceFiles = referencesCreated;

    auto pos = m_calibrations.find(instrument);

    if (pos == m_calibrations.end())
    {
        InstrumentCalibrationCollection newCollection;
        newCollection.calibrationsPerformed.push_back(newCalibration);
        m_calibrations[instrument] = newCollection;
    }
    else
    {
        InstrumentCalibrationCollection& collection = pos->second;
        collection.calibrationsPerformed.push_back(newCalibration);

        // keep the collection sorted
        std::sort(
            begin(collection.calibrationsPerformed),
            end(collection.calibrationsPerformed),
            [](const InstrumentCalibration& c1, const InstrumentCalibration& c2)
            {
                return c1.calibrationTimeStamp < c2.calibrationTimeStamp;
            });
    }
}

int CPostCalibrationStatistics::GetNumberOfCalibrationsPerformedFor(const SpectrometerId& instrument) const
{
    auto pos = m_calibrations.find(instrument);
    if (pos == m_calibrations.end())
    {
        return 0;
    }

    return static_cast<int>(pos->second.calibrationsPerformed.size());
}

void CPostCalibrationStatistics::GetCalibration(
    const SpectrometerId& instrument,
    int index,
    novac::CDateTime& validFrom,
    novac::CDateTime& validTo,
    std::vector<novac::CReferenceFile>& referencesCreated) const
{
    auto pos = m_calibrations.find(instrument);
    if (pos == m_calibrations.end())
    {
        std::stringstream message;
        message << "No calibrations have been performed for device (serial: " << instrument.serial << ", channel: " << instrument.channel << ")";
        return throw std::invalid_argument(message.str());
    }
    if (index < 0 || index >= static_cast<int>(pos->second.calibrationsPerformed.size()))
    {
        std::stringstream message;
        message << "Invalid index '" << index << "' when retrieving calibrations for device (serial: " << instrument.serial << ", channel: " << instrument.channel << ")";
        return throw std::invalid_argument(message.str());
    }

    const size_t idx = static_cast<size_t>(index);
    if (idx == 0)
    {
        validFrom = novac::CDateTime(0, 0, 0, 0, 0, 0);
    }
    else
    {
        // Calculate the midpoint between previous and current
        validFrom = pos->second.calibrationsPerformed[idx - 1].calibrationTimeStamp;
        const auto& current = pos->second.calibrationsPerformed[idx].calibrationTimeStamp;
        const double secondsDifference = novac::CDateTime::Difference(current, validFrom);
        if (secondsDifference < 0.0)
        {
            throw std::logic_error("Incorrect logic discovered in PostCalibrationStatistics, the timestamps were not sorted in increasing order.");
        }

        validFrom.Increment(static_cast<int>(std::round(0.5 * secondsDifference)));
    }

    if (idx < pos->second.calibrationsPerformed.size() - 1)
    {
        // Calculate the midpoint between current and next
        validTo = pos->second.calibrationsPerformed[idx].calibrationTimeStamp;
        const auto& next = pos->second.calibrationsPerformed[idx + 1].calibrationTimeStamp;
        const double secondsDifference = novac::CDateTime::Difference(next, validTo);
        if (secondsDifference < 0.0)
        {
            throw std::logic_error("Incorrect logic discovered in PostCalibrationStatistics, the timestamps were not sorted in increasing order.");
        }

        validTo.Increment(static_cast<int>(std::round(0.5 * secondsDifference)));
    }
    else
    {
        validTo = novac::CDateTime(9999, 12, 31, 23, 59, 59);
    }

    referencesCreated = pos->second.calibrationsPerformed[idx].referenceFiles;
}