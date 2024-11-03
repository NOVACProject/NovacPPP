#pragma once

#include <vector>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>

namespace Configuration
{
class CDarkCorrectionConfiguration
{
public:
    CDarkCorrectionConfiguration(void);
    ~CDarkCorrectionConfiguration(void);

    /** Clears all settings */
    void Clear();

    /** Inserts a new set of settings for correcting dark spectra into the configuration for this spectrometer.
        @param dSettings - dark-current correction settings
        @param validFrom - the time from which this settings is valid, NULL if valid since the beginning of time
        @param validTo - the time to which this settings is valid, NULL if valid until the end of time */
    void InsertDarkCurrentCorrectionSettings(const CDarkSettings& dSettings, const novac::CDateTime& validFrom, const novac::CDateTime& validTo);

    /** Retrieves how the dark-current should be corrected for this spectrometer at the
            given time.
        @param dSettings - will on successfull return be filled with the settings that are
            valid at the given time
        @param time - the time when we want to know the settings for how the
            dark current should be corrected
        @return 0 if sucessful, otherwise 1 */
    int GetDarkSettings(CDarkSettings& dSettings, const novac::CDateTime& time) const;

    /** Retrieves a dark-current settings from the configuration for this spectrometer.
        @param index - the index of the configuration to get. If this is < 0 or
            larger than the number of dark-current settings configured this function
            returns 1 and the dark-current setting is undefined
        @param dSettings - the dark-current settings to get
        @param validFrom - the time from which this fit-window is valid
        @param validTo - the time to which this fit-window is valid
        @return 0 if sucessful, otherwise 1 */
    int GetDarkSettings(size_t index, CDarkSettings& dSettings, novac::CDateTime& validFrom, novac::CDateTime& validTo) const;

    /** Gets the number of dark-current settings configured for this spectrometer */
    size_t GetSettingsNum() const;


private:

    /** Basic structure for keeping track of a fit window and the time-range for which it is valid */
    struct DarkSettingWithTime
    {
        CDarkSettings setting;
        novac::CDateTime validFrom;
        novac::CDateTime validTo;
    };

    /** The array of dark-current correction options for this spectrometer */
    std::vector<DarkSettingWithTime> m_darkSettings;

};
}