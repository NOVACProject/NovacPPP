#pragma once

#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>

namespace Configuration
{
/** Basic structure for keeping track of a fit window and the time-range for which it is valid */
struct FitWindowWithTime
{
    novac::CFitWindow window;
    novac::CDateTime validFrom;
    novac::CDateTime validTo;
};

/** The class CEvaluationConfiguration helps keeping track of all the defined fit-windows for one spectrometer
    (notice that there may be multiple channels on this spectrometer however, some fit windows may be for one channel and
    some may be for the other). */
class CEvaluationConfiguration
{
public:
    CEvaluationConfiguration();
    ~CEvaluationConfiguration();

    /** The serial-number of the spectrometer for which this setup is valid */
    std::string m_serial;

    /** Clears all configurations for this spectrometer */
    void Clear();

    /** Inserts a new fit-window into the configuration for this spectrometer.
        @param window - the fit-window to insert
        @param validFrom - the time from which this fit-window is valid.
        @param validTo - the time to which this fit-window is valid.  */
    void InsertFitWindow(const novac::CFitWindow& window, const novac::CDateTime& validFrom, const novac::CDateTime& validTo);

    /** Sets the properties of the fit-window number 'index'
        @param index - the index of the configuration to set. If this is < 0
            this function returns 1 and nothing is changed
        @param window - the fit-window to set
        @param validFrom - the time from which this fit-window is valid.
        @param validTo - the time to which this fit-window is valid.
        @return 0 if sucessful, otherwise 1 */
    int SetFitWindow(int index, const novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo);

    /** Retrieves a fit-window from the configuration for this spectrometer.
        @param index - the index of the configuration to get. If this is < 0 or
            larger than the number of fit-windows configured this function
            returns 1 and the window is undefined
        @param window - the fit-window to get
        @param validFrom - the time from which this fit-window is valid
        @param validTo - the time to which this fit-window is valid
        @return 0 if sucessful, otherwise 1 */
    int GetFitWindow(int index, novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo) const;

    /** Gets the number of fit-windows configured for this spectrometer */
    int NumberOfFitWindows() const { return static_cast<int>(m_windows.size()); }

    /** This goes through the settings for the fit-windows to test that the settings
        make sense.
        @throw std::invalid_argument with an error message if the settings are not valid */
    void CheckSettings() const;

private:

    /** The array of fit-windows for this spectrometer */
    std::vector<FitWindowWithTime> m_windows;
};
}