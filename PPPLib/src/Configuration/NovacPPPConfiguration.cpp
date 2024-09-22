#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Logging.h>

// The global configuration object. TODO: Remove as a global variable when no longer used by any other classes.
Configuration::CNovacPPPConfiguration g_setup;

using namespace novac;

namespace Configuration
{

const CInstrumentConfiguration* CNovacPPPConfiguration::GetInstrument(const std::string& serial) const
{
    novac::CString errorMessage;

    for (auto& instrument : m_instrument)
    {
        if (Equals(instrument.m_serial, serial))
        {
            return &instrument;
        }
    }

    // nothing found
    errorMessage.Format("Recieved spectrum from not-configured instrument %s. Cannot Evaluate!", serial.c_str());
    ShowMessage(errorMessage);

    return nullptr;
}

CInstrumentLocation CNovacPPPConfiguration::GetInstrumentLocation(const std::string& serial, const CDateTime& day) const
{
    CInstrumentLocation singleLocation;

    // First of all find the instrument 
    const CInstrumentConfiguration* instrumentConf = GetInstrument(serial);
    if (instrumentConf == nullptr)
    {
        novac::CString errorMessage;
        errorMessage.Format("Cannot find configuration for instrument with serial number '%s'", serial.c_str());
        throw PPPLib::NotFoundException(errorMessage.std_str());
    }

    // Next find the instrument location that is valid for this date
    const CLocationConfiguration& locationconf = instrumentConf->m_location;
    for (unsigned int k = 0; k < locationconf.GetLocationNum(); ++k)
    {
        locationconf.GetLocation(k, singleLocation);

        if (singleLocation.m_validFrom < day && (day < singleLocation.m_validTo || day == singleLocation.m_validTo))
        {
            return singleLocation;
        }
    }

    novac::CString errorMessage;
    errorMessage.Format("Recieved spectrum from instrument %s which is does not have a configured location on %04d.%02d.%02d. Cannot Evaluate!", serial.c_str(), day.year, day.month, day.day);
    throw PPPLib::NotFoundException(errorMessage.std_str());
}

novac::CFitWindow CNovacPPPConfiguration::GetFitWindow(
    const std::string& serial,
    int channel,
    const CDateTime& dateAndTime,
    const novac::CString* fitWindowName) const
{
    novac::CString windowName;
    if (fitWindowName != nullptr)
    {
        windowName.Format(*fitWindowName);
    }
    else
    {
        windowName.Format("");
    }

    // First of all find the instrument 
    const CInstrumentConfiguration* instrumentConf = GetInstrument(serial);
    if (instrumentConf == nullptr)
    {
        novac::CString errorMessage;
        errorMessage.Format("Cannot find configuration for instrument with serial number '%s'", serial.c_str());
        throw PPPLib::NotFoundException(errorMessage.std_str());
    }

    // Then find the evaluation fit-window that is valid for this date
    const Configuration::CEvaluationConfiguration& evalConf = instrumentConf->m_eval;
    for (int k = 0; k < evalConf.NumberOfFitWindows(); ++k)
    {
        novac::CFitWindow window;
        CDateTime evalValidFrom, evalValidTo;
        evalConf.GetFitWindow(k, window, evalValidFrom, evalValidTo);

        if (evalValidFrom < dateAndTime && (dateAndTime < evalValidTo || dateAndTime == evalValidTo) && ((channel % 16) == window.channel))
        {
            if (windowName.GetLength() >= 1)
            {
                if (Equals(windowName, window.name))
                {
                    // if we're searching for a specific name of fit-windows
                    //	then the name must also match
                    return window;
                }
            }
            else
            {
                // if we're not searching for any specific name, then anything
                //	which is valid within the given time-range will do.
                return window;
            }
        }
    }

    novac::CString errorMessage;
    if (windowName.GetLength() >= 1)
    {
        errorMessage.Format("Recieved spectrum from instrument %s which is does not have a configured fit-window \"%s\" on %04d.%02d.%02d. Cannot Evaluate!", serial.c_str(), (const char*)windowName, dateAndTime.year, dateAndTime.month, dateAndTime.day);
    }
    else
    {
        errorMessage.Format("Recieved spectrum from instrument %s which is does not have a configured fit-window on %04d.%02d.%02d. Cannot Evaluate!", serial.c_str(), dateAndTime.year, dateAndTime.month, dateAndTime.day);
    }
    throw PPPLib::NotFoundException(errorMessage.std_str());
}

CDarkSettings CNovacPPPConfiguration::GetDarkCorrection(const std::string& serial, const CDateTime& day) const
{
    // First of all find the instrument 
    const CInstrumentConfiguration* instrumentConf = GetInstrument(serial);
    if (instrumentConf == nullptr)
    {
        novac::CString errorMessage;
        errorMessage.Format("Cannot find configuration for instrument with serial number '%s'", serial.c_str());
        throw PPPLib::NotFoundException(errorMessage.std_str());
    }

    // Next find the CDarkCorrectionConfiguration that is valid for this date
    const CDarkCorrectionConfiguration& darkConf = instrumentConf->m_darkCurrentCorrection;
    CDarkSettings settings;
    if (0 == darkConf.GetDarkSettings(settings, day))
    {
        return settings;
    }

    novac::CString errorMessage;
    errorMessage.Format("Recieved spectrum from instrument %s which is does not have a configured dark current configuration on %04d.%02d.%02d. Cannot Evaluate!", serial.c_str(), day.year, day.month, day.day);
    throw PPPLib::NotFoundException(errorMessage.std_str());
}
}

