#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <PPPLib/Logging.h>

// The global configuration object
Configuration::CNovacPPPConfiguration g_setup;

using namespace novac;

namespace Configuration
{
CNovacPPPConfiguration::CNovacPPPConfiguration()
{
}

const CInstrumentConfiguration* CNovacPPPConfiguration::GetInstrument(const novac::CString& serial) const
{
    std::string stdSerial{ (const char*)serial };
    return GetInstrument(stdSerial);
}

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

int CNovacPPPConfiguration::GetInstrumentLocation(const novac::CString& serial, const CDateTime& day, CInstrumentLocation& instrLocation) const
{
    CInstrumentLocation singleLocation;

    // First of all find the instrument 
    const CInstrumentConfiguration* instrumentConf = GetInstrument(serial);
    if (instrumentConf == nullptr)
        return 1;

    // Next find the instrument location that is valid for this date
    const CLocationConfiguration& locationconf = instrumentConf->m_location;
    bool foundValidLocation = false;
    for (unsigned int k = 0; k < locationconf.GetLocationNum(); ++k)
    {
        locationconf.GetLocation(k, singleLocation);

        if (singleLocation.m_validFrom < day && (day < singleLocation.m_validTo || day == singleLocation.m_validTo))
        {
            instrLocation = singleLocation;
            foundValidLocation = true;
            break;
        }
    }

    if (!foundValidLocation)
    {
        novac::CString errorMessage;
        errorMessage.Format("Recieved spectrum from instrument %s which is does not have a configured location on %04d.%02d.%02d. Cannot Evaluate!", (const char*)serial, day.year, day.month, day.day);
        ShowMessage(errorMessage);
        return 1;
    }

    return 0;
}

int CNovacPPPConfiguration::GetFitWindow(
    const novac::CString& serial,
    int channel,
    const CDateTime& dateAndTime,
    novac::CFitWindow& window,
    const novac::CString* fitWindowName) const
{

    CDateTime evalValidFrom, evalValidTo;
    novac::CString errorMessage, windowName;

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
        return 1;

    // Then find the evaluation fit-window that is valid for this date
    const Configuration::CEvaluationConfiguration& evalConf = instrumentConf->m_eval;
    bool foundValidEvaluation = false;
    for (int k = 0; k < evalConf.NumberOfFitWindows(); ++k)
    {
        evalConf.GetFitWindow(k, window, evalValidFrom, evalValidTo);

        if (evalValidFrom < dateAndTime && (dateAndTime < evalValidTo || dateAndTime == evalValidTo) && ((channel % 16) == window.channel))
        {
            if (windowName.GetLength() >= 1)
            {
                if (Equals(windowName, window.name))
                {
                    // if we're searching for a specific name of fit-windows
                    //	then the name must also match
                    foundValidEvaluation = true;
                    break;
                }
            }
            else
            {
                // if we're not searching for any specific name, then anything
                //	which is valid within the given time-range will do.
                foundValidEvaluation = true;
                break;
            }
        }
    }

    if (!foundValidEvaluation)
    {
        if (windowName.GetLength() >= 1)
        {
            errorMessage.Format("Recieved spectrum from instrument %s which is does not have a configured fit-window \"%s\" on %04d.%02d.%02d. Cannot Evaluate!", (const char*)serial, (const char*)windowName, dateAndTime.year, dateAndTime.month, dateAndTime.day);
        }
        else
        {
            errorMessage.Format("Recieved spectrum from instrument %s which is does not have a configured fit-window on %04d.%02d.%02d. Cannot Evaluate!", (const char*)serial, dateAndTime.year, dateAndTime.month, dateAndTime.day);
        }
        ShowMessage(errorMessage);
        return 1;
    }

    return 0;
}

int CNovacPPPConfiguration::GetDarkCorrection(const novac::CString& serial, const CDateTime& dateAndTime, CDarkSettings& settings) const
{
    // First of all find the instrument 
    const CInstrumentConfiguration* instrumentConf = GetInstrument(serial);
    if (instrumentConf == nullptr)
        return 1;

    // Next find the CDarkCorrectionConfiguration that is valid for this date
    const CDarkCorrectionConfiguration& darkConf = instrumentConf->m_darkCurrentCorrection;
    return darkConf.GetDarkSettings(settings, dateAndTime);
}
}

