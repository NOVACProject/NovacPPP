#include <PPPLib/Configuration/DarkCorrectionConfiguration.h>

namespace Configuration
{
CDarkCorrectionConfiguration::CDarkCorrectionConfiguration(void)
{
    Clear();
}

CDarkCorrectionConfiguration::~CDarkCorrectionConfiguration(void)
{
    m_darkSettings.clear();
}

void CDarkCorrectionConfiguration::Clear()
{
    m_darkSettings.clear();
}

void CDarkCorrectionConfiguration::InsertDarkCurrentCorrectionSettings(const CDarkSettings& dSettings, const novac::CDateTime& validFrom, const novac::CDateTime& validTo)
{

    // make a copy of the settings
    DarkSettingWithTime newSetting;
    newSetting.setting = dSettings;
    newSetting.validFrom = validFrom;
    newSetting.validTo = validTo;

    m_darkSettings.push_back(newSetting);
}

int CDarkCorrectionConfiguration::GetDarkSettings(CDarkSettings& dSettings, const novac::CDateTime& time) const
{
    for (const auto& setting : m_darkSettings)
    {
        if ((setting.validFrom < time || setting.validFrom == time) && time < setting.validTo)
        {
            dSettings = setting.setting;
            return 0;
        }
    }

    // always return the default if no special setting can be found.
    DarkSettingWithTime defaultSetting;
    dSettings = defaultSetting.setting;
    return 0;
}

int CDarkCorrectionConfiguration::GetDarkSettings(size_t index, CDarkSettings& dSettings, novac::CDateTime& validFrom, novac::CDateTime& validTo) const
{
    if (index >= GetSettingsNum())
    {
        return 1;
    }

    if (m_darkSettings.size() == 0)
    {
        DarkSettingWithTime defaultSetting;
        dSettings = defaultSetting.setting;
        validFrom = novac::CDateTime::MinValue();
        validTo = novac::CDateTime::MaxValue();
        return 0;
    }

    dSettings = m_darkSettings[index].setting;
    validFrom = m_darkSettings[index].validFrom;
    validTo = m_darkSettings[index].validTo;

    return 0;
}

size_t CDarkCorrectionConfiguration::GetSettingsNum() const
{
    // There's always the default setting..
    return std::max(static_cast<size_t>(1), m_darkSettings.size());
}

}
