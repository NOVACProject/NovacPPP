#include <PPPLib/Configuration/UserConfiguration.h>
#include <algorithm>

namespace Configuration
{

CUserConfiguration::CUserConfiguration()
{
    sky.skyOption = Configuration::SKY_OPTION::MEASURED_IN_SCAN;
    sky.indexInScan = 0;
    sky.skySpectrumFile = "";
}

bool CUserConfiguration::operator==(const CUserConfiguration& settings2)
{

    if (!Equals(m_tempDirectory, settings2.m_tempDirectory))
        return false;

    if (!Equals(m_outputDirectory, settings2.m_outputDirectory))
        return false;

    // the processing mode
    if (m_processingMode != settings2.m_processingMode)
        return false;

    // the molecule
    if (m_molecule != settings2.m_molecule)
        return false;

    // The volcano that we should process
    if (settings2.m_volcano != m_volcano)
        return false;

    // The time range where to search for data
    if (settings2.m_fromDate.year != this->m_fromDate.year)
        return false;
    if (settings2.m_fromDate.month != this->m_fromDate.month)
        return false;
    if (settings2.m_fromDate.day != this->m_fromDate.day)
        return false;
    if (settings2.m_toDate.year != this->m_toDate.year)
        return false;
    if (settings2.m_toDate.month != this->m_toDate.month)
        return false;
    if (settings2.m_toDate.day != this->m_toDate.day)
        return false;

    // The directory where to search for data
    if (!novac::EqualsIgnoringCase(m_LocalDirectory, settings2.m_LocalDirectory))
        return false;
    if (settings2.m_includeSubDirectories_Local != m_includeSubDirectories_Local)
        return false;

    // The FTP-directory where to search for data
    if (!novac::EqualsIgnoringCase(m_FTPDirectory, settings2.m_FTPDirectory))
        return false;
    if (!novac::EqualsIgnoringCase(m_FTPUsername, settings2.m_FTPUsername))
        return false;
    if (!novac::EqualsIgnoringCase(m_FTPPassword, settings2.m_FTPPassword))
        return false;
    if (settings2.m_includeSubDirectories_FTP != m_includeSubDirectories_FTP)
        return false;

    // uploading of results?
    if (settings2.m_uploadResults != m_uploadResults)
        return false;

    // the settings for the fit-windows to use
    if (settings2.m_nFitWindowsToUse != m_nFitWindowsToUse)
        return false;
    if (settings2.m_mainFitWindow != m_mainFitWindow)
        return false;

    for (int k = 0; k < MAX_FIT_WINDOWS; ++k)
    {
        if (!Equals(m_fitWindowsToUse[k], settings2.m_fitWindowsToUse[k]))
            return false;
    }

    // The settings for the sky-spectrum
    if (settings2.sky.skyOption != sky.skyOption)
        return false;
    if (settings2.sky.indexInScan != sky.indexInScan)
        return false;
    if (sky.skySpectrumFile.compare(settings2.sky.skySpectrumFile) != 0)
        return false;

    // the wind field
    if (!Equals(m_windFieldFile, settings2.m_windFieldFile))
        return false;
    if (m_windFieldFileOption != settings2.m_windFieldFileOption)
        return false;

    // The geometry calculations
    if (std::abs(settings2.m_calcGeometry_CompletenessLimit - m_calcGeometry_CompletenessLimit) > 0.01)
        return false;
    if (settings2.m_calcGeometryValidTime != m_calcGeometryValidTime)
        return false;
    if (settings2.m_calcGeometry_MaxTimeDifference != m_calcGeometry_MaxTimeDifference)
        return false;
    if (settings2.m_calcGeometry_MinDistance != m_calcGeometry_MinDistance)
        return false;
    if (settings2.m_calcGeometry_MaxDistance != m_calcGeometry_MaxDistance)
        return false;
    if (settings2.m_calcGeometry_MaxPlumeAltError != m_calcGeometry_MaxPlumeAltError)
        return false;
    if (settings2.m_calcGeometry_MaxWindDirectionError != m_calcGeometry_MaxWindDirectionError)
        return false;

    // the dual-beam calculations
    if (settings2.m_fUseMaxTestLength_DualBeam != m_fUseMaxTestLength_DualBeam)
        return false;
    if (settings2.m_dualBeam_MaxWindSpeedError != m_dualBeam_MaxWindSpeedError)
        return false;
    if (settings2.m_dualBeam_ValidTime != m_dualBeam_ValidTime)
        return false;

    // The quality parameters
    if (std::abs(settings2.m_completenessLimitFlux - m_completenessLimitFlux) > 0.01)
        return false;
    if (std::abs(settings2.m_minimumSaturationInFitRegion - m_minimumSaturationInFitRegion) > 0.01)
        return false;
    if (settings2.m_maxExposureTime_got != m_maxExposureTime_got)
        return false;
    if (settings2.m_maxExposureTime_hei != m_maxExposureTime_hei)
        return false;

    // we've checked every single parameters, these two settings are identical
    return true;
}
}
