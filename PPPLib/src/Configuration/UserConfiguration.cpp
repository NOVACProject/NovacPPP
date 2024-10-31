#include <PPPLib/Configuration/UserConfiguration.h>
#include <algorithm>

// The global user-configuration object
Configuration::CUserConfiguration g_userSettings;

namespace Configuration
{
CUserConfiguration::CUserConfiguration()
{
    Clear();
}

void CUserConfiguration::Clear()
{
    m_maxThreadNum = 2;

    m_fIsContinuation = false;

    // don't know where to put the output data
    m_tempDirectory.Format("");
    m_outputDirectory.Format("");

    m_processingMode = PROCESSING_MODE::PROCESSING_MODE_FLUX;
    m_doEvaluations = true;

    // the molecule
    m_molecule = novac::StandardMolecule::SO2;

    // The volcano that we should process
    m_volcano = 0;

    // The time range where to search for data
    m_fromDate = novac::CDateTime(2005, 10, 01, 00, 00, 00);
    m_toDate.SetToNow();

    // The directory where to search for data
    m_LocalDirectory.Format("C:\\Novac\\Data\\");
    m_includeSubDirectories_Local = 1;

    // The FTP-directory where to search for data
    m_includeSubDirectories_FTP = 1;

    // If we should upload the results to the ftp or not
    m_uploadResults = 0;

    // the settings for the fit-windows to use
    m_nFitWindowsToUse = 1;
    m_fitWindowsToUse[0].Format("SO2");
    for (int k = 1; k < MAX_FIT_WINDOWS; ++k)
    {
        m_fitWindowsToUse[k].Format("");
    }
    m_mainFitWindow = 0;

    // The settings for the sky-spectrum
    sky.skyOption = Configuration::SKY_OPTION::MEASURED_IN_SCAN;
    sky.indexInScan = 0;
    sky.skySpectrumFile = "";

    // the wind field
    m_windFieldFile.Format("");
    m_windFieldFileOption = 0;

    // The geometry calculations
    m_calcGeometry_CompletenessLimit = 0.7;
    m_calcGeometryValidTime = 10 * 60;
    m_calcGeometry_MaxTimeDifference = 900;
    m_calcGeometry_MinDistance = 200;
    m_calcGeometry_MaxDistance = 10000;
    m_calcGeometry_MaxPlumeAltError = 500.0;
    m_calcGeometry_MaxWindDirectionError = 10.0;

    // the dual-beam calculations
    m_fUseMaxTestLength_DualBeam = true;
    m_dualBeam_MaxWindSpeedError = 10.0;
    m_dualBeam_ValidTime = 15 * 60;

    // The quality parameters
    m_completenessLimitFlux = 0.9;
    m_minimumSaturationInFitRegion = 0.05;
    m_maxExposureTime_got = 900;
    m_maxExposureTime_hei = 4000;

}

/** Compares to set of settings */
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
    if (!Equals(m_LocalDirectory, settings2.m_LocalDirectory))
        return false;
    if (settings2.m_includeSubDirectories_Local != m_includeSubDirectories_Local)
        return false;

    // The FTP-directory where to search for data
    if (!Equals(m_FTPDirectory, settings2.m_FTPDirectory))
        return false;
    if (!Equals(m_FTPUsername, settings2.m_FTPUsername))
        return false;
    if (!Equals(m_FTPPassword, settings2.m_FTPPassword))
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
