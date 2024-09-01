#include "FluxCalculator.h"

// Here we need to know about the global settings
#include <PPPLib/Configuration/NovacPPPConfiguration.h>

// ... support for handling the evaluation-log files...
#include "../Common/EvaluationLogFileHandler.h"

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

#include <PPPLib/File/Filesystem.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <Poco/Path.h>

extern Configuration::CNovacPPPConfiguration g_setup;	   // <-- The settings
extern Configuration::CUserConfiguration g_userSettings;// <-- The settings of the user

using namespace Flux;
using namespace novac;


CFluxCalculator::CFluxCalculator(void)
{
}

CFluxCalculator::~CFluxCalculator(void)
{
}

/** Calculates the flux from the scan found in the given evaluation log file
    @param evalLogFileName - the name of the .txt-file that contains
        the result of the evaluation.
    @return 0 on success, else non-zero value
    */
int CFluxCalculator::CalculateFlux(const novac::CString& evalLogFileName, const Meteorology::CWindDataBase& windDataBase, const Geometry::CPlumeHeight& plumeAltitude, CFluxResult& fluxResult)
{
    CDateTime skyStartTime;
    novac::CString errorMessage, shortFileName, serial;
    Geometry::CPlumeHeight relativePlumeHeight;
    Meteorology::CWindField windField;
    int channel;
    MEASUREMENT_MODE mode;

    // 1. Assert that the evaluation-log-file exists
    if (!Filesystem::IsExistingFile(evalLogFileName))
    {
        ShowMessage("Recieved evaluation-log with illegal path. Could not calculate flux.");
        return 1;
    }

    // 2. Get some information about the scan from the file-name
    shortFileName.Format(evalLogFileName);
    Common::GetFileName(shortFileName);
    novac::CFileUtils::GetInfoFromFileName(shortFileName, skyStartTime, serial, channel, mode);

    // 3. Find the location of this instrument
    Configuration::CInstrumentLocation instrLocation;
    if (GetLocation(serial, skyStartTime, instrLocation))
    {
        errorMessage.Format("Failed to find location of instrument %s on %04d.%02d.%02d. Could not calculate flux", (const char*)serial, skyStartTime.year, skyStartTime.month, skyStartTime.day);
        ShowMessage(errorMessage);
        return 3;
    }

    // 4. Get the wind field at the time of the collection of this scan
    if (!windDataBase.GetWindField(skyStartTime, CGPSData(instrLocation.m_latitude, instrLocation.m_longitude, plumeAltitude.m_plumeAltitude), Meteorology::INTERP_NEAREST_NEIGHBOUR, windField))
    {
        ShowMessage("Could not retrieve wind field at time of measurement. Could not calculate flux");
        return 4;
    }

    // 4b. Adjust the altitude of the plume so that it is relative to the altitude of the instrument...
    relativePlumeHeight = plumeAltitude;
    relativePlumeHeight.m_plumeAltitude -= instrLocation.m_altitude;
    if (relativePlumeHeight.m_plumeAltitude <= 0)
    {
        ShowMessage(" Negative plume height obtained when calculating flux. No flux can be calculated.");
        return 6;
    }

    // 5. Read in the evaluation log file 
    FileHandler::CEvaluationLogFileHandler reader;
    reader.m_evaluationLog.Format(evalLogFileName);
    reader.ReadEvaluationLog();
    if (reader.m_scan.size() == 0)
    {
        errorMessage.Format("Recieved evaluation log file (%s) with no scans inside. Cannot calculate flux", (const char*)evalLogFileName);
        ShowMessage(errorMessage);
        return 2;
    }
    else if (reader.m_scan.size() > 1)
    {
        errorMessage.Format("Recieved evaluation log file (%s) with more than one scans inside. Can only calculate flux for the first scan.", (const char*)evalLogFileName);
        ShowMessage(errorMessage);
    }

    // 6. extract the scan
    Evaluation::CScanResult& result = reader.m_scan[0];

    // 6b. Improve on the start-time of the scan...
    result.GetSkyStartTime(skyStartTime);

    // 7. Calculate the offset of the scan
    if (result.CalculateOffset(CMolecule(g_userSettings.m_molecule)))
    {
        ShowMessage("Could not calculate offset for scan. No flux can be calculated.");
        return 7;
    }

    // 8. Check that the completeness is higher than our limit...
    if (!result.CalculatePlumeCentre(CMolecule(g_userSettings.m_molecule)))
    {
        ShowMessage(" - Scan does not see the plume, no flux can be calculated");
        return 8;
    }
    double completeness = result.GetCalculatedPlumeCompleteness();
    if (completeness < (g_userSettings.m_completenessLimitFlux + 0.01))
    {
        errorMessage.Format(" - Scan has completeness = %.2lf which is less than limit of %.2lf. Rejected!", completeness, g_userSettings.m_completenessLimitFlux);
        ShowMessage(errorMessage);
        return 8;
    }

    // 9. Calculate the flux
    if (result.CalculateFlux(CMolecule(g_userSettings.m_molecule), windField, relativePlumeHeight, instrLocation.m_compass, instrLocation.m_coneangle, instrLocation.m_tilt))
    {
        ShowMessage("Flux calculation failed. No flux generated");
        return 9;
    }
    fluxResult = result.GetFluxResult();

    // 10. make a simple estimate of the quality of the flux measurement
    FLUX_QUALITY_FLAG windFieldQuality = FLUX_QUALITY_GREEN;
    FLUX_QUALITY_FLAG plumeHeightQuality = FLUX_QUALITY_GREEN;
    FLUX_QUALITY_FLAG completenessQuality = FLUX_QUALITY_GREEN;

    switch (windField.GetWindSpeedSource())
    {
    case Meteorology::MET_DEFAULT:				windFieldQuality = FLUX_QUALITY_RED; break;
    case Meteorology::MET_USER:					windFieldQuality = FLUX_QUALITY_RED; break;
    case Meteorology::MET_ECMWF_FORECAST:		windFieldQuality = FLUX_QUALITY_GREEN; break;
    case Meteorology::MET_ECMWF_ANALYSIS:		windFieldQuality = FLUX_QUALITY_GREEN; break;
    case Meteorology::MET_DUAL_BEAM_MEASUREMENT:windFieldQuality = FLUX_QUALITY_GREEN; break;
    case Meteorology::MET_MODEL_WRF:			windFieldQuality = FLUX_QUALITY_GREEN; break;
    case Meteorology::MET_NOAA_GDAS:			windFieldQuality = FLUX_QUALITY_GREEN; break;
    case Meteorology::MET_NOAA_FNL:				windFieldQuality = FLUX_QUALITY_GREEN; break;
    default:									windFieldQuality = FLUX_QUALITY_YELLOW; break;
    }
    switch (relativePlumeHeight.m_plumeAltitudeSource)
    {
    case Meteorology::MET_DEFAULT:				plumeHeightQuality = FLUX_QUALITY_RED; break;
    case Meteorology::MET_USER:					plumeHeightQuality = FLUX_QUALITY_RED; break;
    case Meteorology::MET_GEOMETRY_CALCULATION:	plumeHeightQuality = FLUX_QUALITY_GREEN; break;
    default:									plumeHeightQuality = FLUX_QUALITY_YELLOW; break;
    }
    if (fluxResult.m_completeness < 0.7)
    {
        completenessQuality = FLUX_QUALITY_RED;
    }
    else if (fluxResult.m_completeness < 0.9)
    {
        completenessQuality = FLUX_QUALITY_YELLOW;
    }
    else
    {
        completenessQuality = FLUX_QUALITY_GREEN;
    }

    // if any of the parameters has a low quality, then the result is judged to have a low quality...
    if (windFieldQuality == FLUX_QUALITY_RED || plumeHeightQuality == FLUX_QUALITY_RED || completenessQuality == FLUX_QUALITY_RED)
    {
        fluxResult.m_fluxQualityFlag = FLUX_QUALITY_RED;
    }
    else if (windFieldQuality == FLUX_QUALITY_YELLOW || plumeHeightQuality == FLUX_QUALITY_YELLOW || completenessQuality == FLUX_QUALITY_YELLOW)
    {
        fluxResult.m_fluxQualityFlag = FLUX_QUALITY_YELLOW;
    }
    else
    {
        fluxResult.m_fluxQualityFlag = FLUX_QUALITY_GREEN;
    }

    // ok
    return 0;
}

/** Looks in the configuration of the instruments and searches
    for a configured location which is valid for the spectrometer that
    collected the given scan and is also valid at the time when the scan
    was made.
    @return 0 if successful otherwise non-zero
*/
int CFluxCalculator::GetLocation(const novac::CString& serial, const CDateTime& startTime, Configuration::CInstrumentLocation& instrLocation)
{
    CDateTime day, evalValidFrom, evalValidTo;
    Configuration::CInstrumentConfiguration* instrumentConf = nullptr;
    Configuration::CInstrumentLocation singleLocation;
    novac::CString errorMessage;

    // First of all find the instrument 
    for (int k = 0; k < g_setup.NumberOfInstruments(); ++k)
    {
        if (Equals(g_setup.m_instrument[k].m_serial, serial))
        {
            instrumentConf = &g_setup.m_instrument[k];
            break;
        }
    }
    if (instrumentConf == nullptr)
    {
        errorMessage.Format("Recieved spectrum from not-configured instrument %s. Cannot calculate flux!", (const char*)serial);
        ShowMessage(errorMessage);
        return 1;
    }

    // Next find the instrument location that is valid for this date
    Configuration::CLocationConfiguration& locationconf = instrumentConf->m_location;
    bool foundValidLocation = false;
    for (int k = 0; k < (int)locationconf.GetLocationNum(); ++k)
    {
        locationconf.GetLocation(k, singleLocation);

        if (singleLocation.m_validFrom < startTime && (startTime < singleLocation.m_validTo || startTime == singleLocation.m_validTo))
        {
            instrLocation = singleLocation;
            foundValidLocation = true;
            break;
        }
    }
    if (!foundValidLocation)
    {
        errorMessage.Format("Recieved spectrum from instrument %s with not configured location. Cannot calculate flux!", (const char*)serial);
        ShowMessage(errorMessage);
        return 1;
    }

    // we're done! Return!
    return 0;
}

/** Appends the evaluated flux to the appropriate log file.
    @param scan - the scan itself, also containing information about the evaluation and the flux.
    @return SUCCESS if operation completed sucessfully. */
RETURN_CODE CFluxCalculator::WriteFluxResult(const Flux::CFluxResult& fluxResult, const Evaluation::CScanResult* result)
{
    novac::CString string, dateStr, dateStr2, serialNumber;
    novac::CString fluxLogFile, directory;
    novac::CString wdSrc, wsSrc, phSrc;
    novac::CString errorMessage;
    CDateTime dateTime;

    // 0. Get the sources for the wind-field
    fluxResult.m_windField.GetWindSpeedSource(wsSrc);
    fluxResult.m_windField.GetWindDirectionSource(wdSrc);
    Meteorology::MetSourceToString(fluxResult.m_plumeHeight.m_plumeAltitudeSource, phSrc);

    // 1. Output the day and time the scan that generated this measurement started
    result->GetSkyStartTime(dateTime);
    dateStr.Format("%04d-%02d-%02d", dateTime.year, dateTime.month, dateTime.day);
    dateStr2.Format("%04d.%02d.%02d", dateTime.year, dateTime.month, dateTime.day);
    string.Format("%s\t", (const char*)dateStr);
    string.AppendFormat("%02d:%02d:%02d\t", dateTime.hour, dateTime.minute, dateTime.second);

    // 2. Output the time the scan stopped
    result->GetStopTime(result->GetEvaluatedNum() - 1, dateTime);
    string.AppendFormat("%02d:%02d:%02d\t", dateTime.hour, dateTime.minute, dateTime.second);

    // 3. Output the flux result
    string.AppendFormat("%.2lf\t", fluxResult.m_flux);

    // 4. Output the wind speed (and the estimated error in the wind-speed)
    //	that was used for calculating this flux
    string.AppendFormat("%.3lf\t", fluxResult.m_windField.GetWindSpeed());
    string.AppendFormat("%.3lf\t", fluxResult.m_windField.GetWindSpeedError());

    // 5. Output the wind direction (and it's estimated error) 
    //	that was used for calculating this flux
    string.AppendFormat("%.3lf\t", fluxResult.m_windField.GetWindDirection());
    string.AppendFormat("%.3lf\t", fluxResult.m_windField.GetWindDirectionError());

    // 6. Output where we got the wind speed from
    string.AppendFormat("%s\t", (const char*)wsSrc);

    // 7. Output where we got the wind direction from
    string.AppendFormat("%s\t", (const char*)wdSrc);

    // 8. Output the plume height that was used for calculating this flux
    string.AppendFormat("%.1lf\t", fluxResult.m_plumeHeight.m_plumeAltitude);
    string.AppendFormat("%.1lf\t", fluxResult.m_plumeHeight.m_plumeAltitudeError);

    // 9. Output where we got the plume height from 
    string.AppendFormat("%s\t", (const char*)phSrc);

    // 10. Output the compass direction
    string.AppendFormat("%.2lf\t", fluxResult.m_compass);

    // 11. Output where we got the compass direction from
    //if(fabs(spectrometer.m_scanner.compass) > 360)
    //	string.AppendFormat("compassreading\t");
    //else
    //	string.AppendFormat("user\t");

    // 12. Output the plume centre
    string.AppendFormat("%.1lf\t", result->GetCalculatedPlumeCentre());

    // 13. Output the plume completeness
    string.AppendFormat("%.2lf\t", result->GetCalculatedPlumeCompleteness());

    // 14. Output the cone-angle of the scanner
    string.AppendFormat("%.1lf\t", fluxResult.m_coneAngle);

    // 15. Output the tilt of the scanner
    string.AppendFormat("%.1lf\t", fluxResult.m_tilt);

    // 16. Output whether we think this is a good measurement or not
    if (result->IsFluxOk())
        string.AppendFormat("1\t");
    else
        string.AppendFormat("0\t");

    // 17. Output the instrument temperature
    string.AppendFormat("%.1lf\t", result->GetTemperature());

    // 18. Output the instrument battery voltage
    string.AppendFormat("%.1lf\t", result->GetBatteryVoltage());

    // 19. Output the exposure-time
    string.AppendFormat("%.1ld\t", result->GetSkySpectrumInfo().m_exposureTime);

    // 20. Find the name of the flux-log file to write to

    // 20a. Make the directory
    serialNumber.Format("%s", (const char*)result->GetSerial());
    directory.Format("%s%s%c%s%c", (const char*)g_userSettings.m_outputDirectory, (const char*)dateStr2,
        Poco::Path::separator(), (const char*)serialNumber, Poco::Path::separator());
    if (Filesystem::CreateDirectoryStructure(directory))
    {
        Common common;
        directory.Format("%sOutput%c%s%c%s", (const char*)common.m_exePath, Poco::Path::separator(), (const char*)dateStr2, Poco::Path::separator(), (const char*)serialNumber);
        if (Filesystem::CreateDirectoryStructure(directory))
        {
            errorMessage.Format("Could not create storage directory for flux-data. Please check settings and restart.");
            ShowError(errorMessage);
            return RETURN_CODE::FAIL;
        }
    }

    // 20b. Get the file-name
    fluxLogFile.Format("%sFluxLog_%s_%s.txt", (const char*)directory, (const char*)serialNumber, (const char*)dateStr2);

    // 20c. Check if the file exists
    if (!Filesystem::IsExistingFile(fluxLogFile))
    {
        // write the header
        FILE* f = fopen(fluxLogFile, "w");
        if (f != nullptr)
        {
            fprintf(f, "serial=%s\n", (const char*)serialNumber);
            fprintf(f, "volcano=x\n");//,		m_common.SimplifyString(spectrometer.m_scanner.volcano));
            fprintf(f, "site=x\n");//,				m_common.SimplifyString(spectrometer.m_scanner.site));
            fprintf(f, "#scandate\tscanstarttime\tscanstoptime\t");
            fprintf(f, "flux_[kg/s]\t");
            fprintf(f, "windspeed_[m/s]\twinddirection_[deg]\twindspeedsource\twinddirectionsource\t");
            fprintf(f, "plumeheight_[m]\tplumeheightsource\t");
            fprintf(f, "compassdirection_[deg]\tcompasssource\t");
            fprintf(f, "plumecentre_[deg]\tplumecompleteness_[%%]\t");
            fprintf(f, "coneangle\ttilt\tokflux\ttemperature\tbatteryvoltage\texposuretime\n");
            fclose(f);
        }
    }

    // 20d. Write the flux-result to the file
    FILE* f = fopen(fluxLogFile, "a+");
    if (f != nullptr)
    {
        fprintf(f, "%s\n", (const char*)string);
        fclose(f);
    }

    return RETURN_CODE::SUCCESS;
}
