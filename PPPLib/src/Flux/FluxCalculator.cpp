#include <PPPLib/Flux/FluxCalculator.h>
#include <PPPLib/Logging.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Flux/Flux.h>

// ... support for handling the evaluation-log files...
#include <PPPLib/File/EvaluationLogFileHandler.h>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

#include <PPPLib/File/Filesystem.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <Poco/Path.h>
#include <cmath>

using namespace Flux;

CFluxCalculator::CFluxCalculator(
    novac::ILogger& log,
    const Configuration::CNovacPPPConfiguration& setup,
    const Configuration::CUserConfiguration& userSettings)
    : m_log(log), m_setup(setup), m_userSettings(userSettings)
{
}

bool CFluxCalculator::CalculateFlux(
    novac::LogContext context,
    const std::string& evalLogFileName,
    const Meteorology::CWindDataBase& windDataBase,
    const Geometry::CPlumeHeight& plumeAltitude,
    CFluxResult& fluxResult)
{
    novac::CDateTime skyStartTime;
    novac::CString errorMessage, serial;
    Geometry::CPlumeHeight relativePlumeHeight;
    Meteorology::CWindField windField;
    int channel;
    MEASUREMENT_MODE mode;

    context = context.With("file", novac::GetFileName(evalLogFileName));

    if (!Filesystem::IsExistingFile(evalLogFileName))
    {
        m_log.Error(context, "Recieved evaluation-log with illegal path. Could not calculate flux.");
        return false;
    }

    // 2. Get some information about the scan from the file-name
    const std::string shortFileName = novac::GetFileName(evalLogFileName);
    novac::CFileUtils::GetInfoFromFileName(shortFileName, skyStartTime, serial, channel, mode);

    context = context.With("instrument", serial.std_str());

    // 3. Find the location of this instrument
    Configuration::CInstrumentLocation instrLocation;
    if (GetLocation(context, serial, skyStartTime, instrLocation))
    {
        return false;
    }

    // 4. Get the wind field at the time of the collection of this scan
    if (!windDataBase.GetWindField(skyStartTime, novac::CGPSData(instrLocation.m_latitude, instrLocation.m_longitude, plumeAltitude.m_plumeAltitude), Meteorology::INTERP_NEAREST_NEIGHBOUR, windField))
    {
        m_log.Error(context, "Could not retrieve wind field at time of measurement. Could not calculate flux");
        return false;
    }

    // 4b. Adjust the altitude of the plume so that it is relative to the altitude of the instrument...
    relativePlumeHeight = plumeAltitude;
    relativePlumeHeight.m_plumeAltitude -= instrLocation.m_altitude;
    if (relativePlumeHeight.m_plumeAltitude <= 0)
    {
        m_log.Error(context, "Negative plume height obtained when calculating flux. No flux can be calculated.");
        return false;
    }

    // 5. Read in the evaluation log file 
    FileHandler::CEvaluationLogFileHandler reader(m_log, evalLogFileName, m_userSettings.m_molecule);
    reader.ReadEvaluationLog();
    if (reader.m_scan.size() == 0)
    {
        m_log.Error("Recieved evaluation log file with no scans inside. Cannot calculate flux");
        return false;
    }
    else if (reader.m_scan.size() > 1)
    {
        m_log.Error("Recieved evaluation log file with more than one scans inside. Can only calculate flux for the first scan.");
    }

    // 6. extract the scan
    Evaluation::CScanResult& result = reader.m_scan[0];

    // 6b. Improve on the start-time of the scan...
    result.GetSkyStartTime(skyStartTime);

    // 7. Calculate the offset of the scan
    if (result.CalculateOffset(CMolecule(m_userSettings.m_molecule)))
    {
        m_log.Information(context, "Could not calculate offset for scan. No flux can be calculated.");
        return false;
    }

    // 8. Check that the completeness is higher than our limit...
    if (!result.CalculatePlumeCentre(CMolecule(m_userSettings.m_molecule)))
    {
        m_log.Information(context, "Scan does not see the plume, no flux can be calculated");
        return false;
    }
    double completeness = result.GetCalculatedPlumeCompleteness();
    if (completeness < (m_userSettings.m_completenessLimitFlux + 0.01))
    {
        errorMessage.Format("Scan has completeness = %.2lf which is less than limit of %.2lf. Rejected!", completeness, m_userSettings.m_completenessLimitFlux);
        m_log.Information(context, errorMessage.std_str());
        return false;
    }

    // 9. Calculate the flux
    if (result.CalculateFlux(CMolecule(m_userSettings.m_molecule), windField, relativePlumeHeight, instrLocation.m_compass, instrLocation.m_coneangle, instrLocation.m_tilt))
    {
        m_log.Information(context, "Flux calculation failed. No flux generated");
        return false;
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
    return true;
}

int CFluxCalculator::GetLocation(
    novac::LogContext context,
    const novac::CString& serial,
    const novac::CDateTime& startTime,
    Configuration::CInstrumentLocation& instrLocation)
{
    novac::CDateTime day, evalValidFrom, evalValidTo;
    const Configuration::CInstrumentConfiguration* instrumentConf = nullptr;
    Configuration::CInstrumentLocation singleLocation;
    novac::CString errorMessage;

    // First of all find the instrument 
    for (int k = 0; k < m_setup.NumberOfInstruments(); ++k)
    {
        if (Equals(m_setup.m_instrument[k].m_serial, serial))
        {
            instrumentConf = &m_setup.m_instrument[k];
            break;
        }
    }
    if (instrumentConf == nullptr)
    {
        m_log.Error(context, "Recieved spectrum from not-configured instrument. Cannot calculate flux!");
        return 1;
    }

    // Next find the instrument location that is valid for this date
    const Configuration::CLocationConfiguration& locationconf = instrumentConf->m_location;
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
        m_log.Error(context, "Recieved spectrum from not-configured instrument. Cannot calculate flux!");
        return 1;
    }

    // we're done! Return!
    return 0;
}

RETURN_CODE CFluxCalculator::WriteFluxResult(
    novac::LogContext context,
    const Flux::CFluxResult& fluxResult,
    const Evaluation::CScanResult* result)
{
    novac::CString string, dateStr, dateStr2, serialNumber;
    novac::CString fluxLogFile, directory;
    novac::CString wdSrc, wsSrc, phSrc;
    novac::CString errorMessage;
    novac::CDateTime dateTime;

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
    directory.Format("%s%s%c%s%c", (const char*)m_userSettings.m_outputDirectory, (const char*)dateStr2,
        Poco::Path::separator(), (const char*)serialNumber, Poco::Path::separator());
    if (Filesystem::CreateDirectoryStructure(directory))
    {
        m_log.Error(context, "Could not create storage directory for flux-data. Please check settings and restart.");
        return RETURN_CODE::FAIL;
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

// region The actual flux calculations


double CFluxCalculator::CalculateFlux(const double* scanAngle, const double* scanAngle2, const double* column, double offset, int nDataPoints, const Meteorology::CWindField& wind, const Geometry::CPlumeHeight& relativePlumeHeight, double compass, INSTRUMENT_TYPE type, double coneAngle, double tilt)
{
    double windSpeed = wind.GetWindSpeed();
    double windDirection = wind.GetWindDirection();
    double plumeHeight = relativePlumeHeight.m_plumeAltitude;

    if (type == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
    {
        return CalculateFluxHeidelbergScanner(scanAngle, scanAngle2, column, offset, nDataPoints, windSpeed, windDirection, plumeHeight, compass);
    }
    else if (type == INSTRUMENT_TYPE::INSTR_GOTHENBURG)
    {
        // In the NovacPPP, the gas factor isn't used. However the flux-calculation formula, shared with the NovacProgram, requires the gas factor.
        //  This compensation factor is used to compensate for how the gas factor is weighted into the calculation...
        const double gasFactorCompensation = 1e6;
        if (fabs(coneAngle - 90.0) < 1.0)
        {
            return CalculateFluxFlatScanner(scanAngle, column, offset, nDataPoints, windSpeed, windDirection, plumeHeight, compass, gasFactorCompensation);
        }
        else
        {
            return CalculateFluxConicalScanner(scanAngle, column, offset, nDataPoints, windSpeed, windDirection, plumeHeight, compass, coneAngle, tilt, gasFactorCompensation);
        }
    }
    else
    {
        return 0.0; // unsupported instrument-type
    }
}


// endregion The actual flux calculations