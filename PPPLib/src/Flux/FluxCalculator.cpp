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
#include <limits>
#include <assert.h>
#include <sstream>

#undef min
#undef max

namespace Flux
{
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
    novac::CString errorMessage, serial;

    context = context.With("file", novac::GetFileName(evalLogFileName));

    if (!Filesystem::IsExistingFile(evalLogFileName))
    {
        m_log.Error(context, "Recieved evaluation-log with illegal path. Could not calculate flux.");
        return false;
    }

    // Get some information about the scan from the file-name
    novac::CDateTime skyStartTime;
    int channel;
    MEASUREMENT_MODE mode;
    const std::string shortFileName = novac::GetFileName(evalLogFileName);
    if (!novac::CFileUtils::GetInfoFromFileName(shortFileName, skyStartTime, serial, channel, mode))
    {
        m_log.Error(context, "Failed to read instrument information from the name of the evaluation log.");
    }

    context = context.With("instrument", serial.std_str());

    // Find the location of this instrument
    Configuration::CInstrumentLocation instrLocation;
    if (GetLocation(context, serial, skyStartTime, instrLocation))
    {
        return false;
    }
    if (instrLocation.m_coneangle < 45.0)
    {
        m_log.Error(context, "Invalid cone angle in setup. Cannot calculate flux");
        return false;
    }

    // Get the wind field at the time of the collection of this scan
    Meteorology::CWindField windField;
    if (!windDataBase.GetWindField(skyStartTime, novac::CGPSData(instrLocation.m_latitude, instrLocation.m_longitude, plumeAltitude.m_plumeAltitude), Meteorology::INTERP_NEAREST_NEIGHBOUR, windField))
    {
        m_log.Error(context, "Could not retrieve wind field at time of measurement. Could not calculate flux");
        return false;
    }

    // 4b. Adjust the altitude of the plume so that it is relative to the altitude of the instrument...
    Geometry::CPlumeHeight relativePlumeHeight = plumeAltitude;
    relativePlumeHeight.m_plumeAltitude -= instrLocation.m_altitude;
    if (relativePlumeHeight.m_plumeAltitude <= 0)
    {
        m_log.Error(context, "Negative plume height obtained when calculating flux. No flux can be calculated.");
        return false;
    }

    // Read in the evaluation log file 
    FileHandler::CEvaluationLogFileHandler reader(m_log, evalLogFileName, m_userSettings.m_molecule);
    if (RETURN_CODE::SUCCESS != reader.ReadEvaluationLog())
    {
        m_log.Error("Failed to read evaluation log");
        return false;
    }
    if (reader.m_scan.size() == 0)
    {
        m_log.Error("Recieved evaluation log file with no scans inside. Cannot calculate flux");
        return false;
    }
    else if (reader.m_scan.size() > 1)
    {
        m_log.Error("Recieved evaluation log file with more than one scans inside. Can only calculate flux for the first scan.");
    }

    // Extract the scan
    Evaluation::CScanResult& result = reader.m_scan[0];

    // Improve on the start-time of the scan...
    result.GetSkyStartTime(skyStartTime);

    // Calculate the offset of the scan
    if (result.CalculateOffset(CMolecule(m_userSettings.m_molecule)))
    {
        m_log.Information(context, "Could not calculate offset for scan. No flux can be calculated.");
        return false;
    }

    // Check that the completeness is higher than our limit...
    if (!result.CalculatePlumeCentre(CMolecule(m_userSettings.m_molecule)))
    {
        m_log.Information(context, "Scan does not see the plume, no flux can be calculated");
        return false;
    }

    const double completeness = result.GetCalculatedPlumeCompleteness();
    if (completeness < (m_userSettings.m_completenessLimitFlux + 0.01))
    {
        errorMessage.Format("Scan has completeness = %.2lf which is less than limit of %.2lf. Rejected!", completeness, m_userSettings.m_completenessLimitFlux);
        m_log.Information(context, errorMessage.std_str());
        return false;
    }

    // Calculate the flux
    if (!CalculateFlux(context, result, CMolecule(m_userSettings.m_molecule), windField, relativePlumeHeight, instrLocation.m_compass, instrLocation.m_coneangle, instrLocation.m_tilt))
    {
        m_log.Information(context, "Flux calculation failed. No flux generated");
        return false;
    }
    fluxResult = result.GetFluxResult();

    // Make a simple estimate of the quality of the flux measurement
    FluxQuality windFieldQuality = FluxQuality::Green;
    FluxQuality plumeHeightQuality = FluxQuality::Green;
    FluxQuality completenessQuality = FluxQuality::Green;

    switch (windField.GetWindSpeedSource())
    {
    case Meteorology::MET_DEFAULT:				windFieldQuality = FluxQuality::Red; break;
    case Meteorology::MET_USER:					windFieldQuality = FluxQuality::Red; break;
    case Meteorology::MET_ECMWF_FORECAST:		windFieldQuality = FluxQuality::Green; break;
    case Meteorology::MET_ECMWF_ANALYSIS:		windFieldQuality = FluxQuality::Green; break;
    case Meteorology::MET_DUAL_BEAM_MEASUREMENT:windFieldQuality = FluxQuality::Green; break;
    case Meteorology::MET_MODEL_WRF:			windFieldQuality = FluxQuality::Green; break;
    case Meteorology::MET_NOAA_GDAS:			windFieldQuality = FluxQuality::Green; break;
    case Meteorology::MET_NOAA_FNL:				windFieldQuality = FluxQuality::Green; break;
    default:									windFieldQuality = FluxQuality::Yellow; break;
    }
    switch (relativePlumeHeight.m_plumeAltitudeSource)
    {
    case Meteorology::MET_DEFAULT:				plumeHeightQuality = FluxQuality::Red; break;
    case Meteorology::MET_USER:					plumeHeightQuality = FluxQuality::Red; break;
    case Meteorology::MET_GEOMETRY_CALCULATION:	plumeHeightQuality = FluxQuality::Green; break;
    default:									plumeHeightQuality = FluxQuality::Yellow; break;
    }
    if (fluxResult.m_completeness < 0.7)
    {
        completenessQuality = FluxQuality::Red;
    }
    else if (fluxResult.m_completeness < 0.9)
    {
        completenessQuality = FluxQuality::Yellow;
    }
    else
    {
        completenessQuality = FluxQuality::Green;
    }

    // If any of the parameters has a low quality, then the result is judged to have a low quality...
    if (windFieldQuality == FluxQuality::Red || plumeHeightQuality == FluxQuality::Red || completenessQuality == FluxQuality::Red)
    {
        fluxResult.m_fluxQualityFlag = FluxQuality::Red;
    }
    else if (windFieldQuality == FluxQuality::Yellow || plumeHeightQuality == FluxQuality::Yellow || completenessQuality == FluxQuality::Yellow)
    {
        fluxResult.m_fluxQualityFlag = FluxQuality::Yellow;
    }
    else
    {
        fluxResult.m_fluxQualityFlag = FluxQuality::Green;
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

// Moved from CScanResult::CalculateFlux
bool CFluxCalculator::CalculateFlux(
    novac::LogContext context,
    Evaluation::CScanResult& result,
    const CMolecule& specie,
    const Meteorology::CWindField& wind,
    const Geometry::CPlumeHeight& relativePlumeHeight,
    double compass,
    double coneAngle,
    double tilt)
{
    Meteorology::CWindField modifiedWind;

    // If this is a not a flux measurement, then don't calculate any flux
    if (!result.IsFluxMeasurement())
    {
        m_log.Information(context, "Measurement is not a flux measurement, cannot calculate flux.");
        return false;
    }

    // get the specie index
    int specieIndex = result.GetSpecieIndex(specie.name);
    if (specieIndex == -1)
    {
        novac::CString message;
        message.Format("Failed to retrieve specie with name '%s' from measurement, cannot calculate flux.", specie.name.c_str());
        m_log.Information(context, message.std_str());
        return false;
    }

    // pull out the good data points out of the measurement and ignore the bad points
    // at the same time convert to mg/m2
    std::vector<double> scanAngle;
    std::vector<double> scanAngle2;
    std::vector<double> column;
    scanAngle.reserve(result.GetEvaluatedNum());
    scanAngle2.reserve(result.GetEvaluatedNum());
    column.reserve(result.GetEvaluatedNum());
    unsigned int numberOfGoodDataPoints = 0;
    for (unsigned long i = 0; i < result.GetEvaluatedNum(); ++i)
    {
        if (result.IsBad(i) || result.IsDeleted(i))
        {
            continue; // this is a bad measurement
        }
        if (result.m_specInfo[i].m_flag >= 64)
        {
            continue; // this is a direct-sun measurement, don't use it to calculate the flux...
        }

        scanAngle.push_back(result.m_specInfo[i].m_scanAngle);
        scanAngle2.push_back(result.m_specInfo[i].m_scanAngle2);
        column.push_back(specie.Convert_MolecCm2_to_kgM2(result.m_spec[i].m_referenceResult[specieIndex].m_column));

        ++numberOfGoodDataPoints;
    }

    Flux::CFluxResult fluxResult;

    // if there are no good datapoints in the measurement, the flux is assumed to be zero
    if (numberOfGoodDataPoints < 10)
    {
        if (numberOfGoodDataPoints == 0)
        {
            m_log.Information(context, "Could not calculate flux, no good datapoints in measurement");
        }
        else
        {
            m_log.Information(context, "Could not calculate flux, too few good datapoints in measurement");
        }

        result.m_flux = fluxResult;
        return false;
    }

    // Get the times of the scan
    novac::CDateTime startTime1 = result.GetSkyStartTime();
    novac::CDateTime startTime2;
    result.GetStartTime(0, startTime2);
    if (startTime1.year != 0 && startTime1 < startTime2)
    {
        fluxResult.m_startTime = startTime1;
    }
    else
    {
        fluxResult.m_startTime = startTime2;
    }

    result.GetStopTime(result.GetEvaluatedNum() - 1, fluxResult.m_stopTime);

    // and the serial number of the instrument
    fluxResult.m_instrument = result.GetSerial().std_str();

    // Calculate the flux
    fluxResult.m_flux = Flux::CFluxCalculator::CalculateFlux(context, scanAngle.data(), scanAngle2.data(), column.data(), specie.Convert_MolecCm2_to_kgM2(result.m_plumeProperties.offset), numberOfGoodDataPoints, wind, relativePlumeHeight, compass, result.m_instrumentType, coneAngle, tilt);
    fluxResult.m_windField = wind;
    fluxResult.m_plumeHeight = relativePlumeHeight;
    fluxResult.m_compass = compass;
    fluxResult.m_coneAngle = coneAngle;
    fluxResult.m_tilt = tilt;
    fluxResult.m_numGoodSpectra = numberOfGoodDataPoints;

    if (std::isnan(fluxResult.m_flux))
    {
        return false;
    }

    // calculate the completeness and centre of the plume
    result.CalculatePlumeCentre(specie);

    fluxResult.m_scanOffset = result.m_plumeProperties.offset;
    fluxResult.m_completeness = result.m_plumeProperties.completeness;
    fluxResult.m_plumeCentre[0] = result.m_plumeProperties.plumeCenter;
    fluxResult.m_plumeCentre[1] = result.m_plumeProperties.plumeCenter2;
    fluxResult.m_instrumentType = result.m_instrumentType;

    // Try to make an estimation of the error in flux from the
    //  wind field used and from the plume height used

    // 1. the wind field
    modifiedWind = wind;
    modifiedWind.SetWindDirection(wind.GetWindDirection() - wind.GetWindDirectionError(), wind.GetWindDirectionSource());
    double flux1 = Flux::CFluxCalculator::CalculateFlux(context, scanAngle.data(), scanAngle2.data(), column.data(), specie.Convert_MolecCm2_to_kgM2(result.m_plumeProperties.offset), numberOfGoodDataPoints, modifiedWind, relativePlumeHeight, compass, result.m_instrumentType, coneAngle, tilt);

    modifiedWind.SetWindDirection(wind.GetWindDirection() + wind.GetWindDirectionError(), wind.GetWindDirectionSource());
    double flux2 = Flux::CFluxCalculator::CalculateFlux(context, scanAngle.data(), scanAngle2.data(), column.data(), specie.Convert_MolecCm2_to_kgM2(result.m_plumeProperties.offset), numberOfGoodDataPoints, modifiedWind, relativePlumeHeight, compass, result.m_instrumentType, coneAngle, tilt);

    double fluxErrorDueToWindDirection = std::max(std::abs(flux2 - fluxResult.m_flux), std::abs(flux1 - fluxResult.m_flux));

    double fluxErrorDueToWindSpeed = fluxResult.m_flux * wind.GetWindSpeedError() / wind.GetWindSpeed();

    fluxResult.m_fluxError_Wind = sqrt(fluxErrorDueToWindDirection * fluxErrorDueToWindDirection + fluxErrorDueToWindSpeed * fluxErrorDueToWindSpeed);

    // 2. the plume height
    fluxResult.m_fluxError_PlumeHeight = fluxResult.m_flux * relativePlumeHeight.m_plumeAltitudeError / relativePlumeHeight.m_plumeAltitude;

    result.m_flux = fluxResult;

    return true;
}

double CFluxCalculator::CalculateFlux(
    novac::LogContext context,
    const double* scanAngle,
    const double* scanAngle2,
    const double* column,
    double offset,
    int nDataPoints,
    const Meteorology::CWindField& wind,
    const Geometry::CPlumeHeight& relativePlumeHeight,
    double compass,
    INSTRUMENT_TYPE type,
    double coneAngle,
    double tilt)
{
    const double windSpeed = wind.GetWindSpeed();
    const double windDirection = wind.GetWindDirection();
    const double plumeHeight = relativePlumeHeight.m_plumeAltitude;

    assert(!std::isnan(windSpeed));
    assert(!std::isnan(windDirection));
    assert(!std::isnan(plumeHeight));

    if (type == INSTRUMENT_TYPE::INSTR_HEIDELBERG)
    {
        m_log.Debug(context, "Calculating flux for Mark-II type instrument");
        return CalculateFluxHeidelbergScanner(scanAngle, scanAngle2, column, offset, nDataPoints, windSpeed, windDirection, plumeHeight, compass);
    }
    else if (type == INSTRUMENT_TYPE::INSTR_GOTHENBURG)
    {
        if (coneAngle < 10.0 || coneAngle > 91.0)
        {
            m_log.Error(context, "Invalid cone angle in setup");
            return std::numeric_limits<double>::quiet_NaN();
        }

        // In the NovacPPP, the gas factor isn't used, since the colums should have been converted from molec/cm2
        // into kg/m2 already before calling this function. Hence the conversion to mass has already been done.
        // However the flux-calculation formula, shared with the NovacProgram, requires the gas factor.
        // This compensation factor is used to compensate for how the gas factor is weighted into the calculation...
        const double gasFactorCompensation = 1e6;
        if (std::abs(coneAngle - 90.0) < 1.0)
        {
            return CalculateFluxFlatScanner(
                scanAngle,
                column,
                offset,
                nDataPoints,
                windSpeed,
                windDirection,
                plumeHeight,
                compass,
                gasFactorCompensation);
        }
        else
        {
            return CalculateFluxConicalScanner(
                scanAngle,
                column,
                offset,
                nDataPoints,
                windSpeed,
                windDirection,
                plumeHeight,
                compass,
                gasFactorCompensation,
                coneAngle,
                tilt);
        }
    }

    std::stringstream msg;
    msg << "Not supported instrument type: " << (int)type;
    throw std::invalid_argument(msg.str());
}

// endregion The actual flux calculations
}
