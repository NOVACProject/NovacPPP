#include <PPPLib/File/EvaluationLogFileHandler.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/MFC/CSingleLock.h>

#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/StringUtils.h>

#include <cstring>
#include <algorithm>
#include <cmath>
#include <iostream>

// Global variables;
novac::CCriticalSection g_evalLogCritSect; // synchronization access to evaluation-log files

using namespace FileHandler;

CEvaluationLogFileHandler::CEvaluationLogFileHandler(
    novac::ILogger& log,
    std::string evaluationLog,
    novac::Molecule molecule,
    novac::SpectrometerModel* spectrometerModel)
    : m_evaluationLog(evaluationLog), m_molecule(molecule), m_log(log)
{
    // Defining which column contains which information
    m_col.position = 1;
    m_col.position2 = -1; // azimuth does not exist in the typical eval-log files
    m_col.starttime = 2;
    m_col.stoptime = 3;
    m_col.name = -1; // the name does not exist in the original version
    m_col.delta = 4;
    m_col.chiSquare = -1; // Chi-square does not exist in the original version
    m_col.expTime = 5;
    m_col.nSpec = 6;
    m_col.intensity = -1;      // intensity does not exist in the original version
    m_col.fitIntensity = -1;   // fit intensity only exists in the reevaluation logs
    m_col.peakSaturation = -1; // peak and fit saturation does not exist in the original version
    m_col.fitSaturation = -1;  // peak and fit saturation does not exist in the original version
    m_col.offset = -1;         // offset does not exist in the original version
    m_col.nSpecies = 0;
    for (int i = 0; i < MAX_N_REFERENCES; ++i)
    {
        m_col.column[i] = 7;
        m_col.columnError[i] = 8;
        m_col.shift[i] = 9;
        m_col.shiftError[i] = 10;
        m_col.squeeze[i] = 11;
        m_col.squeezeError[i] = 12;
    }

    m_instrumentType = novac::NovacInstrumentType::Gothenburg;

    if (spectrometerModel != nullptr)
    {
        m_spectrometerModel = *spectrometerModel;
    }
    else
    {
        m_spectrometerModel = novac::CSpectrometerDatabase::SpectrometerModel_Unknown();
    }
}

void CEvaluationLogFileHandler::ParseScanHeader(const char szLine[8192])
{
    // reset some old information
    ResetColumns();

    char str[8192];
    if (szLine[0] == '#')
        strncpy(str, szLine + 1, 8191 * sizeof(char));
    else
        strncpy(str, szLine, 8192 * sizeof(char));

    char* szToken = str;
    int curCol = -1;
    char elevation[] = "elevation";
    char scanAngle[] = "scanangle";
    char obsAngle[] = "observationangle";
    char azimuth[] = "azimuth";
    char column[] = "column";
    char columnError[] = "columnerror";
    char intensity[] = "intensity";            // peak intensity
    char fitIntensity[] = "intens(fitregion)"; // fit-region intensity
    char fitIntensity2[] = "fitintensity";     // fit-region intensity
    char peakSat[] = "specsaturation";         // maximum saturation ratio of the whole spectrum
    char fitSat[] = "fitsaturation";           // maximum saturation ratio in the fit region
    char delta[] = "delta";
    char chiSquare[] = "chisquare";
    char shift[] = "shift";
    char shiftError[] = "shifterror";
    char squeeze[] = "squeeze";
    char squeezeError[] = "squeezeerror";
    char exposureTime[] = "exposuretime";
    char numSpec[] = "numSpec";
    char offset[] = "offset";
    char starttime[] = "starttime";
    char stoptime[] = "stoptime";
    char nameStr[] = "name";

    while (nullptr != (szToken = strtok(szToken, "\t")))
    {
        ++curCol;

        // The scan-angle (previously known as elevation)
        if (novac::Equals(szToken, elevation, strlen(elevation)))
        {
            m_col.position = curCol;
            szToken = nullptr;
            continue;
        }

        // The scan-angle (previously known as elevation)
        if (novac::Equals(szToken, scanAngle, strlen(scanAngle)))
        {
            m_col.position = curCol;
            szToken = nullptr;
            continue;
        }

        // The observation-angle (the scan-angle for the heidelberg instrument)
        if (novac::Equals(szToken, obsAngle, strlen(obsAngle)))
        {
            m_col.position = curCol;
            szToken = nullptr;
            continue;
        }

        // The azimuth-angle (defined for the heidelberg instrument)
        if (novac::Equals(szToken, azimuth, strlen(azimuth)))
        {
            m_col.position2 = curCol;
            szToken = nullptr;
            continue;
        }

        // The exposure time
        if (novac::Equals(szToken, exposureTime, strlen(exposureTime)))
        {
            m_col.expTime = curCol;
            szToken = nullptr;
            continue;
        }

        // The start time
        if (novac::Equals(szToken, starttime, strlen(starttime)))
        {
            m_col.starttime = curCol;
            szToken = nullptr;
            continue;
        }

        // The stop time
        if (novac::Equals(szToken, stoptime, strlen(stoptime)))
        {
            m_col.stoptime = curCol;
            szToken = nullptr;
            continue;
        }

        // The name of the spectrum
        if (novac::Equals(szToken, nameStr, strlen(nameStr)))
        {
            m_col.name = curCol;
            szToken = nullptr;
            continue;
        }

        // The number of co-added spectra
        if (novac::Equals(szToken, numSpec, strlen(numSpec)))
        {
            m_col.nSpec = curCol;
            szToken = nullptr;
            continue;
        }

        // The offset
        if (novac::Equals(szToken, offset, strlen(offset)))
        {
            m_col.offset = curCol;
            szToken = nullptr;
            continue;
        }

        // The column error (must be looked for before 'column')
        if (novac::Equals(szToken, columnError, strlen(columnError)))
        {
            m_col.columnError[m_evResult.m_referenceResult.size() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The column
        if (novac::Equals(szToken, column, strlen(column)))
        {
            m_col.column[m_evResult.m_referenceResult.size()] = curCol;
            char* pt = szToken + strlen(column) + 1;
            szToken[strlen(szToken) - 1] = 0;
            std::string specieStr(pt);
            std::transform(specieStr.begin(), specieStr.end(), specieStr.begin(), ::toupper);
            m_evResult.InsertSpecie(specieStr);
            ++m_col.nSpecies;
            szToken = nullptr;
            continue;
        }

        // The shift error (must be checked before 'shift')
        if (novac::Equals(szToken, shiftError, strlen(shiftError)))
        {
            m_col.shiftError[m_evResult.m_referenceResult.size() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The shift
        if (novac::Equals(szToken, shift, strlen(shift)))
        {
            m_col.shift[m_evResult.m_referenceResult.size() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The squeeze error (must be checked before 'squeeze')
        if (novac::Equals(szToken, squeezeError, strlen(squeezeError)))
        {
            m_col.squeezeError[m_evResult.m_referenceResult.size() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The squeeze
        if (novac::Equals(szToken, squeeze, strlen(squeeze)))
        {
            m_col.squeeze[m_evResult.m_referenceResult.size() - 1] = curCol;
            szToken = nullptr;
            continue;
        }

        // The spectrum peak-intensity
        if (novac::Equals(szToken, intensity, strlen(intensity)))
        {
            m_col.intensity = curCol;
            szToken = nullptr;
            continue;
        }

        // The spectrum fit-intensity
        if (novac::Equals(szToken, fitIntensity, strlen(fitIntensity)) ||
            novac::Equals(szToken, fitIntensity2, strlen(fitIntensity2)))
        {

            m_col.fitIntensity = curCol;
            szToken = nullptr;
            continue;
        }

        // The spectrum maximum saturation ratio of the whole spectrum
        if (novac::Equals(szToken, peakSat, strlen(peakSat)))
        {
            m_col.peakSaturation = curCol;
            szToken = nullptr;
            continue;
        }

        // The spectrum maximum saturation ratio in the fit region
        if (novac::Equals(szToken, fitSat, strlen(fitSat)))
        {
            m_col.fitSaturation = curCol;
            szToken = nullptr;
            continue;
        }

        // The delta of the fit
        if (novac::Equals(szToken, delta, strlen(delta)))
        {
            m_col.delta = curCol;
            szToken = nullptr;
            continue;
        }

        // The chi-square of the fit
        if (novac::Equals(szToken, chiSquare, strlen(chiSquare)))
        {
            m_col.chiSquare = curCol;
            szToken = nullptr;
            continue;
        }

        szToken = nullptr;
    }

    for (size_t k = 0; k < m_evResult.m_referenceResult.size(); ++k)
    {
        m_specieName.push_back(m_evResult.m_referenceResult[k].m_specieName);
    }

    return;
}

RETURN_CODE CEvaluationLogFileHandler::ReadEvaluationLog()
{
    char expTimeStr[] = "exposuretime";           // this string only exists in the header line.
    char scanInformation[] = "<scaninformation>"; // this string only exists in the scan-information section before the scan-data
    char fluxInformation[] = "<fluxinfo>";        // this string only exists in the flux-information section before the scan-data
    char spectralData[] = "<spectraldata>";
    char endofSpectralData[] = "</spectraldata>";
    novac::CString str;
    char szLine[8192];
    int measNr = 0;
    double fValue;
    bool fReadingScan = false;
    double flux = 0.0;

    // If no evaluation log selected, quit
    if (m_evaluationLog.size() <= 1)
        return RETURN_CODE::FAIL;

    // First count the number of scans in the file.
    //	This to speed up the initialization of the arrays
    const size_t nScans = CountScansInFile();
    if (nScans > 0)
    {
        m_scan.reserve(nScans);
        m_windField.reserve(nScans + 1);
    }

    Evaluation::CScanResult newResult; // this is the scan we're reading in right now

    // Open the evaluation log
    novac::CSingleLock singleLock(&g_evalLogCritSect);
    singleLock.Lock();
    if (singleLock.IsLocked())
    {
        FILE* f = fopen(m_evaluationLog.c_str(), "r");
        if (nullptr == f)
        {
            singleLock.Unlock();
            return RETURN_CODE::FAIL;
        }

        // Reset the column- and spectrum info
        ResetColumns();
        ResetScanInformation();

        // Read the file, one line at a time
        while (fgets(szLine, 8192, f))
        {

            // ignore empty lines
            if (strlen(szLine) < 2)
            {
                if (fReadingScan)
                {
                    fReadingScan = false;
                    // Reset the column- and spectrum-information
                    ResetColumns();
                    ResetScanInformation();
                }
                continue;
            }

            // convert the string to all lower-case letters
            for (unsigned int it = 0; it < strlen(szLine); ++it)
            {
                szLine[it] = (char)tolower(szLine[it]);
            }

            // find the next scan-information section
            if (nullptr != strstr(szLine, scanInformation))
            {
                ResetScanInformation();
                ParseScanInformation(m_specInfo, flux, f);
                newResult.m_skySpecInfo = m_specInfo;
                continue;
            }

            // find the next flux-information section
            if (nullptr != strstr(szLine, fluxInformation))
            {
                Meteorology::WindField windField;
                ParseFluxInformation(windField, flux, f);
                m_windField.push_back(windField);
                continue;
            }

            if (nullptr != strstr(szLine, spectralData))
            {
                fReadingScan = true;
                continue;
            }
            else if (nullptr != strstr(szLine, endofSpectralData))
            {
                fReadingScan = false;
                continue;
            }

            // find the next start of a scan
            if (nullptr != strstr(szLine, expTimeStr))
            {

                // check so that there was some information in the last scan read
                //	if not the re-use the memory space
                if (measNr > 0)
                {
                    // The current measurement position inside the scan
                    measNr = 0;

                    // before we start the next scan, calculate some information about
                    // the old one

                    // 1. If the sky and dark were specified, remove them from the measurement
                    if (m_scan.size() >= 0 && fabs(m_scan.back().GetScanAngle(1) - 180.0) < 1)
                    {
                        m_scan.back().RemoveResult(0); // remove sky
                        m_scan.back().RemoveResult(0); // remove dark
                    }

                    // start the next scan.
                }

                // This line is the header line which says what each column represents.
                //  Read it and parse it to find out how to interpret the rest of the
                //  file.
                ParseScanHeader(szLine);

                // start parsing the lines
                fReadingScan = true;

                // read the next line, which is the first line in the scan
                continue;
            }

            // ignore comment lines
            if (szLine[0] == '#')
                continue;

            // if we're not reading a scan, let's read the next line
            if (!fReadingScan)
                continue;

            // Split the scan information up into tokens and parse them.
            char* szToken = (char*)szLine;
            int curCol = -1;
            while (nullptr != (szToken = strtok(szToken, " \t")))
            {
                ++curCol;

                // First check the starttime
                if (curCol == m_col.starttime)
                {
                    int fValue1, fValue2, fValue3;
                    if (strstr(szToken, ":"))
                    {
                        sscanf(szToken, "%d:%d:%d", &fValue1, &fValue2, &fValue3);
                    }
                    else
                    {
                        sscanf(szToken, "%d.%d.%d", &fValue1, &fValue2, &fValue3);
                    }
                    m_specInfo.m_startTime.hour = (unsigned char)fValue1;
                    m_specInfo.m_startTime.minute = (unsigned char)fValue2;
                    m_specInfo.m_startTime.second = (unsigned char)fValue3;
                    szToken = nullptr;
                    continue;
                }

                // Then check the stoptime
                if (curCol == m_col.stoptime)
                {
                    int fValue1, fValue2, fValue3;
                    if (strstr(szToken, ":"))
                    {
                        sscanf(szToken, "%d:%d:%d", &fValue1, &fValue2, &fValue3);
                    }
                    else
                    {
                        sscanf(szToken, "%d.%d.%d", &fValue1, &fValue2, &fValue3);
                    }
                    m_specInfo.m_stopTime.hour = (unsigned char)fValue1;
                    m_specInfo.m_stopTime.minute = (unsigned char)fValue2;
                    m_specInfo.m_stopTime.second = (unsigned char)fValue3;
                    szToken = nullptr;
                    continue;
                }

                // Also check the name...
                if (curCol == m_col.name)
                {
                    m_specInfo.m_name = std::string(szToken);
                    szToken = nullptr;
                    continue;
                }

                // ignore columns whose value cannot be parsed into a float
                if (1 != sscanf(szToken, "%lf", &fValue))
                {
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.position)
                {
                    m_specInfo.m_scanAngle = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.position2)
                {
                    m_specInfo.m_scanAngle2 = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.intensity)
                {
                    m_specInfo.m_peakIntensity = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.fitIntensity)
                {
                    m_specInfo.m_fitIntensity = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.fitSaturation)
                {
                    m_specInfo.m_fitIntensity = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.peakSaturation)
                {
                    m_specInfo.m_peakIntensity = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.offset)
                {
                    m_specInfo.m_offset = (float)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.delta)
                {
                    m_evResult.m_delta = fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.chiSquare)
                {
                    m_evResult.m_chiSquare = fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.nSpec)
                {
                    m_specInfo.m_numSpec = (long)fValue;
                    szToken = nullptr;
                    continue;
                }

                if (curCol == m_col.expTime)
                {
                    m_specInfo.m_exposureTime = (long)fValue;
                    szToken = nullptr;
                    continue;
                }

                for (size_t k = 0; k < static_cast<size_t>(m_col.nSpecies); ++k)
                {
                    if (curCol == m_col.column[k])
                    {
                        m_evResult.m_referenceResult[k].m_column = fValue;
                        break;
                    }
                    if (curCol == m_col.columnError[k])
                    {
                        m_evResult.m_referenceResult[k].m_columnError = fValue;
                        break;
                    }
                    if (curCol == m_col.shift[k])
                    {
                        m_evResult.m_referenceResult[k].m_shift = fValue;
                        break;
                    }
                    if (curCol == m_col.shiftError[k])
                    {
                        m_evResult.m_referenceResult[k].m_shiftError = fValue;
                        break;
                    }
                    if (curCol == m_col.squeeze[k])
                    {
                        m_evResult.m_referenceResult[k].m_squeeze = fValue;
                        break;
                    }
                    if (curCol == m_col.squeezeError[k])
                    {
                        m_evResult.m_referenceResult[k].m_squeezeError = fValue;
                        break;
                    }
                }
                szToken = nullptr;
            }

            // start reading the next line in the evaluation log (i.e. the next
            //  spectrum in the scan). Insert the data from this spectrum into the
            //  CScanResult structure

            m_specInfo.m_scanIndex = (short)measNr;
            if (novac::Equals(m_specInfo.m_name, "sky"))
            {
                newResult.SetSkySpecInfo(m_specInfo);
            }
            else if (novac::Equals(m_specInfo.m_name, "dark"))
            {
                newResult.SetDarkSpecInfo(m_specInfo);
            }
            else if (novac::Equals(m_specInfo.m_name, "offset"))
            {
                newResult.SetOffsetSpecInfo(m_specInfo);
            }
            else if (novac::Equals(m_specInfo.m_name, "dark_cur"))
            {
                newResult.SetDarkCurrentSpecInfo(m_specInfo);
            }
            else
            {
                newResult.AppendResult(m_evResult, m_specInfo);
                newResult.SetFlux(flux);
                newResult.SetInstrumentType(m_instrumentType);
            }

            if (m_col.peakSaturation != -1)
            {
                // If the intensity is specified as a saturation ratio...
                // double dynamicRange = CSpectrometerModel::GetMaxIntensity(m_specInfo.m_specModel);
            }

            // Guess the spectrometer model
            if (m_spectrometerModel.IsUnknown())
            {
                m_spectrometerModel = novac::CSpectrometerDatabase::GetInstance().GuessModelFromSerial(m_specInfo.m_device);
            }

            newResult.CheckGoodnessOfFit(m_specInfo, &m_spectrometerModel);
            ++measNr;
        }

        // close the evaluation log
        fclose(f);
    }
    singleLock.Unlock();

    // If the sky and dark were specified, remove them from the measurement
    if (fabs(newResult.GetScanAngle(1) - 180.0) < 1)
    {
        newResult.RemoveResult(0); // remove sky
        newResult.RemoveResult(0); // remove dark
    }

    // Insert the new scan
    m_scan.push_back(newResult);

    newResult = Evaluation::CScanResult();

    // Sort the scans in order of collection
    SortScans();

    return RETURN_CODE::SUCCESS;
}

size_t CEvaluationLogFileHandler::CountScansInFile()
{
    char expTimeStr[] = "exposuretime"; // this string only exists in the header line.
    char szLine[8192];
    size_t nScans = 0;

    // If no evaluation log selected, quit
    if (m_evaluationLog.size() <= 1)
    {
        return 0;
    }

    novac::CSingleLock singleLock(&g_evalLogCritSect);
    singleLock.Lock();
    if (singleLock.IsLocked())
    {

        // Open the evaluation log
        FILE* f = fopen(m_evaluationLog.c_str(), "r");
        if (nullptr == f)
        {
            singleLock.Unlock();
            return 0;
        }

        // Read the file, one line at a time
        while (fgets(szLine, 8192, f))
        {
            // convert the string to all lower-case letters
            for (unsigned int it = 0; it < strlen(szLine); ++it)
            {
                szLine[it] = (char)tolower(szLine[it]);
            }

            // find the next start of a scan
            if (nullptr != strstr(szLine, expTimeStr))
            {
                ++nScans;
            }
        }

        fclose(f);
    }
    singleLock.Unlock();

    // Return the number of scans found in the file
    return nScans;
}

/** Reads and parses the 'scanInfo' header before the scan */
void CEvaluationLogFileHandler::ParseScanInformation(novac::CSpectrumInfo& scanInfo, double& flux, FILE* f)
{
    char szLine[8192];
    char* pt = nullptr;
    double tmpDouble;
    char instrumentType[64];

    // Reset the column- and spectrum info
    ResetColumns();

    // read the additional scan-information, line by line
    while (fgets(szLine, 8192, f))
    {

        // convert to lower-case
        for (unsigned int it = 0; it < strlen(szLine); ++it)
        {
            szLine[it] = (char)tolower(szLine[it]);
        }

        pt = strstr(szLine, "</scaninformation>");
        if (nullptr != pt)
        {
            break;
        }

        pt = strstr(szLine, "compiledate=");
        if (nullptr != pt)
        {
            continue;
        }

        pt = strstr(szLine, "site=");
        if (nullptr != pt)
        {
            scanInfo.m_site = std::string(pt + 5);
            Remove(scanInfo.m_site, '\n'); // Remove newline characters
            continue;
        }

        pt = strstr(szLine, "date=");
        if (nullptr != pt)
        {
            int tmpInt[3];
            if (3 == sscanf(pt + 5, "%d.%d.%d", &tmpInt[0], &tmpInt[1], &tmpInt[2]))
            {
                scanInfo.m_startTime.year = (unsigned short)tmpInt[2];
                scanInfo.m_startTime.month = (unsigned char)tmpInt[1];
                scanInfo.m_startTime.day = (unsigned char)tmpInt[0];
                scanInfo.m_stopTime.year = (unsigned short)tmpInt[2];
                scanInfo.m_stopTime.month = (unsigned char)tmpInt[1];
                scanInfo.m_stopTime.day = (unsigned char)tmpInt[0];
            }
            continue;
        }

        pt = strstr(szLine, "starttime=");
        if (nullptr != pt)
        {
            int tmpInt[4];
            if (3 == sscanf(pt + strlen("starttime="), "%d.%d.%d", &tmpInt[0], &tmpInt[1], &tmpInt[2]))
            {
                scanInfo.m_startTime.hour = (unsigned char)tmpInt[0];
                scanInfo.m_startTime.minute = (unsigned char)tmpInt[1];
                scanInfo.m_startTime.second = (unsigned char)tmpInt[2];
            }
            else if (3 == sscanf(pt + strlen("starttime="), "%d:%d:%d", &tmpInt[0], &tmpInt[1], &tmpInt[2]))
            {
                scanInfo.m_startTime.hour = (unsigned char)tmpInt[0];
                scanInfo.m_startTime.minute = (unsigned char)tmpInt[1];
                scanInfo.m_startTime.second = (unsigned char)tmpInt[2];
            }
            continue;
        }

        pt = strstr(szLine, "stoptime=");
        if (nullptr != pt)
        {
            int tmpInt[3];
            if (3 == sscanf(pt + 9, "%d.%d.%d", &tmpInt[0], &tmpInt[1], &tmpInt[2]))
            {
                scanInfo.m_stopTime.hour = (unsigned char)tmpInt[0];
                scanInfo.m_stopTime.minute = (unsigned char)tmpInt[1];
                scanInfo.m_stopTime.second = (unsigned char)tmpInt[2];
            }
            else if (3 == sscanf(pt + 9, "%d:%d:%d", &tmpInt[0], &tmpInt[1], &tmpInt[2]))
            {
                scanInfo.m_stopTime.hour = (unsigned char)tmpInt[0];
                scanInfo.m_stopTime.minute = (unsigned char)tmpInt[1];
                scanInfo.m_stopTime.second = (unsigned char)tmpInt[2];
            }
            continue;
        }

        pt = strstr(szLine, "compass=");
        if (nullptr != pt)
        {
            if (1 == sscanf(pt + 8, "%lf", &tmpDouble))
            {
                scanInfo.m_compass = (float)fmod(tmpDouble, 360.0);
            }
            continue;
        }

        pt = strstr(szLine, "tilt=");
        if (nullptr != pt)
        {
            if (1 == sscanf(pt + 5, "%lf", &tmpDouble))
            {
                scanInfo.m_pitch = (float)tmpDouble;
            }
        }

        pt = strstr(szLine, "lat=");
        if (nullptr != pt)
        {
            if (1 == sscanf(pt + 4, "%lf", &tmpDouble))
            {
                scanInfo.m_gps.m_latitude = tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "long=");
        if (nullptr != pt)
        {
            if (1 == sscanf(pt + 5, "%lf", &tmpDouble))
            {
                scanInfo.m_gps.m_longitude = tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "alt=");
        if (nullptr != pt)
        {
            if (1 == sscanf(pt + 4, "%lf", &tmpDouble))
            {
                scanInfo.m_gps.m_altitude = (long)tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "serial=");
        if (nullptr != pt)
        {
            scanInfo.m_device = std::string(pt + 7);
            Remove(scanInfo.m_device, '\n'); // remove remaining strange things in the serial-number
            MakeUpper(scanInfo.m_device);    // Convert the serial-number to all upper case letters

            // Extract the spectrometer-model from the serial-number of the spectrometer
            if (m_spectrometerModel.IsUnknown())
            {
                m_spectrometerModel = novac::CSpectrometerDatabase::GetInstance().GuessModelFromSerial(scanInfo.m_device);
            }
            scanInfo.m_specModelName = m_spectrometerModel.modelName;

            continue;
        }

        pt = strstr(szLine, "spectrometer=");
        if (nullptr != pt)
        {
            scanInfo.m_specModelName = std::string(pt + 13);
            continue;
        }

        pt = strstr(szLine, "volcano=");
        if (nullptr != pt)
        {
            scanInfo.m_volcano = std::string(pt + 8);
            Remove(scanInfo.m_volcano, '\n'); // Remove newline characters
            continue;
        }

        pt = strstr(szLine, "observatory=");
        if (nullptr != pt)
        {
            scanInfo.m_observatory = std::string(pt + 12);
            Remove(scanInfo.m_observatory, '\n'); // Remove newline characters
            continue;
        }

        pt = strstr(szLine, "channel=");
        if (nullptr != pt)
        {
            if (1 == sscanf(pt + 8, "%lf", &tmpDouble))
            {
                scanInfo.m_channel = (unsigned char)tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "coneangle=");
        if (nullptr != pt)
        {
            if (1 == sscanf(pt + 10, "%lf", &tmpDouble))
            {
                scanInfo.m_coneAngle = (float)tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "flux=");
        if (nullptr != pt)
        {
            if (1 == sscanf(pt + 5, "%lf", &tmpDouble))
            {
                flux = tmpDouble;
            }
            continue;
        }

        pt = strstr(szLine, "battery=");
        if (nullptr != pt)
        {
            float tmpFloat = 0.0F;
            if (1 == sscanf(pt + 8, "%f", &tmpFloat))
            {
                scanInfo.m_batteryVoltage = tmpFloat;
            }
            continue;
        }

        pt = strstr(szLine, "temperature");
        if (nullptr != pt)
        {
            float tmpFloat = 0.0F;
            if (1 == sscanf(pt + 12, "%f", &tmpFloat))
            {
                scanInfo.m_temperature = tmpFloat;
            }
            continue;
        }

        pt = strstr(szLine, "instrumenttype=");
        if (nullptr != pt)
        {
            sscanf(pt + 15, "%s", instrumentType);
            if (novac::Equals(instrumentType, "heidelberg"))
            {
                m_instrumentType = novac::NovacInstrumentType::Heidelberg;
            }
            else
            {
                m_instrumentType = novac::NovacInstrumentType::Gothenburg;
            }
        }
    }
}

void CEvaluationLogFileHandler::ParseFluxInformation(Meteorology::WindField& windField, double& flux, FILE* f)
{
    char szLine[8192];
    char* pt = nullptr;
    double windSpeed = 10, windDirection = 0, plumeHeight = 1000;
    Meteorology::MeteorologySource windSpeedSource = Meteorology::MeteorologySource::User;
    Meteorology::MeteorologySource windDirectionSource = Meteorology::MeteorologySource::User;
    // Meteorology::MeteorologySource plumeHeightSource = Meteorology::MeteorologySource::User;
    char source[512];

    // read the additional scan-information, line by line
    while (fgets(szLine, 8192, f))
    {
        pt = strstr(szLine, "</fluxinfo>");
        if (nullptr != pt)
        {
            // save all the values
            //			windField.SetPlumeHeight(plumeHeight, plumeHeightSource);
            windField.SetWindDirection(windDirection, windDirectionSource);
            windField.SetWindSpeed(windSpeed, windSpeedSource);
            break;
        }

        pt = strstr(szLine, "flux=");
        if (nullptr != pt)
        {
            sscanf(pt + 5, "%lf", &flux);
            continue;
        }

        pt = strstr(szLine, "windspeed=");
        if (nullptr != pt)
        {
            sscanf(pt + 10, "%lf", &windSpeed);
            continue;
        }

        pt = strstr(szLine, "winddirection=");
        if (nullptr != pt)
        {
            sscanf(pt + 14, "%lf", &windDirection);
            continue;
        }

        pt = strstr(szLine, "plumeheight=");
        if (nullptr != pt)
        {
            sscanf(pt + 12, "%lf", &plumeHeight);
            continue;
        }

        pt = strstr(szLine, "windspeedsource=");
        if (nullptr != pt)
        {
            sscanf(pt + 16, "%s", source);
            windSpeedSource = Meteorology::StringToMetSource(source);
            continue;
        }

        pt = strstr(szLine, "winddirectionsource=");
        if (nullptr != pt)
        {
            sscanf(pt + 20, "%s", source);
            windDirectionSource = Meteorology::StringToMetSource(source);
            continue;
        }

        /* pt = strstr(szLine, "plumeheightsource=");
        if (nullptr != pt)
        {
            sscanf(pt + 18, "%s", source);
            plumeHeightSource = Meteorology::StringToMetSource(source);
            continue;
            } */
    }
}

void CEvaluationLogFileHandler::ResetColumns()
{
    for (int k = 0; k < MAX_N_REFERENCES; ++k)
    {
        m_col.column[k] = -1;
        m_col.columnError[k] = -1;
        m_col.shift[k] = -1;
        m_col.shiftError[k] = -1;
        m_col.squeeze[k] = -1;
        m_col.squeezeError[k] = -1;
    }
    m_col.delta = m_col.intensity = m_col.position = m_col.position2 = -1;
    m_col.nSpecies = 0;
    m_evResult.m_referenceResult.clear();
    m_col.expTime = m_col.nSpec = -1;
    m_col.name = -1;
}

void CEvaluationLogFileHandler::ResetScanInformation()
{
    m_specInfo.m_channel = 0;
    m_specInfo.m_compass = 0.0;
    m_specInfo.m_scanAngle = 0.0;
    m_col.starttime = -1;
    m_col.stoptime = -1;
}

void CEvaluationLogFileHandler::SortScans()
{
    // If the scans are already in order then we don't need to sort them
    //	or if there is only one scan then we don't need to fix it.
    if (IsSorted() || m_scan.size() <= 1)
    {
        return;
    }

    // Then sort the array
    std::sort(begin(m_scan), end(m_scan), [&](const Evaluation::CScanResult& r1, const Evaluation::CScanResult& r2)
              { return r1.GetSkyStartTime() < r2.GetSkyStartTime(); });

    //	CEvaluationLogFileHandler::SortScans(m_scan);
    // CEvaluationLogFileHandler::BubbleSortScans(m_scan);

    return;
}

/** Returns true if the scans are already ordered */
bool CEvaluationLogFileHandler::IsSorted()
{
    novac::CDateTime time1, time2;

    for (size_t k = 0; k < m_scan.size() - 1; ++k)
    {
        // Get the start-times
        m_scan[k].GetStartTime(0, time1);
        m_scan[k + 1].GetStartTime(0, time2);

        // If the second scan has started before the first,
        //	then change the order, otherwise don't do anything
        if (time2 < time1)
        {
            return false;
        }
        else
        {
            continue;
        }
    }

    return true; // no un-ordered scans were found
}

RETURN_CODE CEvaluationLogFileHandler::WriteEvaluationLog(const std::string& fileName, novac::SpectrometerModel spectrometerModel, int softwareMajorNumber, int softwareMinorNumber)
{
    novac::CString string;
    novac::CString wsSrc, wdSrc, phSrc;
    novac::CDateTime startTime;

    // 1. Test if the file already exists, if so then return false
    if (Filesystem::IsExistingFile(fileName))
        return RETURN_CODE::FAIL;

    // 2. Write the file
    FILE* f = fopen(fileName.c_str(), "w");

    for (size_t scanIndex = 0; scanIndex < this->m_scan.size(); ++scanIndex)
    {
        Evaluation::CScanResult& scan = this->m_scan[scanIndex];
        Meteorology::WindField& wind = this->m_windField[scanIndex];

        scan.GetStartTime(0, startTime);

        // ----------------- Create the additional scan-information ----------------------
        string.Format("\n<scaninformation>\n");
        string.AppendFormat("\tdate=%02d.%02d.%04d\n", startTime.day, startTime.month, startTime.year);
        string.AppendFormat("\tstarttime=%02d:%02d:%02d\n", startTime.hour, startTime.minute, startTime.second);
        string.AppendFormat("\tcompass=%.1lf\n", scan.GetCompass());
        string.AppendFormat("\ttilt=%.1lf\n", scan.GetPitch());
        string.AppendFormat("\tlat=%.6lf\n", scan.GetLatitude());
        string.AppendFormat("\tlong=%.6lf\n", scan.GetLongitude());
        string.AppendFormat("\talt=%ld\n", scan.GetAltitude());

        string.AppendFormat("\tvolcano=%s\n", m_specInfo.m_volcano.c_str());
        string.AppendFormat("\tsite=%s\n", m_specInfo.m_site.c_str());
        string.AppendFormat("\tobservatory=%s\n", m_specInfo.m_observatory.c_str());

        string.AppendFormat("\tserial=%s\n", scan.GetSerial().c_str());

        string.AppendFormat("\tspectrometer=%s\n", m_specInfo.m_specModelName.c_str());

        if (!spectrometerModel.IsUnknown())
        {
            string.AppendFormat("\tspectrometer_maxIntensity=%lf\n", spectrometerModel.maximumIntensityForSingleReadout);
            string.AppendFormat("\tspectrometer_numPixels=%d\n", spectrometerModel.numberOfPixels);
        }

        string.AppendFormat("\tchannel=%d\n", m_specInfo.m_channel);
        string.AppendFormat("\tconeangle=%.1lf\n", scan.GetConeAngle());
        string.AppendFormat("\tinterlacesteps=%d\n", m_specInfo.m_interlaceStep);
        string.AppendFormat("\tstartchannel=%d\n", m_specInfo.m_startChannel);
        string.AppendFormat("\tspectrumlength=%d\n", 2048);
        string.AppendFormat("\tflux=%.2lf\n", scan.GetFlux());
        string.AppendFormat("\tbattery=%.2f\n", scan.GetBatteryVoltage());
        string.AppendFormat("\ttemperature=%.2f\n", scan.GetTemperature());

        // The mode
        string.AppendFormat("\tmode=%s\n", novac::ToString(scan.m_measurementMode).c_str());

        // The type of instrument used...
        string.AppendFormat("\tinstrumenttype=%s\n", novac::ToString(scan.m_instrumentType).c_str());

        // Finally, the version of the file and the version of the program
        string.Append("\tversion=2.0\n");
        string.AppendFormat("\tsoftwareversion=%d.%02d\n", softwareMajorNumber, softwareMinorNumber);
        string.AppendFormat("\tcompiledate=%s\n", __DATE__);

        string.Append("</scaninformation>\n\n");
        fprintf(f, "%s", string.c_str());

        // ----------------- Create the flux-information ----------------------
        wind.GetWindSpeedSource(wsSrc);
        wind.GetWindDirectionSource(wdSrc);
        //		wind.GetPlumeHeightSource(phSrc);
        string.Format("<fluxinfo>\n");
        string.AppendFormat("\tflux=%.4lf\n", scan.GetFlux());
        string.AppendFormat("\twindspeed=%.4lf\n", wind.GetWindSpeed());
        string.AppendFormat("\twinddirection=%.4lf\n", wind.GetWindDirection());
        //		string.AppendFormat("\tplumeheight=%.2lf\n",			wind.GetPlumeHeight());
        string.AppendFormat("\twindspeedsource=%s\n", (const char*)wsSrc);
        string.AppendFormat("\twinddirectionsource=%s\n", (const char*)wdSrc);
        string.AppendFormat("\tplumeheightsource=%s\n", (const char*)phSrc);
        // if(fabs(spectrometer.m_scanner.compass) > 360.0)
        //	string.AppendFormat("\tcompasssource=compassreading\n");
        // else
        //	string.AppendFormat("\tcompasssource=user\n");
        string.Append("</fluxinfo>\n");
        fprintf(f, "%s", string.c_str());

        // ----------------------- write the header --------------------------------
        if (m_instrumentType == novac::NovacInstrumentType::Gothenburg)
        {
            string.Format("#scanangle\t");
        }
        else if (m_instrumentType == novac::NovacInstrumentType::Heidelberg)
        {
            string.Format("#observationangle\tazimuth\t");
        }
        string.Append("starttime\tstoptime\tname\tspecsaturation\tfitsaturation\tdelta\tchisquare\texposuretime\tnumspec\t");

        for (size_t itSpecie = 0; itSpecie < scan.GetSpecieNum(0); ++itSpecie)
        {
            const std::string specieName = scan.GetSpecieName(0, itSpecie);
            string.AppendFormat("column(%s)\tcolumnerror(%s)\t", specieName.c_str(), specieName.c_str());
            string.AppendFormat("shift(%s)\tshifterror(%s)\t", specieName.c_str(), specieName.c_str());
            string.AppendFormat("squeeze(%s)\tsqueezeerror(%s)\t", specieName.c_str(), specieName.c_str());
        }
        string.Append("isgoodpoint\toffset\tflag");
        string.Append("\n<spectraldata>\n");

        fprintf(f, "%s", string.c_str());

        // ------------------- Then write the parameters for each spectrum ---------------------------
        for (size_t itSpectrum = 0; itSpectrum < scan.GetEvaluatedNum(); ++itSpectrum)
        {
            // 3a. Pretty print the result and the spectral info into a string
            FormatEvaluationResult(&scan.GetSpectrumInfo(itSpectrum), scan.GetResult(itSpectrum), m_instrumentType, 0.0, scan.GetSpecieNum(itSpectrum), string);

            // 3b. Write it to the evaluation log file
            fprintf(f, "%s", string.c_str());
            fprintf(f, "\n");
        }
        fprintf(f, "</spectraldata>\n");
    }

    // Remember to close the file
    fclose(f);

    return RETURN_CODE::SUCCESS;
}

RETURN_CODE CEvaluationLogFileHandler::FormatEvaluationResult(const novac::CSpectrumInfo* info, const novac::CEvaluationResult* result, novac::NovacInstrumentType iType, double maxIntensity, size_t nSpecies, novac::CString& string)
{
    if (result != nullptr && result->m_referenceResult.size() < nSpecies)
    {
        return RETURN_CODE::FAIL; // something's wrong here!
    }

    // 1. The Scan angle
    string.Format("%.0lf\t", info->m_scanAngle);

    // 2. The azimuth angle
    if (iType == novac::NovacInstrumentType::Heidelberg)
    {
        string.AppendFormat("%.0lf\t", info->m_scanAngle2);
    }

    // 3. The start time
    string.AppendFormat("%02d:%02d:%02d\t", info->m_startTime.hour, info->m_startTime.minute, info->m_startTime.second);

    // 4. The stop time
    string.AppendFormat("%02d:%02d:%02d\t", info->m_stopTime.hour, info->m_stopTime.minute, info->m_stopTime.second);

    // 5 The name of the spectrum
    novac::CString simpleName;
    SimplifyString(info->m_name, simpleName);
    string.AppendFormat("%s\t", (const char*)simpleName);

    // 6. The (maximum) saturation ratio of the whole spectrum,
    //      the (maximum) saturation ratio in the fit-region
    //      and the normalized maximum intensity of the whole spectrum
    if (maxIntensity > 0.0)
    {
        string.AppendFormat("%.2lf\t", info->m_peakIntensity / maxIntensity);
        string.AppendFormat("%.2lf\t", info->m_fitIntensity / maxIntensity);
    }
    else
    {
        string.AppendFormat("%.2lf\t", info->m_peakIntensity);
        string.AppendFormat("%.2lf\t", info->m_fitIntensity);
    }
    string.AppendFormat("%.2lf\t", (info->m_peakIntensity - info->m_offset) / info->m_exposureTime);

    // 7. The delta of the fit
    if (result != nullptr)
        string.AppendFormat("%.2e\t", result->m_delta);
    else
        string.AppendFormat("%.2e\t", 0.0);

    // 8. The chi-square of the fit
    if (result != nullptr)
        string.AppendFormat("%.2e\t", result->m_chiSquare);
    else
        string.AppendFormat("%.2e\t", 0.0);

    // 9. The exposure time and the number of spectra averaged
    string.AppendFormat("%ld\t%ld\t", info->m_exposureTime, info->m_numSpec);

    // 10. The column/column error for each specie
    for (size_t itSpecie = 0; itSpecie < nSpecies; ++itSpecie)
    {
        if (result != nullptr)
        {
            string.AppendFormat("%.2e\t%.2e\t", result->m_referenceResult[itSpecie].m_column, result->m_referenceResult[itSpecie].m_columnError);
            string.AppendFormat("%.2lf\t%.2lf\t", result->m_referenceResult[itSpecie].m_shift, result->m_referenceResult[itSpecie].m_shiftError);
            string.AppendFormat("%.2lf\t%.2lf\t", result->m_referenceResult[itSpecie].m_squeeze, result->m_referenceResult[itSpecie].m_squeezeError);
        }
        else
        {
            string.Append("0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t");
        }
    }

    // 11. The quality of the fit
    if (result != nullptr)
        string.AppendFormat("%d\t", result->IsOK());
    else
        string.Append("1\t");

    // 12. The offset
    string.AppendFormat("%.0lf\t", info->m_offset);

    // 13. The 'flag' in the spectra
    string.AppendFormat("%d", info->m_flag);

    return RETURN_CODE::SUCCESS;
}
