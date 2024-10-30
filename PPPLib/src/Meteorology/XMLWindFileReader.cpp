#include <PPPLib/Meteorology/XMLWindFileReader.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/Communication/FTPServerConnection.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <PPPLib/MFC/CList.h>
#include <PPPLib/Logging.h>
#include <SpectralEvaluation/Exceptions.h>
#include <Poco/Glob.h>
#include <Poco/Path.h>
#include <string.h>
#include <memory>
#include <cmath>

#ifndef MAX_PATH
#define MAX_PATH 260
#endif

using namespace FileHandler;
using namespace Meteorology;
using namespace novac;

CXMLWindFileReader::CXMLWindFileReader(ILogger& logger, const Configuration::CUserConfiguration& userSettings)
    : CXMLFileReader(logger), m_userSettings(userSettings)
{
}

CXMLWindFileReader::~CXMLWindFileReader(void)
{
}

void CXMLWindFileReader::ReadWindFile(novac::LogContext context, const novac::CString& fileName, Meteorology::CWindDataBase& dataBase)
{
    novac::CString localFileName, userMessage;

    // 0. If the file is on the server, then download it first
    if (Equals(fileName.Left(6), "ftp://"))
    {
        Communication::CFTPServerConnection ftp(m_log, m_userSettings);

        novac::CString tmpFileName;
        tmpFileName.Format(fileName);
        novac::CFileUtils::GetFileName(tmpFileName); // this is the name of the file, without the path...
        localFileName.Format("%s%s", (const char*)m_userSettings.m_tempDirectory, (const char*)tmpFileName);

        // make sure that the tmp-directory exists
        if (Filesystem::CreateDirectoryStructure(m_userSettings.m_tempDirectory))
        {
            userMessage.Format("Could not create temp directory: '%s' required for downloading wind field from ftp", (const char*)m_userSettings.m_tempDirectory);
            throw novac::FileIoException(userMessage.c_str());
        }

        if (ftp.DownloadFileFromFTP(context, fileName, localFileName, m_userSettings.m_FTPUsername, m_userSettings.m_FTPPassword))
        {
            throw novac::FileIoException("Failed to download wind file from FTP server");
        }
    }
    else
    {
        localFileName.Format("%s", (const char*)fileName);
    }

    // 1. Open the file
    if (!Open(localFileName))
    {
        std::string message{ std::string("Failed to open wind field file for reading: '") + localFileName.std_str() };
        throw novac::FileIoException(message.c_str());
    }

    // parse the file
    while (nullptr != (szToken = NextToken()))
    {
        // no use to parse empty lines
        if (strlen(szToken) < 3)
        {
            continue;
        }

        novac::CString trimmedToken(szToken);
        trimmedToken.Trim(); // remove blanks in the beginning and in the end

        // If this is a wind-field item then parse this. 
        if (novac::Equals(trimmedToken, "windfield", 9))
        {
            Parse_WindField(dataBase);
            continue;
        }

        // if this is the beginning of the wind-section,  extract the name of the 
        //	database
        if (novac::Equals(trimmedToken, "Wind", 4))
        {
            const char* volcanoName = GetAttributeValue("volcano");
            if (volcanoName != nullptr)
            {
                dataBase.m_dataBaseName = std::string(volcanoName);
            }
            continue;
        }
    }

    Close();
}

void CXMLWindFileReader::ReadWindDirectory(novac::LogContext context, const novac::CString& directory, Meteorology::CWindDataBase& dataBase, const CDateTime* dateFrom, const CDateTime* dateTo)
{
    novac::CStdioFile file;
    novac::CString localFileName, remoteFileName, userMessage, ftpDir;
    std::vector<std::string> remoteFileList; // the list of wind field files
    novac::CList <novac::CString, novac::CString&> localFileList; // the list of files on the local computer

    if (Equals(directory.Left(6), "ftp://"))
    {
        // If the directory is on the server, then this is how to check the files
        std::unique_ptr<Communication::CFTPServerConnection> ftp;
        ftp.reset(new Communication::CFTPServerConnection(m_log, m_userSettings));

        // Make sure that the directory does end with a trailing '/'
        ftpDir.Format(directory);
        if (!Equals(ftpDir.Right(1), "/"))
        {
            ftpDir.Append("/");
        }

        // Get the list of files on the server
        if (ftp->DownloadFileListFromFTP(ftpDir, remoteFileList, m_userSettings.m_FTPUsername, m_userSettings.m_FTPPassword))
        {
            throw novac::FileIoException("Failed to download list of wind files from FTP server");
        }

        // make sure that the tmp-directory exists
        if (Filesystem::CreateDirectoryStructure(m_userSettings.m_tempDirectory))
        {
            userMessage.Format("Could not create temp directory: '%s' required to download wind fields from FTP server", (const char*)m_userSettings.m_tempDirectory);
            throw novac::FileIoException(userMessage.c_str());
        }

        // Download the files, one at a time
        for (std::string item : remoteFileList)
        {
            novac::CString name(item);

            // only download .wxml files
            if (!Equals(name.Right(5), ".wxml"))
            {
                continue;
            }

            // if this file has a name that matches XXXX_YYYYMMDD.wxml then use this info to see 
            //	weather we actually should download this file
            int rpos = name.ReverseFind('_');
            if (rpos > 0 && ((name.GetLength() - rpos) == 14))
            {
                novac::CString dateStr(name.Right(13).Left(8));
                CDateTime t;
                if (CDateTime::ParseDate(dateStr, t))
                {
                    if (t < m_userSettings.m_fromDate || t > m_userSettings.m_toDate)
                    {
                        continue;
                    }
                }
            }

            localFileName.Format("%s%s", (const char*)m_userSettings.m_tempDirectory, (const char*)name);

            if (Filesystem::IsExistingFile(localFileName))
            {
                userMessage.Format("File %s is already downloaded", (const char*)localFileName);
                m_log.Information(userMessage.std_str());
            }
            else
            {
                remoteFileName.Format("%s%s", (const char*)ftpDir, (const char*)name);

                if (ftp->DownloadFileFromFTP(context, remoteFileName, localFileName, m_userSettings.m_FTPUsername, m_userSettings.m_FTPPassword))
                {
                    throw novac::FileIoException("Failed to download wind file from FTP server");
                }
            }

            localFileList.AddTail(localFileName);
        }
    }
    else
    {
        // If the directory is on the local computer, then this is how to check the files
        char fileToFind[MAX_PATH];
        sprintf(fileToFind, "%s/*.wxml", (const char*)directory);
        m_log.Information(context, "Searching for wind files in: " + directory.std_str());

        // Search for the files
        std::set<std::string> filesFound;
        Poco::Glob::glob(fileToFind, filesFound);

        if (filesFound.size() == 0)
        {
            throw std::invalid_argument("No wind files found");
        }

        for (const std::string& fName : filesFound)
        {
            localFileList.AddTail(fName);
        }
    }

    // Now we got a list of files on the local computer. Read them in!
    auto p = localFileList.GetHeadPosition();
    int nFilesRead = 0;
    while (p != nullptr)
    {
        localFileName.Format("%s", (const char*)localFileList.GetNext(p));

        // make sure that this file falls in the appropriate date-range
        if (dateFrom != nullptr)
        {

        }
        if (dateTo != nullptr)
        {

        }

        try
        {
            ReadWindFile(context, localFileName, dataBase);
        }
        catch (const std::exception& e)
        {
            ShowMessage(e.what());
        }
    }

    // Tell the user what we've done
    if (nFilesRead > 0)
    {
        userMessage.Format("Successfully read in %d wind field files", nFilesRead);
        m_log.Information(context, userMessage.std_str());
        return;
    }

    throw std::invalid_argument("Failed to find any wind field files in the specified directory");;
}

int CXMLWindFileReader::Parse_WindField(Meteorology::CWindDataBase& dataBase)
{
    novac::CString sourceStr, userMessage;
    Meteorology::WindField w;
    CDateTime validFrom, validTo;
    double latitude = 0.0;
    double longitude = 0.0;
    double altitude = 0.0;
    double windspeed = 0.0, windspeederror = 0.0;
    double winddirection = 0.0, winddirectionerror = 0.0;
    MeteorologySource windSource = MeteorologySource::Default;

    // parse the file
    while (nullptr != (szToken = NextToken()))
    {
        // no use to parse empty lines
        if (strlen(szToken) < 2)
        {
            continue;
        }

        // ignore comments
        if (Equals(szToken, "!--", 3))
        {
            continue;
        }

        // end of fit-window section
        if (IsClosingTag("/windfield", szToken))
        {
            return 0;
        }

        if (Equals(szToken, "source"))
        {
            Parse_StringItem("/source", sourceStr);
            windSource = Meteorology::StringToMetSource(sourceStr);
            continue;
        }

        if (Equals(szToken, "altitude"))
        {
            Parse_FloatItem("/altitude", altitude);
            continue;
        }

        if (Equals(szToken, "valid_from"))
        {
            Parse_Date("/valid_from", validFrom);
            continue;
        }

        if (Equals(szToken, "valid_to"))
        {
            Parse_Date("/valid_to", validTo);
            continue;
        }

        if (Equals(szToken, "item", 4))
        {
            // find the wind-speed
            const char* str = GetAttributeValue("ws");
            if (str != nullptr)
            {
                windspeed = atof(str);
            }

            // find the wind-speed-error
            str = GetAttributeValue("wse");
            if (str != nullptr)
            {
                windspeederror = atof(str);
            }

            // find the wind-direction
            str = GetAttributeValue("wd");
            if (str != nullptr)
            {
                winddirection = atof(str);
            }
            // find the wind-direction-error
            str = GetAttributeValue("wde");
            if (str != nullptr)
            {
                winddirectionerror = atof(str);
            }

            // find the latitude
            str = GetAttributeValue("lat");
            if (str != nullptr)
            {
                latitude = atof(str);
            }

            // find the longitude
            str = GetAttributeValue("lon");
            if (str != nullptr)
            {
                longitude = atof(str);
            }
            else
            {
                str = GetAttributeValue("long");
                if (str != nullptr)
                {
                    longitude = atof(str);
                }
            }

            // check that the latitude is within -90 to +90 degrees...
            latitude = (latitude > 90.0) ? latitude - std::floor(latitude / 90.0) * 90.0 : latitude;

            // check that the longitude is within -180 to +180 degrees...
            longitude = (longitude > 180.0) ? longitude - (1 + std::floor(longitude / 360.0)) * 360.0 : longitude;

            // check the reasonability of the values
            if (winddirection < -360.0 || winddirection > 360.0)
            {
                userMessage.Format("Received wind-field file with invalid wind direction (%lf degrees) in file %s", winddirection, (const char*)m_filename);
                ShowMessage(userMessage);
            }
            if (windspeed < 0.0 || windspeed > 50.0)
            {
                userMessage.Format("Received wind-field file with invalid wind speed (%lf m/s) in file %s", windspeed, (const char*)m_filename);
                ShowMessage(userMessage);
            }

            // we have now enough information to make a wind-field and insert it into the database
            w = WindField(windspeed, windspeederror, windSource, winddirection, winddirectionerror, windSource, validFrom, validTo, latitude, longitude, altitude);

            dataBase.InsertWindField(w);
        }
    }

    return 1;
}

int CXMLWindFileReader::WriteWindFile(const novac::CString& fileName, const Meteorology::CWindDataBase& dataBase)
{
    return dataBase.WriteToFile(fileName);
}
