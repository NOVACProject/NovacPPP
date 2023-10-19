#include <PPPLib/Meteorology/XMLWindFileReader.h>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

// we need to be able to download data from the FTP-server
#include <PPPLib/Communication/FTPServerConnection.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <PPPLib/MFC/CList.h>
#include <PPPLib/Logging.h>
#include <Poco/Glob.h>
#include <Poco/Path.h>
#include <string.h>
#include <memory>

extern Configuration::CUserConfiguration g_userSettings;// <-- The settings of the user

using namespace FileHandler;
using namespace Meteorology;
using namespace novac;

CXMLWindFileReader::CXMLWindFileReader(ILogger& logger)
    : CXMLFileReader(logger)
{
}

CXMLWindFileReader::~CXMLWindFileReader(void)
{
}

/** Reads in an wind-field file.
    In the format specified for the NovacPostProcessingProgram (NPPP)
    @return 0 on sucess */
int CXMLWindFileReader::ReadWindFile(const novac::CString& fileName, Meteorology::CWindDataBase& dataBase) {
    novac::CString localFileName, userMessage;

    // 0. If the file is on the server, then download it first
    if (Equals(fileName.Left(6), "ftp://")) {
        Communication::CFTPServerConnection ftp;

        novac::CString tmpFileName;
        tmpFileName.Format(fileName);
        novac::CFileUtils::GetFileName(tmpFileName); // this is the name of the file, without the path...
        localFileName.Format("%s%c%s", (const char*)g_userSettings.m_tempDirectory, Poco::Path::separator(), (const char*)tmpFileName);

        // make sure that the tmp-directory exists
        if (Filesystem::CreateDirectoryStructure(g_userSettings.m_tempDirectory)) {
            userMessage.Format("Could not create temp directory: %s", (const char*)g_userSettings.m_tempDirectory);
            ShowMessage(userMessage);
            return 1;
        }

        if (ftp.DownloadFileFromFTP(fileName, localFileName, g_userSettings.m_FTPUsername, g_userSettings.m_FTPPassword)) {
            ShowMessage("Failed to download wind file from FTP server");
            return 1;
        }
    }
    else {
        localFileName.Format("%s", (const char*)fileName);
    }


    // 1. Open the file
    if (!Open(localFileName)) {
        ShowMessage(std::string("Failed to open wind field file for reading: '") + localFileName.std_str());
        return 1;
    }

    // parse the file
    while (nullptr != (szToken = NextToken())) {

        // no use to parse empty lines
        if (strlen(szToken) < 3) {
            continue;
        }

        // If this is a wind-field item then parse this. 
        if (novac::Equals(szToken, "windfield", 9)) {
            Parse_WindField(dataBase);
            continue;
        }

        // if this is the beginning of the wind-section,  extract the name of the 
        //	database
        if (novac::Equals(szToken, "Wind", 4)) {
            const char* volcanoName = GetAttributeValue("volcano");
            if (volcanoName != NULL) {
                dataBase.m_dataBaseName.Format("%s", volcanoName);
            }
            continue;
        }
    }

    Close();

    return 0;
}

/** Reads in all the wind-field files that are found in a given directory
    The directory can be on the local computer or on the FTP-server
    @param directory - the full path to the directory where the files are
    @param dataBase - this will on successfull return be filled with the wind
        information found in the wind field files
    @param dateFrom - if not null then only files which contain a wind field after (and including)
        the date 'dateFrom' will be read in.
    @param dateTo - if not null then only file which contain a wind field before (and including)
        the date 'dateTo' will be read in.
    @return 0 on success */
int CXMLWindFileReader::ReadWindDirectory(const novac::CString& directory, Meteorology::CWindDataBase& dataBase, const CDateTime* dateFrom, const CDateTime* dateTo) {
    novac::CStdioFile file;
    novac::CString localFileName, remoteFileName, userMessage, ftpDir;
    std::vector<std::string> remoteFileList; // the list of wind field files
    novac::CList <novac::CString, novac::CString&> localFileList; // the list of files on the local computer

    if (Equals(directory.Left(6), "ftp://")) {
        // If the directory is on the server, then this is how to check the files
        std::unique_ptr<Communication::CFTPServerConnection> ftp;
        ftp.reset(new Communication::CFTPServerConnection());

        // Make sure that the directory does end with a trailing '/'
        ftpDir.Format(directory);
        if (!Equals(ftpDir.Right(1), "/"))
            ftpDir.Append("/");

        // Get the list of files on the server
        if (ftp->DownloadFileListFromFTP(ftpDir, remoteFileList, g_userSettings.m_FTPUsername, g_userSettings.m_FTPPassword)) {
            ShowMessage("Failed to download list of wind files from FTP server");
            return 1;
        }

        // make sure that the tmp-directory exists
        if (Filesystem::CreateDirectoryStructure(g_userSettings.m_tempDirectory)) {
            userMessage.Format("Could not create temp directory: %s", (const char*)g_userSettings.m_tempDirectory);
            ShowMessage(userMessage);
            return 1;
        }

        // Download the files, one at a time
        for (std::string item : remoteFileList)
        {
            novac::CString name(item);

            // only download .wxml files
            if (!Equals(name.Right(5), ".wxml"))
                continue;

            // if this file has a name that matches XXXX_YYYYMMDD.wxml then use this info to see 
            //	weather we actually should download this file
            int rpos = name.ReverseFind('_');
            if (rpos > 0 && ((name.GetLength() - rpos) == 14)) {
                novac::CString dateStr(name.Right(13).Left(8));
                CDateTime t;
                if (CDateTime::ParseDate(dateStr, t)) {
                    if (t < g_userSettings.m_fromDate || t > g_userSettings.m_toDate)
                        continue;
                }
            }

            localFileName.Format("%s%s", (const char*)g_userSettings.m_tempDirectory, (const char*)name);

            if (Filesystem::IsExistingFile(localFileName)) {
                userMessage.Format("File %s is already downloaded", (const char*)localFileName);
                ShowMessage(userMessage);
            }
            else {
                remoteFileName.Format("%s%s", (const char*)ftpDir, (const char*)name);

                if (ftp->DownloadFileFromFTP(remoteFileName, localFileName, g_userSettings.m_FTPUsername, g_userSettings.m_FTPPassword)) {
                    ShowMessage("Failed to download wind file from FTP server");
                    return 1;
                }
            }
            localFileList.AddTail(localFileName);
        }
    }
    else {
        // If the directory is on the local computer, then this is how to check the files
        char fileToFind[MAX_PATH];
        sprintf(fileToFind, "%s/*.wxml", (const char*)directory);

        // Search for the files
        std::set<std::string> filesFound;
        Poco::Glob::glob(fileToFind, filesFound);

        if (filesFound.size() == 0)
        {
            return 1;
        }

        for (const std::string& fName : filesFound) {
            // localFileName.Format("%s%c%s", (const char*)directory, Poco::Path::separator(), fName.c_str());
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
        if (dateFrom != NULL) {

        }
        if (dateTo != NULL) {

        }

        if (0 == ReadWindFile(localFileName, dataBase))
            ++nFilesRead;
    }

    // Tell the user what we've done
    if (nFilesRead > 0) {
        userMessage.Format("Successfully read in %d wind field files", nFilesRead);
        ShowMessage(userMessage);
        return 0;
    }
    else {
        ShowMessage("Failed to find any wind field files in the specified directory");
        return 1;
    }
}

/** Reads a 'windfield' section */
int CXMLWindFileReader::Parse_WindField(Meteorology::CWindDataBase& dataBase) {
    novac::CString sourceStr, userMessage;
    Meteorology::CWindField w;
    CDateTime validFrom, validTo;
    double latitude = 0.0;
    double longitude = 0.0;
    double altitude = 0.0;
    double windspeed = 0.0, windspeederror = 0.0;
    double winddirection = 0.0, winddirectionerror = 0.0;
    MET_SOURCE windSource = MET_DEFAULT;

    // parse the file
    while (nullptr != (szToken = NextToken())) {
        // no use to parse empty lines
        if (strlen(szToken) < 2) {
            continue;
        }

        // ignore comments
        if (Equals(szToken, "!--", 3)) {
            continue;
        }

        // end of fit-window section
        if (Equals(szToken, "/windfield")) {
            return 0;
        }

        if (Equals(szToken, "source")) {
            Parse_StringItem("/source", sourceStr);
            windSource = Meteorology::StringToMetSource(sourceStr);
            continue;
        }

        if (Equals(szToken, "altitude")) {
            Parse_FloatItem("/altitude", altitude);
            continue;
        }

        if (Equals(szToken, "valid_from")) {
            Parse_Date("/valid_from", validFrom);
            continue;
        }

        if (Equals(szToken, "valid_to")) {
            Parse_Date("/valid_to", validTo);
            continue;
        }

        if (Equals(szToken, "item", 4)) {
            // find the wind-speed
            const char* str = GetAttributeValue("ws");
            if (str != NULL)
                windspeed = atof(str);

            // find the wind-speed-error
            str = GetAttributeValue("wse");
            if (str != NULL)
                windspeederror = atof(str);

            // find the wind-direction
            str = GetAttributeValue("wd");
            if (str != NULL)
                winddirection = atof(str);

            // find the wind-direction-error
            str = GetAttributeValue("wde");
            if (str != NULL)
                winddirectionerror = atof(str);

            // find the latitude
            str = GetAttributeValue("lat");
            if (str != NULL)
                latitude = atof(str);

            // find the longitude
            str = GetAttributeValue("lon");
            if (str != NULL)
                longitude = atof(str);
            else {
                str = GetAttributeValue("long");
                if (str != NULL)
                    longitude = atof(str);
            }

            // check that the latitude is within -90 to +90 degrees...
            latitude = (latitude > 90.0) ? latitude - floor(latitude / 90.0) * 90.0 : latitude;

            // check that the longitude is within -180 to +180 degrees...
            longitude = (longitude > 180.0) ? longitude - (1 + floor(longitude / 360.0)) * 360.0 : longitude;

            // check the reasonability of the values
            if (winddirection < -360.0 || winddirection > 360.0) {
                userMessage.Format("Received wind-field file with invalid wind direction (%lf degrees) in file %s", winddirection, (const char*)m_filename);
                ShowMessage(userMessage);
            }
            if (windspeed < 0.0 || windspeed > 50.0) {
                userMessage.Format("Received wind-field file with invalid wind speed (%lf m/s) in file %s", windspeed, (const char*)m_filename);
                ShowMessage(userMessage);
            }

            // we have now enough information to make a wind-field and insert it into the database
            w = CWindField(windspeed, windspeederror, windSource, winddirection, winddirectionerror, windSource, validFrom, validTo, latitude, longitude, altitude);

            ShowMessage("Inserting wind data into wind database");
            dataBase.InsertWindField(w);
        }
    }

    return 1;
}

/** Writes an wind-field file in the NPPP-format
    @return 0 on success */
int CXMLWindFileReader::WriteWindFile(const novac::CString& fileName, const Meteorology::CWindDataBase& dataBase) {
    return dataBase.WriteToFile(fileName);
}


