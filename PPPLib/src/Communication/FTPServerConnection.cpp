#include <PPPLib/Communication/FTPServerConnection.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/Logging.h>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

// This is the global list of volcanoes
#include <PPPLib/VolcanoInfo.h>

#include <PPPLib/MFC/CCriticalSection.h>
#include <PPPLib/MFC/CSingleLock.h>
#include <PPPLib/MFC/CFtpUtils.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <PPPLib/ThreadUtils.h>
#include <Poco/Net/FTPClientSession.h>
#include <Poco/Net/NetException.h>
#include <Poco/Path.h>
#include <Poco/StreamCopier.h>

#include <fstream>
#include <chrono>
#include <iostream>


extern Configuration::CUserConfiguration g_userSettings;// <-- The settings of the user
extern novac::CVolcanoInfo g_volcanoes; // <-- A list of all known volcanoes

using namespace Communication;
using namespace novac;


struct ftpLogin
{
    std::string server;
    std::string userName;
    std::string password;
    std::string workingDirectory;
    int port = 21;
};

void LoginAndDownloadDataFromDir(ftpLogin login, std::string directory, novac::GuardedList<std::string>& downloadedFiles);

/** Downloads .pak files from the provided list of files on the already opened connection.
    Items are consumed from the downloadedQueue one at a time and appended to the 'downladedFiles' when ready.
    The files which have been downloaded are appended to the list 'downloadedFiles' */
void DownloadData(Poco::Net::FTPClientSession& ftp, novac::GuardedList<novac::CFileInfo>& downloadedQueue, novac::GuardedList<std::string>& downloadedFiles);

/** Downloads a specific file from the ftp session. The file is appended to the list of downloaded files upon success. */
void DownloadFile(Poco::Net::FTPClientSession& ftp, const novac::CFileInfo& fileInfo, novac::GuardedList<std::string>& downloadedFiles);

/** Downloads .pak files from the provided directory on the already opened connection.
    The files which have been downloaded are appended to the list 'downloadedFiles' */
void DownloadDataFromDir(Poco::Net::FTPClientSession& ftp, std::string directory, novac::GuardedList<std::string>& downloadedFiles);

CFTPServerConnection::CFTPServerConnection(void)
{
}

CFTPServerConnection::~CFTPServerConnection(void)
{
}

volatile double nMbytesDownloaded = 0.0;
double nSecondsPassed = 0.0;


static novac::GuardedValue nFTPThreadsRunning;


int CFTPServerConnection::DownloadDataFromFTP(const novac::CString& serverDir, const novac::CString& username,
    const novac::CString& password, std::vector<std::string>& pakFileList) {

    unsigned int nRounds = 0;

    ftpLogin login;
    login.userName = username.std_str();
    login.password = password.std_str();

    nFTPThreadsRunning.Zero();

    // Extract the name of the server and each of the sub-directories specified
    std::string directory;
    novac::CFtpUtils ftpUtil{ g_volcanoes, g_userSettings.m_volcano };
    ftpUtil.SplitPathIntoServerAndDirectory(serverDir, login.server, directory);


    // Make sure thath the temporary directory exists
    if (Filesystem::CreateDirectoryStructure(g_userSettings.m_tempDirectory)) {
        novac::CString userMessage;
        userMessage.Format("Could not create temp directory: %s", (const char*)g_userSettings.m_tempDirectory);
        ShowMessage(userMessage);
        return 1;
    }

    // This is (a thread safe) list of files which have been downloaded so far.
    novac::GuardedList<std::string> downloadedFiles;

    // download the data in this directory
    nFTPThreadsRunning.IncrementValue();
    std::thread downloadThread{ LoginAndDownloadDataFromDir, login, directory, std::ref(downloadedFiles) };

    // wait for all threads to terminate
    std::this_thread::sleep_for(std::chrono::milliseconds{ 500 });

    clock_t cStart = clock();
    while (nFTPThreadsRunning.GetValue() > 0) {
        if (++nRounds % 10 == 0) {
            nSecondsPassed = (double)(clock() - cStart) / CLOCKS_PER_SEC;

            novac::CString userMessage;
            if (nSecondsPassed < 100.0 * nMbytesDownloaded) {
                userMessage.Format("  %.0lf MBytes downloaded (<=> %.2lf MBytes/second)", nMbytesDownloaded, nMbytesDownloaded / nSecondsPassed);
            }
            else {
                userMessage.Format("  %.0lf MBytes downloaded (<=> %.2lf kBytes/second)", nMbytesDownloaded, 1024 * nMbytesDownloaded / nSecondsPassed);
            }
            ShowMessage(userMessage);
        }

        std::this_thread::sleep_for(std::chrono::milliseconds{ 500 });
    }
    downloadThread.join();

    // copy the data to the output list
    downloadedFiles.CopyTo<std::string>(pakFileList);

    return 0;
}


bool DownloadAFile(Poco::Net::FTPClientSession& ftp, const std::string& fullRemoteFileName, const std::string& localFileName) {
    try
    {
        std::istream& srcStream = ftp.beginDownload(fullRemoteFileName);
        std::ofstream dstStream{ localFileName, std::ios::binary | std::ios::out };

        Poco::StreamCopier::copyStream(srcStream, dstStream, 524288);
    }
    catch (std::exception& e)
    {
        ftp.endDownload();
        ShowMessage(e.what());
        return false;
    }
    ftp.endDownload();
    return true;
}

bool UploadAFile(Poco::Net::FTPClientSession& ftp, const std::string& localFileName, const std::string& fullRemoteFileName) {
    std::ifstream srcStream{ localFileName, std::ios::binary | std::ios::in };
    std::ostream& dstStream = ftp.beginUpload(fullRemoteFileName);

    while (!srcStream.eof())
    {
        dstStream.put((char)srcStream.get());
    }

    ftp.endDownload();
    return true;
}

bool DownloadFileList(Poco::Net::FTPClientSession& ftp, const std::string& directory, std::vector<novac::CFileInfo>& filesFound)
{
    novac::CFtpUtils ftpHelper{ g_volcanoes, g_userSettings.m_volcano };

    try
    {
        ftp.setTimeout(Poco::Timespan(10, 0));
        std::istream& fileListStream = ftp.beginList(directory, true);
        for (std::string line; std::getline(fileListStream, line); )
        {
            novac::CFileInfo result;
            if (ftpHelper.ReadFtpDirectoryListing(line, result))
            {
                result.path = directory;
                filesFound.push_back(result);
            }
        }
    }
    catch (std::exception& e)
    {
        std::cout << "Failed to download file list: " << e.what() << std::endl;
        ftp.endList();
        return false;
    }
    ftp.endList();

    return true;
}

bool DownloadFileList(Poco::Net::FTPClientSession& ftp, const std::string& directory, std::vector<std::string>& filesFound)
{
    novac::CFtpUtils ftpHelper{ g_volcanoes, g_userSettings.m_volcano };

    try
    {
        std::istream& fileListStream = ftp.beginList(directory, true);
        for (std::string line; std::getline(fileListStream, line); )
        {
            novac::CFileInfo result;
            if (ftpHelper.ReadFtpDirectoryListing(line, result))
            {
                filesFound.push_back(result.fileName);
            }
        }
    }
    catch (std::exception& e)
    {
        std::cout << "Failed to download file list: " << e.what() << std::endl;
        ftp.endList();
        return false;
    }
    ftp.endList();
    return true;
}

//  Creates a new connection to the provided server and downloads a file list will all files and sub-directories to the provided directory.
std::vector<novac::CFileInfo> ListContentsOfDir(ftpLogin login, std::string directory) {
    std::vector<novac::CFileInfo> filesFound;

    try
    {
        Poco::Net::FTPClientSession ftp;
        ftp.open(login.server, (Poco::UInt16)login.port, login.userName, login.password);

        if (!DownloadFileList(ftp, directory, filesFound))
        {
            ShowMessage("Failed to retrieve file list in directory: " + directory);
        }
        ftp.close();
    }
    catch (Poco::Net::FTPException& ex)
    {
        ShowMessage("Failed to connect: " + ex.message());
    }
    catch (std::exception& e)
    {
        ShowMessage(e.what());
    }

    return filesFound;
}

/** Downloads all the data-files from the provided directory on the given FTP server */
void LoginAndDownloadDataFromDir(ftpLogin login, std::string directory, novac::GuardedList<std::string>& downloadedFiles) {

    try
    {
        novac::CString userMessage;
        userMessage.Format("Started thread #%d", nFTPThreadsRunning.GetValue());
        ShowMessage(userMessage);

        // Search for files in the specified directory
        std::vector<novac::CFileInfo> filesFound = ListContentsOfDir(login, directory);

        if (filesFound.size() == 0) {
            nFTPThreadsRunning.DecrementValue();
            return; // nothing more to do..
        }

        // Generate the download queue
        novac::GuardedList<novac::CFileInfo> downloadQueue;
        for (novac::CFileInfo file : filesFound) {
            downloadQueue.AddItem(file);
        }

        // Fork into a number of threads and start downloading the files
        std::vector<std::shared_ptr<std::thread>> downloadThreads;
        std::vector<std::unique_ptr<Poco::Net::FTPClientSession>> connections{ g_userSettings.m_maxThreadNum };

        for (int threadIdx = 0; threadIdx < g_userSettings.m_maxThreadNum; ++threadIdx) {

            // Create the connection
            try {
                connections[threadIdx].reset(new Poco::Net::FTPClientSession());
                connections[threadIdx]->open(login.server, (Poco::UInt16)login.port, login.userName, login.password);
                connections[threadIdx]->setTimeout(Poco::Timespan(60, 0)); // 60 seconds timeout

                auto t = std::make_shared<std::thread>(DownloadData, std::ref(*connections[threadIdx]), std::ref(downloadQueue), std::ref(downloadedFiles));

                downloadThreads.push_back(t);
            }
            catch (Poco::Net::FTPException& ex) {
                ShowMessage("Failed to connect to ftp server. Message was: " + ex.message());
            }
        }

        // Wait for all the threads to terminate
        for (int threadIdx = 0; threadIdx < g_userSettings.m_maxThreadNum; ++threadIdx) {
            downloadThreads[threadIdx]->join();
        }

        // Close the connections
        for (int threadIdx = 0; threadIdx < g_userSettings.m_maxThreadNum; ++threadIdx) {
            connections[threadIdx]->close();
        }

        nFTPThreadsRunning.DecrementValue();
    }
    catch (std::exception& ex)
    {
        ShowMessage(ex.what());
        nFTPThreadsRunning.DecrementValue();
    }
}

bool ParseDate(const std::string& str, CDateTime& result)
{
    int year = 0;
    int month = 0;
    int day = 0;

    int nNumbers = sscanf(str.c_str(), "%d.%d.%d", &year, &month, &day);
    if (nNumbers == 3) {
        result = CDateTime(year, month, day, 0, 0, 0);
        return true;
    }
    else {
        return false;
    }
}

bool IsPakFile(const novac::CFileInfo& item) {
    return novac::Equals(novac::Right(item.fileName, 4), ".pak");
}

void DownloadData(Poco::Net::FTPClientSession& ftp, novac::GuardedList<novac::CFileInfo>& downloadedQueue, novac::GuardedList<std::string>& downloadedFiles) {
    novac::CFileInfo nextDownloadItem;
    while (downloadedQueue.PopFront(nextDownloadItem)) {
        if (nextDownloadItem.isDirectory) {
            CDateTime date;
            if (ParseDate(nextDownloadItem.fileName, date)) {
                if (date <= g_userSettings.m_toDate && g_userSettings.m_fromDate <= date) {
                    DownloadDataFromDir(ftp, nextDownloadItem.path + nextDownloadItem.fileName, downloadedFiles);
                }
                else {
                    // ShowMessage("Ignoring directory '" + nextDownloadItem.path + nextDownloadItem.fileName + "'. The date is out of scope");
                }
            }
            else {
                DownloadDataFromDir(ftp, nextDownloadItem.path + nextDownloadItem.fileName, downloadedFiles);
            }
        }
        else if (IsPakFile(nextDownloadItem)) {
            DownloadFile(ftp, nextDownloadItem, downloadedFiles);
        }
    }
}

void DownloadFile(Poco::Net::FTPClientSession& ftp, const novac::CFileInfo& fileInfo, novac::GuardedList<std::string>& downloadedFiles) {
    novac::CString localFileName, serial;
    novac::CString userMessage;
    CDateTime start;
    int channel;
    MEASUREMENT_MODE mode;

    // if this is a .pak-file then check the date when it was created
    if (novac::CFileUtils::GetInfoFromFileName(fileInfo.fileName, start, serial, channel, mode)) {
        if (start <= g_userSettings.m_toDate && g_userSettings.m_fromDate <= start) {
            // the creation date is between the start and the stop dates. Download the file
            localFileName.Format("%s%c%s", (const char*)g_userSettings.m_tempDirectory, Poco::Path::separator(), fileInfo.fileName.c_str());
            if (Filesystem::IsExistingFile(localFileName)) {
                userMessage.Format("File %s is already downloaded", (const char*)localFileName);
                ShowMessage(userMessage);
                downloadedFiles.AddItem(localFileName.std_str());
            }
            else {
                // download the file
                if (DownloadAFile(ftp, fileInfo.path + "/" + fileInfo.fileName, localFileName.std_str())) {
                    nMbytesDownloaded += fileInfo.fileSize / 1048576.0;
                    downloadedFiles.AddItem(localFileName.std_str());
                }
            }
        }
    }
}

void DownloadDataFromDir(Poco::Net::FTPClientSession& ftp, std::string directory, novac::GuardedList<std::string>& downloadedFiles) {

    std::vector<novac::CFileInfo> filesFound;

    // Search for files in the specified directory
    if (!DownloadFileList(ftp, directory, filesFound))
    {
        return;
    }

    // download each of the files .pak-files found or enter the all the sub-directories
    for (const novac::CFileInfo& fileInfo : filesFound) {
        if (IsPakFile(fileInfo)) {
            DownloadFile(ftp, fileInfo, downloadedFiles);
        }
        else if (g_userSettings.m_includeSubDirectories_FTP && fileInfo.isDirectory) {
            CDateTime date;
            if (ParseDate(fileInfo.fileName, date)) {
                if (date < g_userSettings.m_fromDate) {
                    continue;
                }
            }

            // start downloading using the same thread as we're running in
            novac::CString subDir = directory + fileInfo.fileName + "/";
            DownloadDataFromDir(ftp, subDir.std_str(), downloadedFiles);
        }
    }

    return;
}

/** Retrieves the list of files in a given directory on the FTP-server
    @return 0 on successful connection and completion of the download */
int CFTPServerConnection::DownloadFileListFromFTP(const novac::CString& serverDir, std::vector<std::string>& fileList, const novac::CString& username, const novac::CString& password) {
    // Extract the name of the server and each of the sub-directories specified
    std::string server, directory;
    novac::CFtpUtils ftpUtil{ g_volcanoes, g_userSettings.m_volcano };
    ftpUtil.SplitPathIntoServerAndDirectory(serverDir, server, directory);

    // create a new connection
    Poco::Net::FTPClientSession ftp;

    // connect to the server
    try
    {
        ftp.open(server, 21, username.std_str(), password.std_str());
    }
    catch (Poco::Net::FTPException& ex)
    {
        ShowMessage("Failed to connect: " + ex.message());
        return 1;
    }

    DownloadFileList(ftp, serverDir.std_str(), fileList);

    return 0;
}

/** Downloads a single file from the given FTP-server
    @return 0 on successful connection and completion of the download
*/
int CFTPServerConnection::DownloadFileFromFTP(const novac::CString& remoteFileName, const novac::CString& localFileName,
    const novac::CString& username, const novac::CString& password) {
    CDateTime now;

    // Extract the name of the server and the login
    ftpLogin login;
    login.userName = username.std_str();
    login.password = password.std_str();

    // Extract the name of the server and each of the sub-directories specified
    novac::CString subString;
    novac::CString directory;
    int indexOfSlash[128];
    int nSlashesFound = 0; // the number of slashes found in the 'serverDir' - path
    if (remoteFileName.Find("ftp://") != -1) {
        indexOfSlash[0] = 5;
    }
    else {
        indexOfSlash[0] = 0;
    }
    while (-1 != (indexOfSlash[nSlashesFound + 1] = remoteFileName.Find('/', indexOfSlash[nSlashesFound] + 1))) {
        ++nSlashesFound;
    }
    subString.Format(remoteFileName.Left(indexOfSlash[1]));
    directory.Format(remoteFileName.Right(remoteFileName.GetLength() - indexOfSlash[1] - 1));
    login.server = subString.Right(subString.GetLength() - indexOfSlash[0] - 1).std_str();

    // create a new connection
    Poco::Net::FTPClientSession ftp;

    // connect to the server
    try
    {
        ftp.open(login.server, (Poco::UInt16)login.port, login.userName, login.password);
    }
    catch (Poco::Net::FTPException& ex)
    {
        ShowMessage("Failed to connect: " + ex.message());
        return 1; // failed to connect!
    }

    // Download the file
    if (!DownloadAFile(ftp, remoteFileName.std_str(), localFileName.std_str())) {
        novac::CString errorMessage;
        errorMessage.Format("Failed to download remote file %s from FTP server", (const char*)remoteFileName);
        ShowMessage(errorMessage);
    }

    // disconnect
    ftp.close();

    return 0;
}

/** Uploads result-files to the given FTP-server
    @return 0 on success otherwise non-zero
*/
int CFTPServerConnection::UploadResults(const novac::CString& server, const novac::CString& username,
    const novac::CString& password, novac::CList <novac::CString, novac::CString&>& fileList) {
    CDateTime now;

    // Extract the name of the server and the login
    ftpLogin login;
    login.userName = username.std_str();
    login.password = password.std_str();
    login.server = server.std_str();

    novac::CString volcanoName, directoryName, remoteFile, errorMessage;

    // create a new connection
    Poco::Net::FTPClientSession ftp;

    // connect to the server
    try
    {
        ftp.open(login.server, (Poco::UInt16)login.port, login.userName, login.password);
    }
    catch (Poco::Net::FTPException& ex)
    {
        ShowMessage("Failed to connect: " + ex.message());
        return 1;
    }

    // Enter the volcano's directory
    volcanoName.Format("%s", (const char*)g_volcanoes.GetSimpleVolcanoName(g_userSettings.m_volcano));
    ftp.setWorkingDirectory(volcanoName.std_str());

    // Enter the upload-directory
    now.SetToNow();
    directoryName.Format("PostProcessed_BETA_%04d%02d%02d", now.year, now.month, now.day);
    try {
        ftp.createDirectory(directoryName.std_str());
    }
    catch (Poco::Net::FTPException& ex)
    {
        ShowMessage("Failed to create directory '" + directoryName + "' on ftp-server. Error was: " + ex.message());
        return 2;
    }
    ftp.setWorkingDirectory(directoryName.std_str());

    // Upload the files
    auto p = fileList.GetHeadPosition();
    while (p != nullptr) {
        // Get the local name and path of the file to upload
        novac::CString& localFile = fileList.GetNext(p);

        // Get the file-name to upload the file to...
        remoteFile.Format(localFile);
        CFileUtils::GetFileName(remoteFile);

        if (!UploadAFile(ftp, localFile.std_str(), remoteFile.std_str())) {
            errorMessage.Format("Failed to upload local file %s to FTP server", (const char*)localFile);
            ShowMessage(errorMessage);
        }
    }

    // disconnect
    ftp.close();

    return 0;
}
