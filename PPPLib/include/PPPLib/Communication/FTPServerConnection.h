#pragma once

#include <vector>
#include <PPPLib/MFC/CList.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/Configuration/UserConfiguration.h>

#include <SpectralEvaluation/Log.h>

namespace Communication
{

class CFTPServerConnection
{
public:

    CFTPServerConnection(novac::ILogger& log, const Configuration::CUserConfiguration& userSettings)
        : m_log(log), m_userSettings(userSettings)
    {
    }

    // -----------------------------------------------------------
    // ---------------------- PUBLIC DATA ------------------------
    // -----------------------------------------------------------


    // -----------------------------------------------------------
    // --------------------- PUBLIC METHODS ----------------------
    // -----------------------------------------------------------

    /** Downloads .pak - files from the given FTP-server
        @return 0 on successful connection and completion of the list
    */
    int DownloadDataFromFTP(
        novac::LogContext context,
        const std::string& server,
        const std::string& username,
        const std::string& password,
        std::vector<std::string>& pakFileList);

    /** Downloads a single file from the given FTP-server
        @return 0 on successful connection and completion of the download
    */
    int DownloadFileFromFTP(
        novac::LogContext context,
        const novac::CString& remoteFileName,
        const novac::CString& localFileName,
        const novac::CString& username,
        const novac::CString& password);

    /** Retrieves the list of files (but no directories) in a given directory on the FTP-server
        @return 0 on successful connection and completion of the download
    */
    int DownloadFileListFromFTP(const novac::CString& serverDir, std::vector<std::string>& fileList,
        const novac::CString& username, const novac::CString& password);

    /** Uploads result-files to the given FTP-server
        @return 0 on success otherwise non-zero
    */
    int UploadResults(novac::LogContext context, const novac::CString& server, const novac::CString& username,
        const novac::CString& password, const std::vector<std::string>& fileList);

private:
    // -----------------------------------------------------------
    // ---------------------- PRIVATE DATA -----------------------
    // -----------------------------------------------------------

    novac::ILogger& m_log;

    const Configuration::CUserConfiguration& m_userSettings;

    // -----------------------------------------------------------
    // --------------------- PRIVATE METHODS ---------------------
    // -----------------------------------------------------------

};
}
