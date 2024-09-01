#ifndef FTPSERVER_DATADOWNLOAD_H
#define FTPSERVER_DATADOWNLOAD_H

#include <vector>
#include <PPPLib/MFC/CList.h>
#include <PPPLib/MFC/CString.h>

namespace Communication
{

class CFTPServerConnection
{
public:

    // -----------------------------------------------------------
    // ---------------------- PUBLIC DATA ------------------------
    // -----------------------------------------------------------


    // -----------------------------------------------------------
    // --------------------- PUBLIC METHODS ----------------------
    // -----------------------------------------------------------

    /** Downloads .pak - files from the given FTP-server
        @return 0 on successful connection and completion of the list
    */
    int DownloadDataFromFTP(const novac::CString& server, const novac::CString& username,
        const novac::CString& password, std::vector<std::string>& pakFileList);

    /** Downloads a single file from the given FTP-server
        @return 0 on successful connection and completion of the download
    */
    int DownloadFileFromFTP(const novac::CString& remoteFileName, const novac::CString& localFileName,
        const novac::CString& username, const novac::CString& password);

    /** Retrieves the list of files (but no directories) in a given directory on the FTP-server
        @return 0 on successful connection and completion of the download
    */
    int DownloadFileListFromFTP(const novac::CString& serverDir, std::vector<std::string>& fileList,
        const novac::CString& username, const novac::CString& password);

    /** Uploads result-files to the given FTP-server
        @return 0 on success otherwise non-zero
    */
    int UploadResults(const novac::CString& server, const novac::CString& username,
        const novac::CString& password, novac::CList <novac::CString, novac::CString&>& fileList);

private:
    // -----------------------------------------------------------
    // ---------------------- PRIVATE DATA -----------------------
    // -----------------------------------------------------------


    // -----------------------------------------------------------
    // --------------------- PRIVATE METHODS ---------------------
    // -----------------------------------------------------------

};
}

#endif