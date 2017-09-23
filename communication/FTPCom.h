#pragma once

#include <afxinet.h>
#include <afxtempl.h>

#include "../FileInfo.h"
#include <PPPLib/CString.h>

namespace Communication
{
	class CFTPCom
	{
	public:
		CFTPCom(void);
		~CFTPCom(void);

		/**Connect to one FTP server
		*@param siteName - address of the FTP server
		*@param userName - user name for login
		*@param password - password for login
		*@param mode		 - Specifies passive(TRUE) or active mode(FASLE) for this FTP session. 
		*						 If set to TRUE, it sets the Win32 API dwFlag to INTERNET_FLAG_PASSIVE. 
		*						 Passive mode is for client behind a firewall; it is safer comparing 
		*						 with active mode.
		*/
		int Connect(LPCTSTR siteName, LPCTSTR userName, LPCTSTR password, BOOL mode= FALSE);

		int Disconnect();

		int UploadFile(LPCTSTR localFile, LPCTSTR remoteFile);

		/** Tries to download a file from the remote FTP-server to the local computer
			@return TRUE if successful
		*/
		BOOL DownloadAFile(LPCTSTR remoteFile, LPCTSTR fileFullName);

		int UpdateFile(LPCTSTR localFile, LPCTSTR remoteFile);

		/** Creates a direcotry in the current working-directory
			on the FTP-server. 
			@return non-zero if successful */
		int CreateDirectory(LPCTSTR remoteDirectory);

		/**find file in the current ftp folder
			@return TRUE if exists */
		int FindFile(novac::CString& fileName);

		/** Retrieves the list of files in the current directory
			@return 0 on success */
		int GetFileList(CList <novac::CString, novac::CString &> &fileNames);

		/** Retrieves the list of files in the current directory
			@return 0 on success */
		int GetFileList(const novac::CString &directory, CList <CFileInfo, CFileInfo &> &fileInfos);

		/**Set current directory
		*@param curDirName current directory name
		*@return TRUE  - success
		*/
		BOOL SetCurDirectory(LPCTSTR curDirName);

		/** Remove a folder*/
		BOOL DeleteFolder(const novac::CString& folder);

		/** Enter a folder*/
		BOOL EnterFolder(const novac::CString& folder);

		/** Go to top directory "/" */
		BOOL GotoTopDirectory();

		/** read response from the ftp server*/
		void ReadResponse(CInternetFile* file);

	private:
	
		CInternetSession* m_InternetSession;
		CFtpConnection* m_FtpConnection;
		novac::CString m_FTPSite;
		
		novac::CString m_ErrorMsg;
	};
}