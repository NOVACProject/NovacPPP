#pragma once

class CFileInfo
{
public:
	CFileInfo(void);
	~CFileInfo(void);
	
	/** constructs a new file-info with a given file name 
		and file size */
	CFileInfo(CString fileName, long fileSize, bool isDirectory = false);

	/** the name of the file */
	CString  m_fileName;
	
	/** The full path and file-name of the file */
	CString	 m_fullFileName;

	/** The size of the file, in bytes */
	long	m_fileSize;
	
	/** True if this is a directory */
	bool	m_isDirectory;
};
