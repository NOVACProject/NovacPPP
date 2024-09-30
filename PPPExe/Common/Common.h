#pragma once

#include <PPPLib/PPPLib.h>
#include <PPPLib/Logging.h>
#include <PPPLib/MFC/CString.h>


// ---------------------------------------------------------------
// --------------------------- LOGGING ---------------------------
// ---------------------------------------------------------------

class PocoLogger : public novac::ILogger
{
public:
    virtual void Debug(const std::string& message) override;

    virtual void Information(const std::string& message) override;

    virtual void Error(const std::string& message) override;
};

/** Update the top line of list box */
void UpdateMessage(const novac::CString& message);


// -------------------------------------------------------
// ---------------- CLASS COMMON.H -----------------------
// -------------------------------------------------------

/** The class <b>Common</b> contains misc. functions that we need but
    do not really fit anywhere else.
*/

class Common
{

public:
    Common();

    // --------------------------------------------------------------------
    // ------------------------- FILE -------------------------------------
    // --------------------------------------------------------------------

    /** Get file size in bytes.
        @param - The file name (including path)
        @return - The file size (in bytes)
    */
    static long RetrieveFileSize(novac::CString& fileName);

    /** Compares two files to see if their contents are the same */
    static bool AreIdenticalFiles(const novac::CString& fileName1, const novac::CString& fileName2);

    /** If there's a file with the given input name, then it will be renamed to
        PATH\\FILENAME_creationDate_creationTime.FILEENDING */
    static bool ArchiveFile(const novac::CString& fileName);

    // --------------------------------------------------------------------
    // ------------------------- PATH -------------------------------------
    // --------------------------------------------------------------------

    /** Take out the file name from a long path containing both a directory name and a file name.
        The path separator can be either '/' or '\'.
        @param fileName path of the file, will be set to only contain the filename, without the path . */
    static void GetFileName(novac::CString& fileName);

    /** Take out the directory from a long path name.
        @param fileName - the complete path of the file */
    static void GetDirectory(novac::CString& fileName);

    /** Copies the file to the new file location */
    static void CopyFile(const novac::CString& oldName, const novac::CString& newName);

    /** m_exePath will be set to the path where the application resides,
        This is set once at application startup and should never be changed. */
    const novac::CString m_exePath;

    /** m_exeFileName will be set to the filename of the program.
        This is set once at application startup and should never be changed. */
    const novac::CString m_exeFileName;

};// end of class Common
