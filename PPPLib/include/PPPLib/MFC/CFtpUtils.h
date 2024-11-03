#ifndef NOVAC_PPPLIB_CFTP_UTILS_H
#define NOVAC_PPPLIB_CFTP_UTILS_H

#include <fstream>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/VolcanoInfo.h>

namespace novac
{
struct CFileInfo
{
    // The file name (or name of the directory)
    std::string fileName = "";

    // True if this is a directory
    bool isDirectory = false;

    // The size of the file. Zero if this is a directory
    size_t fileSize = 0U;

    // The path to this file or directory (the parent items)
    std::string path = "";
};

class CFtpUtils
{
public:
    CFtpUtils();

    CFtpUtils(const CVolcanoInfo& volcanoes, unsigned int currentVolcano)
        : m_volcanoes(volcanoes), m_currentVolcano(currentVolcano)
    {
    }

    /** Splits a full directory path on the server (e.g. ftp://127.0.0.1/some/directory/path/ into the server
            component (ftp://127.0.0.1/) and the directory component (some/directory/path/) */
    void SplitPathIntoServerAndDirectory(const novac::CString& fullServerPath, std::string& server, std::string& directory);

    /** Reads file list info item from the provided string, parses the data and fills in the provided CFileInfo item
        @return true if the parse was successful. */
    bool ReadFtpDirectoryListing(const std::string& item, CFileInfo& result) const;

    /**  @return true if the provided string is a typical unix-file permission listing */
    bool IsFilePermissions(const std::string& item) const;

private:
    const CVolcanoInfo m_volcanoes;
    const unsigned int m_currentVolcano;
};
}  // namespace novac

#endif  // NOVAC_PPPLIB_CFTP_UTILS_H