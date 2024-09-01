#ifndef NOVAC_PPPLIB_CSTDIOFILE_H
#define NOVAC_PPPLIB_CSTDIOFILE_H

#include <fstream>
#include <PPPLib/MFC/CString.h>

namespace novac
{
struct CFileException
{
};

class CStdioFile
{
public:
    CStdioFile();
    ~CStdioFile();

    CStdioFile(const CStdioFile&) = delete;
    CStdioFile& operator=(const CStdioFile&) = delete;

    enum ModeFlag
    {
        modeRead = std::ios_base::in,
        modeWrite = std::ios_base::out,
        modeReadWrite = std::ios_base::in | std::ios_base::out,
        typeBinary = std::ios_base::binary,
        typeText = 0,
        modeCreate = 0
    };

    bool Open(const char* fileName, unsigned nOpenFlags);

    bool Open(const char* fileName, unsigned nOpenFlags, CFileException* ex);

    void Close();

    // Reads text data into a buffer, up to a limit of nMax–1 characters, from the file associated with the CStdioFile object.
    // This reads until a new-line is encountered OR nMax-1 characters have been read.
    // Return FALSE if end-of-file was reached without reading any data.
    bool ReadString(char* destination, unsigned int nMax);

    // Reads text data into CString receiver, from the file associated with the CStdioFile object.
    // This reads until a new-line is encountered.
    bool ReadString(CString& destination);

private:
    std::ifstream m_f;


};
}

#endif  // NOVAC_PPPLIB_CSTDIOFILE_H