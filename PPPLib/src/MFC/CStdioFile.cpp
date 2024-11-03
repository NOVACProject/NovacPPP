#include "PPPLib/MFC/CStdioFile.h"
#include <vector>

namespace novac
{
    CStdioFile::CStdioFile()
    {
    }

    CStdioFile::~CStdioFile()
    {
        Close();
    }

    bool CStdioFile::Open(const char* fileName, unsigned nOpenFlags)
    {
        return Open(fileName, nOpenFlags, nullptr);
    }

    bool CStdioFile::Open(const char* fileName, unsigned nOpenFlags, CFileException* /*ex*/)
    {
        m_f.open(fileName, (std::ios_base::openmode)nOpenFlags);

        return m_f.is_open();
    }

    void CStdioFile::Close()
    {
        m_f.close();
    }

    bool CStdioFile::ReadString(char* destination, unsigned int nMax)
    {
        if (!m_f.is_open() || m_f.eof())
        {
            return false;
        }
        else
        {
            m_f.getline(destination, nMax, '\n');
            return true;
        }
    }

    bool CStdioFile::ReadString(CString& destination)
    {
        if (!m_f.is_open() || m_f.eof())
        {
            return false;
        }
        else
        {
            // TODO: This is formally a bit of cheating, this doesn't read the full string....
            const size_t maxSize = 65535;

            std::vector<char> tmpBuffer(maxSize);
            m_f.getline(tmpBuffer.data(), maxSize);

            destination.SetData(tmpBuffer.data());

            return true;
        }
    }
}