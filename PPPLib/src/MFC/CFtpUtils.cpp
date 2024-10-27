#include "PPPLib/MFC/CFtpUtils.h"
#include <sstream>
#include <algorithm>

namespace novac
{
CFtpUtils::CFtpUtils()
    : m_volcanoes(CVolcanoInfo()), m_currentVolcano(0)
{}

void CFtpUtils::SplitPathIntoServerAndDirectory(const novac::CString& fullServerPath, std::string& server, std::string& directory)
{
    novac::CString subString;

    int indexOfSlash[128];
    int nSlashesFound = 0; // the number of slashes found in the 'serverDir' - path

    if (fullServerPath.Find("ftp://") != -1)
    {
        indexOfSlash[0] = 5;
    }
    else
    {
        indexOfSlash[0] = 0;
    }

    while (-1 != (indexOfSlash[nSlashesFound + 1] = fullServerPath.Find('/', indexOfSlash[nSlashesFound] + 1)))
    {
        ++nSlashesFound;
    }
    subString.Format(fullServerPath.Left(indexOfSlash[1]));

    server = (subString.Right(subString.GetLength() - indexOfSlash[0] - 1)).std_str();

    if (nSlashesFound == 1)
    {
        directory = m_volcanoes.GetSimpleVolcanoName(m_currentVolcano).std_str() + "/";
    }
    else
    {
        directory = (fullServerPath.Right(fullServerPath.GetLength() - indexOfSlash[1] - 1)).std_str();
    }
}

void split(const std::string& s, char delim, std::vector<std::string>& result)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        if (item.size() == 0 || (item.size() == 1 && item.at(0) == delim))
        {
            continue;
        }
        else
        {
            result.push_back(item);
        }
    }
}

bool CFtpUtils::ReadFtpDirectoryListing(const std::string& item, CFileInfo& result) const
{
    std::vector<std::string> parts;
    split(item, ' ', parts);

    if (parts.size() < 3 || !IsFilePermissions(parts[0]))
    {
        return false;
    }

    std::string fileName = parts[parts.size() - 1];

    // remove trailing newline characters
    fileName.erase(std::remove(fileName.begin(), fileName.end(), '\r'), fileName.end());

    result.fileName = fileName;

    result.isDirectory = (parts[0].at(0) == 'd');

    result.fileSize = std::atol(parts[4].c_str());

    return true;
}

bool CFtpUtils::IsFilePermissions(const std::string& item) const
{
    if (item.size() != 10)
    {
        return false;
    }
    if (item.find_first_not_of("drwx-") != std::string::npos)
    {
        return false;
    }
    if (item[0] != 'd' && item[0] != '-')
    {
        return false;
    }
    if (item[1] != 'r' && item[1] != '-')
    {
        return false;
    }
    if (item[2] != 'w' && item[2] != '-')
    {
        return false;
    }
    if (item[3] != 'x' && item[3] != '-')
    {
        return false;
    }

    return true;
}
}