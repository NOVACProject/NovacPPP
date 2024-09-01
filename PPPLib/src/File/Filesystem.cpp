#include <PPPLib/File/Filesystem.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <Poco/DirectoryIterator.h>
#include <Poco/Exception.h>

using namespace novac;

void ShowMessage(const char message[]);

namespace Filesystem
{

void SearchDirectoryForFiles(const novac::CString& path, bool includeSubdirectories, std::vector<std::string>& fileList, FileSearchCriterion* criteria)
{
    try
    {
        Poco::DirectoryIterator dir{ path.std_str() };
        Poco::DirectoryIterator end;

        while (dir != end) {
            novac::CString fileName, fullFileName;
            fileName.Format("%s", dir.name().c_str());
            fullFileName.Format("%s/%s", (const char*)path, dir.name().c_str());
            const bool isDirectory = dir->isDirectory();

            ++dir; // go to next file in the directory

            if (novac::Equals(dir.name(), ".") || novac::Equals(dir.name(), "..")) {
                continue;
            }

            // if this is a directory...
            if (isDirectory && includeSubdirectories)
            {
                SearchDirectoryForFiles(fullFileName, includeSubdirectories, fileList, criteria);
                continue;
            }

            // check that this file is in the time-interval that we should evaluate spectra.
            if (nullptr != criteria)
            {
                if (criteria->endTime > criteria->startTime)
                {
                    int channel;
                    CDateTime startTime;
                    MEASUREMENT_MODE mode;
                    novac::CString serial;
                    novac::CFileUtils::GetInfoFromFileName(fileName, startTime, serial, channel, mode);

                    if (startTime < criteria->startTime || criteria->endTime < startTime)
                    {
                        continue;
                    }
                }

                if (criteria->fileExtension.size() > 0)
                {
                    if (!novac::Equals(fileName.Right(criteria->fileExtension.size()), criteria->fileExtension))
                    {
                        continue;
                    }
                    if (novac::Equals(criteria->fileExtension, ".pak"))
                    {
                        if (novac::CFileUtils::IsIncompleteFile(fileName))
                        {
                            continue;
                        }
                    }
                }

                // We've passed all the tests for the .pak-file.
                // Append the found file to the list of files to split and evaluate...
                std::string filenameStr = fullFileName.ToStdString();
                fileList.push_back(filenameStr);
            }
        }
    }
    catch (Poco::PathNotFoundException& e)
    {
        ShowMessage(e.what());
    }
}


bool IsExistingFile(const novac::CString& fileName)
{
    try
    {
        Poco::File file(fileName.c_str());
        return file.exists();
    }
    catch (const std::exception& e)
    {
        novac::CString message;
        message.Format("Exception happened when searching for file: '%s', message: '%s'", fileName.c_str(), e.what());
        ShowMessage(message);
        return false;
    }
}

int CreateDirectoryStructure(const novac::CString& path)
{
    try
    {
        Poco::File directory(path.c_str());
        directory.createDirectories();

        if (directory.exists()) {
            return 0;
        }
        else {
            return 1; // error
        }
    }
    catch (std::exception& e)
    {
        novac::CString message = "Failed to create directory: ";
        message.AppendFormat("%s", e.what());
        ShowMessage(message);
        return 1;
    }
}

std::string GetAbsolutePathFromRelative(const std::string& path, const std::string& baseDirectory)
{
    Poco::Path p{ path };
    Poco::Path exePath{ baseDirectory };
    Poco::Path absolutePath = p.absolute(exePath);
    return absolutePath.toString();
}

novac::CString AppendPathSeparator(novac::CString path)
{
    if (path.Right(1) == "/" || path.Right(1) == "\\")
    {
        return path;
    }

    char separator[8];
#ifdef _MSC_VER
    sprintf_s(separator, "%c", Poco::Path::separator());
#else
    sprintf(separator, "%c", Poco::Path::separator());
#endif

    return path.Append(separator);
}

std::string AppendPathSeparator(std::string path)
{
    CString p{ path };
    p = AppendPathSeparator(p);
    return p.std_str();
}
}