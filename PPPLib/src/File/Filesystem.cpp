#include <PPPLib/File/Filesystem.h>
#include <SpectralEvaluation/File/File.h>
#include <PPPLib/MFC/CFileUtils.h>
#include <Poco/DirectoryIterator.h>
#include <Poco/Exception.h>

void ShowMessage(const char message[]);

namespace Filesystem
{

void SearchDirectoryForFiles(const std::string& path, bool includeSubdirectories, std::vector<std::string>& fileList, FileSearchCriterion* criteria)
{
    try
    {
        Poco::DirectoryIterator dir{ path };
        Poco::DirectoryIterator end;

        while (dir != end)
        {
            const std::string filename = dir.name();
            const std::string filenameIncludingPath = path + "/" + dir.name();
            const bool isDirectory = dir->isDirectory();

            ++dir; // go to next file in the directory

            if (novac::EqualsIgnoringCase(dir.name(), ".") || novac::EqualsIgnoringCase(dir.name(), ".."))
            {
                continue;
            }

            // if this is a directory...
            if (isDirectory && includeSubdirectories)
            {
                SearchDirectoryForFiles(filenameIncludingPath, includeSubdirectories, fileList, criteria);
                continue;
            }

            // check that this file is in the time-interval that we should evaluate spectra.
            if (nullptr != criteria)
            {
                if (criteria->endTime > criteria->startTime)
                {
                    int channel;
                    novac::CDateTime startTime;
                    novac::MeasurementMode mode;
                    novac::CString serial;
                    novac::CFileUtils::GetInfoFromFileName(filename, startTime, serial, channel, mode);

                    if (startTime < criteria->startTime || criteria->endTime < startTime)
                    {
                        continue;
                    }
                }

                if (criteria->fileExtension.size() > 0)
                {
                    if (filename.size() <= criteria->fileExtension.size())
                    {
                        continue;
                    }
                    const std::string currentFileExtension = novac::GetFileExtension(filename);
                    if (!novac::EqualsIgnoringCase(currentFileExtension, criteria->fileExtension))
                    {
                        continue;
                    }
                    if (novac::EqualsIgnoringCase(criteria->fileExtension, ".pak"))
                    {
                        if (novac::CFileUtils::IsIncompleteFile(filename))
                        {
                            continue;
                        }
                    }
                }

                // We've passed all the tests for the .pak-file.
                // Append the found file to the list of files to split and evaluate...
                fileList.push_back(filenameIncludingPath);
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

        if (directory.exists())
        {
            return 0;
        }
        else
        {
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
    const char pathSeparatorChar = '/'; // always use the forward slash for path separation (works on both windows and linux)

    char separator[8];
#ifdef _MSC_VER
    sprintf_s(separator, "%c", pathSeparatorChar);
#else
    sprintf(separator, "%c", pathSeparatorChar);
#endif

    return path.Append(separator);
}

std::string AppendPathSeparator(std::string path)
{
    novac::CString p{ path };
    p = AppendPathSeparator(p);
    return p.std_str();
}
}