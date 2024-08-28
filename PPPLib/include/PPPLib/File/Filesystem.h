#ifndef NOVACPPP_FILESYSTEM_FILESYSTEM_H
#define NOVACPPP_FILESYSTEM_FILESYSTEM_H

#include <string>
#include <vector>
#include <PPPLib/MFC/CString.h>
#include <SpectralEvaluation/DateTime.h>

namespace Filesystem
{
    struct FileSearchCriterion
    {
        FileSearchCriterion()
            : fileExtension("") {
        }

        novac::CDateTime startTime;
        novac::CDateTime endTime;
        std::string fileExtension;
    };

    /** Scans through the given directory in search for files with the given criteria.
        @param path - the directory (on the local computer) where to search for files.
        @param includeSubdirectories If set to true then sub-directories of the provided path will also be searched.
        @param fileList Will be appended with the path's and file-names of the found .pak-files.
        @param criteria If not null, then this is used to filter the list of files. */
    void SearchDirectoryForFiles(const novac::CString& path, bool includeSubdirectories, std::vector<std::string>& fileList, FileSearchCriterion* criteria = nullptr);

    /** A simple function to find out whether a given file exists or not.
        @param - The filename (including path) to the file.
        @return 0 if the file does not exist.
        @return 1 if the file exist. */
    bool IsExistingFile(const novac::CString& fileName);

    /** Creates a directory structure according to the given path.
            @return 0 on success. */
    int CreateDirectoryStructure(const novac::CString& path);

    /** AppendPathSeparator returns a string which does end with the path-separator character of the current system. */
    novac::CString AppendPathSeparator(novac::CString path);
}

#endif // !NOVACPPP_FILESYSTEM_FILESYSTEM_H
