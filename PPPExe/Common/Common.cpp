#include "Common.h"

#include <algorithm>
#include <iostream>
#include <time.h>

#include <Poco/File.h>
#include <Poco/DirectoryIterator.h>
#include <Poco/DateTime.h>
#include <Poco/Message.h>
#include <Poco/Logger.h>

// include the global settings
#include <PPPLib/VolcanoInfo.h>

#include <SpectralEvaluation/DateTime.h>

extern novac::CVolcanoInfo g_volcanoes; // <-- the list of volcanoes

extern std::string s_exePath;
extern std::string s_exeFileName;

#undef min
#undef max

void UpdateMessage(const novac::CString& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message.std_str());
}

void ShowMessage(const novac::CString& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message.std_str());
}
void ShowMessage(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message);
}
void ShowMessage(const novac::CString& message, novac::CString connectionID)
{
    novac::CString msg;

    msg.Format("<%s> : %s", (const char*)connectionID, (const char*)message);

    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(msg.std_str());
}

void ShowMessage(const char message[])
{
    novac::CString msg;
    msg.Format("%s", message);
    ShowMessage(msg);
}

void ShowError(const novac::CString& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.fatal(message.std_str());
}
void ShowError(const char message[])
{
    novac::CString msg;
    msg.Format("%s", message);
    ShowError(msg);
}

void PocoLogger::Debug(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.debug(message);
}

void PocoLogger::Information(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message);
}

void PocoLogger::Error(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.fatal(message);
}

Common::Common()
    :m_exePath(s_exePath), m_exeFileName(s_exeFileName)
{

}

void Common::GetFileName(novac::CString& fileName)
{
    // look for slashes in the path
    int position = std::max(fileName.ReverseFind('\\'), fileName.ReverseFind('/'));
    int length = static_cast<int>(fileName.GetLength());
    fileName = fileName.Right(length - position - 1);
}

/** Take out the directory from a long path name.
    @param fileName - the complete path of the file */
void Common::GetDirectory(novac::CString& fileName)
{
    int position = fileName.ReverseFind('\\');
    if (position >= 0)
    {
        fileName = fileName.Left(position + 1);
    }
}

void Common::CopyFile(const novac::CString& oldName, const novac::CString& newName)
{
    Poco::File oldFile(oldName.std_str());

    oldFile.copyTo(newName.std_str());
}

long Common::RetrieveFileSize(novac::CString& fileName)
{
    Poco::File file(fileName.std_str());
    return static_cast<long>(file.getSize());
}


/** Compares two files to see if their contents are the same */
bool Common::AreIdenticalFiles(const novac::CString& fileName1, const novac::CString& fileName2)
{
    if (Equals(fileName1, fileName2))
        return true; // a file is always identical to itself

    FILE* f1 = fopen(fileName1, "r");
    if (f1 == NULL)
        return false;

    FILE* f2 = fopen(fileName2, "r");
    if (f2 == NULL)
        return false;

    while (1)
    {
        int c1 = fgetc(f1);
        int c2 = fgetc(f2);

        if (c1 == EOF)
        {
            fclose(f1); fclose(f2);
            if (c2 == EOF)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        if (c1 != c2)
            return false;
    }


    fclose(f1); fclose(f2);

    // should never reach this point...
    return false;
}

/** If there's a file with the given input name, then it will be renamed to
    PATH\\FILENAME_creationDate_creationTime.FILEENDING */
bool Common::ArchiveFile(const novac::CString& fileName)
{
    novac::CString newFileName, errorMsg;

    // Search for the file
    Poco::File oldFile(fileName.std_str());

    if (!oldFile.exists())
    {
        return false; // file does not exist
    }

    // Get the time the file was created
    // TODO: Is this the local time or the UTC time stamp??
    Poco::DateTime creationTime = Poco::DateTime(oldFile.created());

    // build the new file-name
    int lastDot = fileName.ReverseFind('.');
    if (lastDot == -1)
    {
        newFileName.Format("%s_%04d%02d%02d_%02d%02d", (const char*)fileName,
            creationTime.year(), creationTime.month(), creationTime.day(), creationTime.hour(), creationTime.minute());
    }
    else
    {
        newFileName.Format("%s_%04d%02d%02d_%02d%02d%s", (const char*)fileName.Left(lastDot),
            creationTime.year(), creationTime.month(), creationTime.day(), creationTime.hour(), creationTime.minute(), (const char*)fileName.Right(fileName.GetLength() - lastDot));
    }

    // move the file
    oldFile.moveTo(newFileName.std_str());

    return true;
}

