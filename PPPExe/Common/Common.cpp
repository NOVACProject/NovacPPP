#include "Common.h"

#include <Poco/File.h>
#include <Poco/DateTime.h>
#include <Poco/Message.h>
#include <Poco/Logger.h>

#include <sstream>

extern std::string s_exePath;
extern std::string s_exeFileName;

#undef min
#undef max

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
void PocoLogger::Debug(const novac::LogContext& c, const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    if (log.debug())
    {
        std::stringstream s;
        s << c << message;
        log.debug(s.str());
    }
}

void PocoLogger::Information(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.information(message);
}
void PocoLogger::Information(const novac::LogContext& c, const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    if (log.information())
    {
        std::stringstream s;
        s << c << message;
        log.information(s.str());
    }
}

void PocoLogger::Error(const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    log.fatal(message);
}
void PocoLogger::Error(const novac::LogContext& c, const std::string& message)
{
    Poco::Logger& log = Poco::Logger::get("NovacPPP");
    if (log.error())
    {
        std::stringstream s;
        s << c << message;
        log.error(s.str());
    }
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

void Common::CopyFile(const novac::CString& oldName, const novac::CString& newName)
{
    Poco::File oldFile(oldName.std_str());

    oldFile.copyTo(newName.std_str());
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

