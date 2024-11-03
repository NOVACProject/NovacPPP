#include <PPPLib/ContinuationOfProcessing.h>
#include <PPPLib/File/Filesystem.h>
#include <SpectralEvaluation/StringUtils.h>

// This is the settings for how to do the procesing
#include <PPPLib/Configuration/UserConfiguration.h>

#include <Poco/Path.h>

#include <string.h>

CContinuationOfProcessing::CContinuationOfProcessing(const Configuration::CUserConfiguration& userSettings)
{
    ScanStatusLogFileForOldScans(userSettings);
}

void CContinuationOfProcessing::ScanStatusLogFileForOldScans(const Configuration::CUserConfiguration& userSettings)
{
    novac::CString oldStatusLogfile;
    novac::CString fileName;

    m_previouslyIgnoredFiles.clear();

    if (userSettings.m_fIsContinuation == false)
        return;

    oldStatusLogfile.Format("%s%cStatusLog.txt", (const char*)userSettings.m_outputDirectory, Poco::Path::separator());
    if (!Filesystem::IsExistingFile(oldStatusLogfile))
        return;

    FILE* f = fopen(oldStatusLogfile, "r");
    if (f == NULL)
    {
        return;
    }

    char* buffer = new char[16384];

    while (NULL != fgets(buffer, 16383, f))
    {
        char* pt = strstr(buffer, " does not see the plume");
        if (NULL != pt)
        {
            // if this line corresponds to an ignored scan
            pt[0] = '\0';
            fileName.Format("%s", buffer + 8);
            m_previouslyIgnoredFiles.push_back(fileName.std_str());
            continue;
        }
    }
    fclose(f);
}

bool CContinuationOfProcessing::IsPreviouslyIgnored(const std::string& pakFileName) const
{
    for (const std::string& fileName : this->m_previouslyIgnoredFiles)
    {
        if (EqualsIgnoringCase(fileName, pakFileName))
        {
            return true;
        }
    }

    return false;
}