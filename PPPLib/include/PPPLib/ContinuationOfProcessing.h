#pragma once

#include <string>
#include <vector>

namespace Configuration
{
class CUserConfiguration;
}

// #include <afxtempl.h>

/** The class <b>CContinuationOfProcessing</b> is used to keep track of
    what has already been done when continuing an old processing run.
    This class is only used if userSettings.m_fIsContinuation == true */
class CContinuationOfProcessing
{
public:
    CContinuationOfProcessing(const Configuration::CUserConfiguration& userSettings);

    /** This returns true if the given pak-file has previously been
        processed at an earlier processing round. */
    bool IsPreviouslyIgnored(const std::string& pakFileName) const;

private:

    /** if userSettings.m_fIsContinuation == true then this will scan through an old
            StatusLog file (if found) for names of .pak-files that have been processed
            at an earlier point.
        if userSettings.m_fIsContinuation == false then this will just return without
            doing anything.	 */
    void ScanStatusLogFileForOldScans(const Configuration::CUserConfiguration& userSettings);

    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    /** This is used to store the names of the .pak-files that we have already processed and ignored */
    std::vector<std::string> m_previouslyIgnoredFiles;
};
