#include "StdAfx.h"
#include "ContinuationOfProcessing.h"

// This is the settings for how to do the procesing
#include "Configuration/UserConfiguration.h"

extern Configuration::CUserConfiguration			g_userSettings;// <-- The settings of the user

CContinuationOfProcessing					g_continuation;		// <-- Information on what has already been done when continuing an old processing round

CContinuationOfProcessing::CContinuationOfProcessing(void)
{
	ScanStatusLogFileForOldScans();
}

CContinuationOfProcessing::~CContinuationOfProcessing(void)
{
}

/** if g_userSettings.m_fIsContinuation == true then this will scan through an old 
		StatusLog file (if found) for names of .pak-files that have been processed 
		at an earlier point.
	if g_userSettings.m_fIsContinuation == false then this will just return without 
		doing anything.	 */
void CContinuationOfProcessing::ScanStatusLogFileForOldScans(){
	CString oldStatusLogfile;
	CString fileName;

	m_previouslyIgnoredFiles.RemoveAll();

	if(g_userSettings.m_fIsContinuation == false)
		return;

	oldStatusLogfile.Format("%s\\StatusLog.txt", g_userSettings.m_outputDirectory);
	if(!IsExistingFile(oldStatusLogfile))
		return;

	FILE *f = fopen(oldStatusLogfile, "r");
	if(f == NULL){
		return;
	}

	char *buffer = new char[16384];

	while(NULL != fgets(buffer, 16383, f)){
		char *pt = strstr(buffer, " does not see the plume");
		if(NULL != pt){
			// if this line corresponds to an ignored scan
			pt[0] = '\0';
			fileName.Format("%s", buffer + 8);
			m_previouslyIgnoredFiles.AddTail(fileName);
			continue;
		}	
	}
	fclose(f);
}

bool CContinuationOfProcessing::IsPreviouslyIgnored(const CString &pakFileName){
	POSITION p = m_previouslyIgnoredFiles.GetHeadPosition();
	while(p != NULL){
		CString &fileName = m_previouslyIgnoredFiles.GetNext(p);
		if(Equals(fileName, pakFileName))
			return true;
	}
	return false; // not found
}