#include "StdAfx.h"
#include "XMLFileReader.h"
#include "../Common/Common.h"

using namespace FileHandler;

CXMLFileReader::CXMLFileReader(void)
{
	this->m_File = NULL;
	this->nLinesRead = 0;
	this->szToken = NULL;
	this->m_filename.Format("");
	
	m_tokenPt = NULL;
}

CXMLFileReader::~CXMLFileReader(void)
{
	m_tokenPt = NULL;
}

void CXMLFileReader::SetFile(novac::CStdioFile *file)
{
	m_File = file;
}

char *CXMLFileReader::NextToken(){
	char separators[] = "<>\t";
	szToken = NULL;

	if(nLinesRead == 0){
		// if this is the first call to this function
		m_File->ReadString(szLine, 4095);
		szToken = (char*)(LPCSTR)szLine;
	  
		m_tokenPt = strtok(szToken, separators);
		++nLinesRead;
		return m_tokenPt;
	}else{
		// this is not the first call to this function
		if(m_tokenPt = strtok(szToken, separators)){
			return m_tokenPt;
		}

		szLine[0] = 0;
		while(0 == strlen(szLine)){
		if(!m_File->ReadString(szLine, 4095))
			return NULL;
		}

		szToken = (char*)(LPCSTR)szLine;

		m_tokenPt = strtok(szToken, separators);
		++nLinesRead;
	}
	
	if(strlen(m_tokenPt) < 2)
		return NextToken();
	else
		return m_tokenPt;
}

/** Retrieves the value of the given attribute from the current token 
	@return NULL if this attribute does not exist or the current token is not
		a valid element */
const char *CXMLFileReader::GetAttributeValue(const novac::CString &label){
	novac::CString toSearchFor;

	if(nLinesRead == 0 || szToken == NULL)
		return NULL; // we haven't started reading the file yet...
		
	// make a copy of szToken, so that we can make changes to it...
	char copyOfToken[4096];
	sprintf(copyOfToken, "%s", szToken);

	// search for the attribute in szToken
	toSearchFor.Format("%s=\"", label);
	char *pt_start = strstr(copyOfToken, toSearchFor);
	if(pt_start == NULL)
		return NULL;
	pt_start += toSearchFor.GetLength(); // point to the character after the double-quote
	
	// search for the ending double-quote
	char *pt_end = strstr(pt_start, "\"");
	if(pt_end == NULL)
		return NULL;
	*pt_end = '\0'; // make the string end at the position of the double-quote
	
	// now we have the value of the attribute	
	sprintf(attributeValue, "%s", pt_start);
	
	return attributeValue;
}

/** General parsing of a single, simple string item */
int CXMLFileReader::Parse_StringItem(const novac::CString &label, novac::CString &string){
	string.Format("");

	while(szToken = NextToken()){

		if(Equals(szToken, label))
			return 1;

		string.Format(szToken);
	}

	return 0;
}
int CXMLFileReader::Parse_LongItem(const novac::CString &label, long &number){

	while(szToken = NextToken()){

		if(Equals(szToken, label))
		return 1;

	   
		number = _tstol(szToken);
	}

	return 0;
}
/** General parsing of a single, simple float item */
int CXMLFileReader::Parse_FloatItem(const novac::CString &label, double &number){
	while(szToken = NextToken()){

		if(Equals(szToken, label))
		return 1;

		number = _tstof(szToken);
	}

	return 0;
}

/** General parsing of a single, simple integer item */
int CXMLFileReader::Parse_IntItem(const novac::CString &label, int &number){
	while(szToken = NextToken()){

		if(Equals(szToken, label))
		return 1;

	   
		number = _tstoi(szToken);
	}

	return 0;
}

/** General parsing of a single, simple long integer item */
int CXMLFileReader::Parse_IPNumber(const novac::CString &label, BYTE &ip0, BYTE &ip1, BYTE &ip2, BYTE &ip3){
	while(szToken = NextToken()){
		int i0, i1, i2,i3;

		if(Equals(szToken, label))
		return 1;

		sscanf(szToken, "%d.%d.%d.%d", &i0, &i1, &i2, &i3);
		ip0 = i0;
		ip1 = i1;
		ip2 = i2;
		ip3 = i3;
	}

	return 0;
}
/** General parsing of a date */
int CXMLFileReader::Parse_Date(const novac::CString &label, CDateTime &datum){
	int nFields = 0;

	while(szToken = NextToken()){
		int i0 = 0;
		int i1 = 0;
		int i2 = 0;
		int i3 = 0;
		int i4 = 0;
		int i5 = 0;

		if(Equals(szToken, label))
			return 1;

		char *pt = strstr(szToken, "T");

		if(pt == NULL){
			nFields = sscanf(szToken, "%d.%d.%d", &i0, &i1, &i2);
			datum.year = i0;
			datum.month = i1;
			datum.day = i2;
		}else{
			nFields = sscanf(szToken, "%d.%d.%dT%d:%d:%d", &i0, &i1, &i2, &i3, &i4, &i5);
			datum.year = i0;
			datum.month = i1;
			datum.day = i2;
			datum.hour = i3;
			datum.minute = i4;
			datum.second = i5;
		}
		
		if(nFields == 0){
			// if the normal parsing didn't work, then try also to parse functional expressions...
			CDateTime::ParseDate(szToken, datum);
		}
	}

	return 0;
}