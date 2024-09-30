#pragma once

#include <SpectralEvaluation/DateTime.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/MFC/CStdioFile.h>
#include <PPPLib/Logging.h>

namespace FileHandler
{
class CXMLFileReader
{
public:
    CXMLFileReader(novac::ILogger& log);
    virtual ~CXMLFileReader();

    // Non copyable object, since we are managing a file pointer
    CXMLFileReader(const CXMLFileReader&) = delete;
    CXMLFileReader& operator=(const CXMLFileReader&) = delete;

    /** Retrieve the next token from the xml file.
        Returns nullptr if no more tokens are available.  */
    char* NextToken();

    /** Retrieves the value of the given attribute from the current token
        @return NULL if this attribute does not exist or the current token is not
            a valid element */
    const char* GetAttributeValue(const novac::CString& label);

    /** General parsing of a single, simple string item */
    int Parse_StringItem(const novac::CString& label, novac::CString& string);
    int Parse_StringItem(const novac::CString& label, std::string& string);

    /** General parsing of a single, simple path.
        This both parses the path and converts it to a full, absolute path by expanding environment variables */
    int Parse_PathItem(const novac::CString& label, novac::CString& path);
    int Parse_PathItem(const novac::CString& label, std::string& path);

    /** General parsing of a single, simple float item */
    int Parse_FloatItem(const novac::CString& label, double& number);

    /** General parsing of a single, simple integer item */
    int Parse_IntItem(const novac::CString& label, int& number);

    /** General parsing of a single, simple long integer item */
    int Parse_LongItem(const novac::CString& label, long& number);

    /** General parsing of a single, simple long integer item */
    int Parse_IPNumber(const novac::CString& label, std::uint8_t& ip0, std::uint8_t& ip1, std::uint8_t& ip2, std::uint8_t& ip3);

    /** General parsing of a date */
    int Parse_Date(const novac::CString& label, novac::CDateTime& datum);

protected:
    novac::ILogger& m_log;

    /** The tokenizer */
    char* szToken = nullptr;

    /** The name of the currently opened file. For debugging reasons */
    novac::CString m_filename = "";

    /** Opens the provided file for reading.
        @return true if successful */
    bool Open(const novac::CString& fileName);

    /** Closes m_File */
    void Close();

    /** Returns true if the provided token equals the closing xml tag */
    static bool IsClosingTag(const novac::CString& endTag, const char* token);

private:
    /** Pointer to the next token. Should only be modified by 'NextToken()' */
    char* m_tokenPt = nullptr;

    /** A handle to the file to read from. */
    novac::CStdioFile* m_File = nullptr;

    /** The string that was read from the file */
    char szLine[4096];

    /** String representing the value of the last retrieved attribute */
    char attributeValue[4096];

    /** The number of lines that has been read from the file */
    long nLinesRead = 0;
};
}