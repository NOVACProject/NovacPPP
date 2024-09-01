#include "stdafx.h"
#include "StringTokenizer.h"
#include <string>
#include <algorithm>
#include <cassert>
#include <string.h>

#undef min
#undef max

CStringTokenizer::CStringTokenizer(void)
{
    m_stringToInterpret = nullptr;
    m_strLen = 0;
    this->m_stringPos = 0;
    memset(m_separators, 0, MAX_N_SEPARATORS * sizeof(char));
    m_separatorNum = 0;
    m_lastToken = nullptr;
}

CStringTokenizer::CStringTokenizer(const char* string, const char* separators)
{
    if (string != nullptr)
    {
        m_strLen = (long)strlen(string);

        // allocate a new string as a copy of the string to tokenize
        this->m_stringToInterpret = new char[m_strLen + 2];
        sprintf(m_stringToInterpret, "%s", string);
    }
    else
    {
        m_stringToInterpret = nullptr;
        m_strLen = 0;
    }

    memset(m_separators, 0, MAX_N_SEPARATORS * sizeof(char));
    if (separators != nullptr)
    {
        m_separatorNum = std::min(MAX_N_SEPARATORS, (int)strlen(separators));
        memcpy(m_separators, separators, m_separatorNum * sizeof(char));
    }
    else
    {
        m_separators[0] = ' ';
        m_separators[1] = '\t';
        m_separators[2] = '\n';
        m_separatorNum = 3;
    }

    m_lastToken = nullptr;
    m_stringPos = 0;
}

CStringTokenizer::~CStringTokenizer(void)
{
    if (m_stringToInterpret != nullptr)
    {
        delete m_stringToInterpret;
        m_stringToInterpret = nullptr;
    }
    if (m_lastToken != nullptr)
    {
        delete m_lastToken;
        m_lastToken = nullptr;
    }
    m_stringPos = 0;
    m_strLen = 0;
}

const char* CStringTokenizer::NextToken()
{
    // check if we're out of the string
    if (m_stringPos >= m_strLen)
    {
        return nullptr;
    }
    if (m_lastToken != nullptr)
    {
        delete m_lastToken;
        m_lastToken = nullptr;
    }

    size_t from = m_stringPos;

    // find the next quotation mark (this is nullptr if none found)
    size_t nextQuote = m_strLen;
    char* pt = strchr(m_stringToInterpret + from, '\"');
    if (pt != nullptr)
    {
        nextQuote = (size_t)(pt - m_stringToInterpret);
    }

    // find the next separator character
    size_t nextSep = m_strLen;
    for (int k = 0; k < m_separatorNum; ++k)
    {
        char* strPt = strchr(m_stringToInterpret + from, m_separators[k]);
        if (strPt != nullptr)
        {
            nextSep = std::min(nextSep, size_t(strPt - m_stringToInterpret));
        }
    }

    // check if the next separator is within a quotation
    if (nextQuote != m_strLen && (nextQuote < nextSep))
    {
        size_t previousQuote = nextQuote;
        const char* strPt = strchr(m_stringToInterpret + previousQuote + 1, '\"');
        if (strPt != nullptr)
        {
            nextQuote = (size_t)(strPt - m_stringToInterpret);
        }

        nextSep = nextQuote + 1;
        for (int k = 0; k < m_separatorNum; ++k)
        {
            const char* strPt2 = strchr(m_stringToInterpret + nextQuote, m_separators[k]);
            if (strPt2 != nullptr)
            {
                nextSep = std::min(nextSep, size_t(strPt2 - m_stringToInterpret));
            }
        }
    }

    m_lastToken = new char[nextSep - m_stringPos + 5]; // allocate a few characters extra, to be safe...
    memcpy(m_lastToken, m_stringToInterpret + from, (nextSep - m_stringPos) * sizeof(char));
    m_lastToken[nextSep - m_stringPos] = '\0';
    size_t N = strlen(m_lastToken);
    assert(N <= (nextSep - m_stringPos + 1));

    // if there's any quotation signs inside this token then remove them...
    while (nullptr != (pt = strchr(m_lastToken, '"')))
    {
        char* last = strchr(m_lastToken, '\0');
        while (pt < last)
        {
            pt[0] = pt[1];
            ++pt;
        }
    }

    // update the counter and return the string
    m_stringPos = nextSep + 1;
    return m_lastToken;
}
