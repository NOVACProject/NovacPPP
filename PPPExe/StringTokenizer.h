#pragma once

#include <stddef.h>

/** The class CStringTokenizer takes care of tokenizing strings the way we want them
    This is primarily used when passing parameters to the program through the command-line
        since the normal CRT tokenizer will not be able to keep together strings containing spaces
    This class will return any string surrounded by ":s as one token (even if it contains spaces)
*/

#define MAX_N_SEPARATORS 256

class CStringTokenizer
{
public:
    CStringTokenizer(void);
    CStringTokenizer(const char* string, const char* separators);
    ~CStringTokenizer(void);

    /** Returns the next token */
    const char* NextToken();

private:
    /** The string that we are to interpret */
    char* m_stringToInterpret;

    /** The length of the string to interpret */
    size_t	m_strLen;

    /** Out current position in the string to interpret */
    size_t	m_stringPos;

    /** Local buffer holding the last token returned */
    char* m_lastToken;

    /** The separator characters */
    char	m_separators[MAX_N_SEPARATORS];
    int		m_separatorNum;
};
