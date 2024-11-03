#include "PPPLib/MFC/CStringTokenizer.h"
#include <set>

namespace novac
{
CStringTokenizer::CStringTokenizer(const char* text, const char* separators)
    : m_data(text), m_curToken(0)
{
    std::set<char> separatorSet;
    const char* pt = separators;
    while (*pt != 0)
    {
        separatorSet.insert(*pt);
        ++pt;
    }

    // split the input text into tokens
    auto pt1 = m_data.cbegin();
    auto pt2 = m_data.cbegin();
    while (pt1 != m_data.cend())
    {
        // remove separators in the beginning...
        while (pt2 != m_data.cend() && separatorSet.find(*pt2) != separatorSet.end())
        {
            ++pt2;
        }

        // go to the next separator char
        while (pt2 != m_data.cend() && separatorSet.find(*pt2) == separatorSet.end())
        {
            ++pt2;
        }

        m_tokens.push_back(std::string(pt1, pt2));

        // remove separators in the end...
        while (pt2 != m_data.cend() && separatorSet.find(*pt2) != separatorSet.end())
        {
            ++pt2;
        }
        pt1 = pt2;
    }

}

const char* CStringTokenizer::NextToken()
{
    if (m_curToken >= m_tokens.size())
    {
        return nullptr;
    }
    else
    {
        return m_tokens[m_curToken++].c_str();
    }
}
}
