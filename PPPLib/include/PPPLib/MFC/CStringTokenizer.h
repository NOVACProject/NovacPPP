#ifndef NOVAC_PPPLIB_CSTRING_TOKENIZER_H
#define NOVAC_PPPLIB_CSTRING_TOKENIZER_H

#include <string>
#include <vector>

namespace novac
{
class CStringTokenizer
{
public:
    CStringTokenizer(const char* text, const char* separators);

    const char* NextToken();

private:

    const std::string m_data;
    std::vector<std::string> m_tokens;
    size_t m_curToken = 0;
};
}  // namespace novac

#endif  // NOVAC_PPPLIB_CSTRING_TOKENIZER_H