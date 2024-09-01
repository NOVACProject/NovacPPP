#ifndef NOVAC_PPPLIB_CSTRING_H
#define NOVAC_PPPLIB_CSTRING_H

#include <string>

namespace novac
{

class CString
{
public:
    // ---------------------- Construction -----------------------
    CString();

    CString(const CString& other);

    CString(const char* other);

    CString(const std::string& other);

    // --------------------- Setting the contents -----------------------

    void SetData(const char* data);
    void SetData(const CString& other);
    void SetData(const std::string& other);

    // --------------------- Simple Operators -----------------------
    bool operator==(const CString& other) const { return this->m_data == other.m_data; }
    bool operator==(const char* other) const { return this->m_data == std::string(other); }
    bool operator!=(const CString& other) const { return this->m_data != other.m_data; }
    bool operator!=(const char* other) const { return this->m_data != std::string(other); }


    // --------------------- Properties -----------------------

    size_t GetLength() const { return m_data.size(); }

    /** @return 0 if this equals the other string */
    int Compare(const CString& other) const;

    // --------------------- Formatting -----------------------


    /** Constructs a new string using the common printf formatting. */
    static CString FormatString(const char* format, ...);

    /** Constructs the contents of this string using the common printf formatting. */
    void Format(const char* format, ...);

    /** Appends the contents of this string using the common printf formatting. */
    CString& AppendFormat(const char* format, ...);

    /** Appends the contents of the other string to this. */
    CString& Append(const CString& other);
    CString& Append(const char* other);
    CString& Append(const std::string& other);

    // ---------------------- Extracting substrings -----------------------

    CString Left(int nChars) const;
    CString Left(size_t nChars) const;
    CString Right(int nChars) const;
    CString Right(size_t nChars) const;

    /** Finds the next token in a target string
        @return A CStringT object containing the current token value. */
    CString Tokenize(const char* delimiters, int& iStart) const;

    /** @return the character at the specified index.
        @throws invalid_argument if index < 0 or index >= GetLength() */
    char GetAt(int index) const;
    char GetAt(size_t index) const;

    // ---------------------- Searching -----------------------

    /** Finds the first occurrence of the given character, -1 if the character is not found. */
    int Find(char ch) const;

    /** Finds the first occurrence of the given character starting at the given position, -1 if the character is not found. */
    int Find(char ch, int pos) const;

    /** Finds the first occurrence of the given string, -1 if the character is not found. */
    int Find(const char* str) const;

    /** Finds the last occurrence of the given character, -1 if the character is not found. */
    int ReverseFind(char ch) const;

    // ---------------------- Changing the String -----------------------

    void Trim();
    void Trim(const char* characters);

    /** Converts all characters in this string to lower-case.
        @return a reference to this object. */
    CString& MakeLower();

    /** Converts all characters in this string to upper-case.
        @return a reference to this object. */
    CString& MakeUpper();

    /** Call this member function to remove instances of ch from the string. Comparisons for the character are case-sensitive. */
    void Remove(char character);

    /** Call this member function to replace instances of the given character in the string with with the other character.
        Comparisons for the character are case-sensitive. */
    void Replace(char character, char replacement);

    // ---------------------- Conversion -----------------------

    // explicit conversion to const char*
    operator const char* () const { return m_data.c_str(); }

    // explicit conversion to std::string
    std::string ToStdString() const { return std::string{ m_data }; }

    const char* c_str() const { return m_data.c_str(); }

    const std::string std_str() const { return std::string(m_data); }

private:
    std::string m_data;
};

// Appending CStrings
CString operator+(const CString& str1, CString& other);
CString operator+(const CString& str1, const char* other);
CString operator+(const CString& str1, std::string other);

// -------------------- Common utility functions for handling CString:s --------------------

/** This function takes a string and removes any 'special' (ASCII code < 32)
    characters in it */
void CleanString(const CString& in, CString& out);
void CleanString(const char* in, CString& out);

/** This function takes a string and simplifies it so that it is more easily readable
    by a machine. changes made are:
    1 - Spaces are replaced with '_' (underscore)
    2 - All characters are converted to lower-case
    3 - accents are removed ('ó' -> 'o') */
void SimplifyString(const CString& in, CString& out);

/** Compares two strings without regard to case.
    @return 1 if the strings are equal. @return 0 if the strings are not equal. */
int Equals(const CString& str1, const CString& str2);

/** Compares two strings without regard to case.
    @return true if the strings are equal, otherwise false. */
bool EqualsIgnoringCase(const std::string& str1, const std::string& str2);

/** Compares at most 'nCharacters' of two strings without regard to case.
    @param nCharacters - The number of characters to compare
    @return 1 if the strings are equal. @return 0 if the strings are not equal. */
int Equals(const CString& str1, const CString& str2, size_t nCharacters);

/** Helper util for extracting the right-most or left-most characters in a std::string */
std::string Right(const std::string& input, size_t nChars);
std::string Left(const std::string& input, size_t nChars);

}  // namespace novac

#endif // !NOVAC_PPPLIB_CSTRING_H
