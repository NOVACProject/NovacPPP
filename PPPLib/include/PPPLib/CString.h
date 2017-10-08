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

		~CString()
		{
		}

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

		/** Constructs the contents of this string using the common printf formatting. */
		void Format(const char* format, ...);

		/** Appends the contents of this string using the common printf formatting. */
		CString& AppendFormat(const char* format, ...);

		/** Appends the contents of the other string to this. */
		CString& Append(const CString& other);

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

		/** Call this member function to remove instances of ch from the string. Comparisons for the character are case-sensitive.
			@return the number of characters removed. */
		void Remove(char character);

		// ---------------------- Conversion -----------------------

		// explicit conversion to const char*
		operator const char*() const { return m_data.c_str(); }

		// explicit conversion to std::string
		std::string ToStdString() const { return std::string{m_data}; }

		const char* c_str() const { return m_data.c_str(); }

		const std::string std_str() const { return std::string(m_data); }

	private:
		std::string m_data;
	};

	// Appending CStrings
	CString operator+(const CString& str1, CString& other);
	CString operator+(const CString& str1, const char* other);
	CString operator+(const CString& str1, std::string other);

}

#endif // !NOVAC_PPPLIB_CSTRING_H
