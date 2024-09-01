#ifndef NOVAC_PPPLIB_CARRAY_H
#define NOVAC_PPPLIB_CARRAY_H

#include <vector>

namespace novac
{
template<class TYPE, class ARG_TYPE = const TYPE&>
class CArray
{
public:

    // ---------------------- Construction -----------------------
    CArray()
    {
    }

    ~CArray()
    {
    }

    /** Copies all the elements in the src array here, overwriting the contents */
    void Copy(const CArray& src)
    {
        this->m_data = src.m_data;
    }

    int GetSize() const
    {
        return (int)m_data.size();
    }

    int GetCount() const
    {
        return (int)m_data.size();
    }

    void RemoveAll()
    {
        m_data.clear();
    }

    void FreeExtra()
    {
        // TODO:
    }

    void RemoveAt(int index)
    {
        m_data.erase(m_data.begin() + index);
    }

    void RemoveAt(int index, int number)
    {
        m_data.erase(m_data.begin() + index, m_data.begin() + index + number);
    }

    void SetSize(int newSize)
    {
        m_data.resize(newSize);
    }

    void SetSize(int newSize, int /*growBy*/)
    {
        m_data.resize(newSize);
    }

    TYPE& GetAt(int index)
    {
        return m_data[index];
    }

    const TYPE& GetAt(int index) const
    {
        return m_data[index];
    }

    /* Sets the value at the provided index. The index must be smaller than the size of the CArray. */
    void SetAt(int index, TYPE newValue)
    {
        m_data[index] = newValue;
    }

    /* Sets the value at the provided index. If the index is larger than the number of elements in the CArray
        then the array will grow to accomodate the new item.. */
    void SetAtGrow(int index, TYPE newValue)
    {
        if ((std::size_t)index >= m_data.size())
        {
            this->SetSize(index + 1);
        }
        SetAt(index, newValue);
    }

    // ---------------------- Operations -----------------------

    const TYPE& operator[](int index) const { return m_data[index]; }
    TYPE& operator[](int index) { return m_data[index]; }

private:
    std::vector<TYPE> m_data;
};
}

#endif // !NOVAC_PPPLIB_CARRAY_H
