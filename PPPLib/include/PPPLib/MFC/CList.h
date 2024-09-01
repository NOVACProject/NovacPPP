#ifndef NOVAC_PPPLIB_CLIST_H
#define NOVAC_PPPLIB_CLIST_H

#include <list>
#include <cassert>

namespace novac
{
template<class TYPE>
struct POSITION
{
public:
    POSITION()
        : m_emptyList(), m_data(&m_emptyList)
    {
        m_position = m_data->begin();
    }

    POSITION(std::list<TYPE>& data)
        : m_emptyList(), m_data(&data)
    {
        m_position = m_data->begin();
    }

    POSITION(const TYPE*)
        : m_emptyList(), m_data(&m_emptyList)
    {
        m_position = m_data->begin();
    }

    POSITION(const POSITION& other)
        : m_emptyList(), m_data(other.m_data)
    {
        m_position = other.m_position;
    }

    POSITION& operator=(const POSITION& other)
    {
        this->m_data = other.m_data;
        this->m_position = other.m_position;

        return *this;
    }

    POSITION& operator=(const TYPE*)
    {
        this->m_data = &m_emptyList;
        this->m_position = m_emptyList.begin();

        return *this;
    }

    bool HasNext() const
    {
        if (nullptr == m_data || 0 == m_data->size())
        {
            return false;
        }
        else
        {
            return m_position != m_data->end();
        }
    }

    TYPE& GetNext()
    {
        TYPE& data = *this->m_position;
        ++(this->m_position);
        return data;
    }

    TYPE& GetAt() const
    {
        return *(this->m_position);
    }

    void InsertBefore(TYPE item)
    {
        m_data->insert(m_position, item);
    }

    bool operator==(void* data)
    {
        return (nullptr == data) ? (!this->HasNext()) : false;
    }

    bool operator!=(void* data)
    {
        return (nullptr == data) ? (this->HasNext()) : true;
    }

private:
    typename std::list<TYPE>::iterator m_position;
    std::list<TYPE> m_emptyList;
    std::list<TYPE>* m_data;
};

template<class TYPE>
struct CONST_POSITION : public POSITION<TYPE>
{
public:
    CONST_POSITION()
        : m_emptyList(), m_data(&m_emptyList)
    {
        m_position = m_data->begin();
    }

    CONST_POSITION(const std::list<TYPE>& data)
        : m_emptyList(), m_data(&data)
    {
        m_position = m_data->begin();
    }

    CONST_POSITION(const TYPE*)
        : m_emptyList(), m_data(&m_emptyList)
    {
        m_position = m_data->begin();
    }

    CONST_POSITION(const CONST_POSITION& other)
        : m_emptyList(), m_data(other.m_data)
    {
        m_position = other.m_position;
    }

    POSITION<TYPE>& operator=(const TYPE*)
    {
        this->m_position = m_data->end();
        return *this;
    }

    bool HasNext() const
    {
        if (nullptr == m_data || 0 == m_data->size())
        {
            return false;
        }
        else
        {
            return m_position != m_data->end();
        }
    }

    const TYPE& GetAt() const
    {
        return *(this->m_position);
    }

    const TYPE& GetNext()
    {
        const TYPE& data = *this->m_position;
        ++(this->m_position);
        return data;
    }

    bool operator==(void* data)
    {
        return (nullptr == data) ? (!this->HasNext()) : false;
    }

    bool operator!=(void* data)
    {
        return (nullptr == data) ? (this->HasNext()) : true;
    }

private:
    typename std::list<TYPE>::const_iterator m_position;
    const std::list<TYPE> m_emptyList;
    const std::list<TYPE>* m_data;
};

template<class TYPE>
struct REVERSE_POSITION
{
public:
    REVERSE_POSITION()
        : m_emptyList(), m_data(&m_emptyList)
    {
        m_position = m_data->rbegin();
    }

    REVERSE_POSITION(std::list<TYPE>& data)
        : m_emptyList(), m_data(&data)
    {
        m_position = data.rbegin();
    }

    REVERSE_POSITION(void*)
        : m_emptyList(), m_data(&m_emptyList)
    {
        m_position = m_data->rbegin();
    }

    bool HasPrevious() const
    {
        if (m_data->size() == 0)
        {
            return false;
        }
        else
        {
            return m_position != m_data->rend();
        }
    }

    TYPE& GetAt()
    {
        return *(this->m_position);
    }

    TYPE& GetPrev()
    {
        TYPE& data = *m_position;
        ++(m_position);
        return data;
    }

    bool operator==(void* data)
    {
        return (nullptr == data) ? (!this->HasPrevious()) : false;
    }

    bool operator!=(void* data)
    {
        return (nullptr == data) ? (this->HasPrevious()) : true;
    }

private:
    typename std::list<TYPE>::reverse_iterator  m_position;
    std::list<TYPE> m_emptyList;
    std::list<TYPE>* m_data;
};

template<class TYPE>
struct CONST_REVERSE_POSITION
{
public:
    CONST_REVERSE_POSITION()
        : m_emptyList(), m_data(&m_emptyList)
    {
        m_position = m_data->crbegin();
    }

    CONST_REVERSE_POSITION(const std::list<TYPE>& data)
        : m_emptyList(), m_data(data)
    {
        m_position = data.crbegin();
    }

    CONST_REVERSE_POSITION(void*)
        : m_emptyList(), m_data(&m_emptyList)
    {
        m_position = m_data->crbegin();
    }

    bool HasPrevious() const
    {
        if (m_data->size() == 0)
        {
            return false;
        }
        else
        {
            return m_position != m_data->crend();
        }
    }

    TYPE& GetAt()
    {
        return *(this->m_position);
    }

    const TYPE& GetPrev() const
    {
        const TYPE& data = *m_position;
        ++(m_position);
        return data;
    }
    bool operator==(void* data)
    {
        return (nullptr == data) ? (!this->HasPrevious()) : false;
    }

    bool operator!=(void* data)
    {
        return (nullptr == data) ? (this->HasPrevious()) : true;
    }

private:
    const typename std::list<TYPE>::const_reverse_iterator  m_position;
    const std::list<TYPE> m_emptyList;
    const std::list<TYPE>* m_data;
};

template<class TYPE, class ARG_TYPE = const TYPE&>
class CList
{
public:

    // ---------------------- Construction -----------------------
    CList()
    {
    }

    ~CList()
    {
    }

    int GetSize() const
    {
        return (int)m_data.size();
    }

    int GetCount() const
    {
        return (int)m_data.size();
    }

    POSITION<TYPE> GetHeadPosition()
    {
        POSITION<TYPE> p(m_data);
        return p;
    }

    const CONST_POSITION<TYPE> GetHeadPosition() const
    {
        CONST_POSITION<TYPE> p(m_data);
        return p;
    }

    REVERSE_POSITION<TYPE> GetTailPosition()
    {
        REVERSE_POSITION<TYPE> p(m_data);
        return p;
    }

    ARG_TYPE GetAt(POSITION<TYPE>& p) const
    {
        return p.GetAt();
    }

    ARG_TYPE GetAt(REVERSE_POSITION<TYPE>& p) const
    {
        return p.GetAt();
    }

    const ARG_TYPE GetAt(CONST_POSITION<TYPE>& p) const
    {
        return p.GetAt();
    }

    ARG_TYPE GetNext(POSITION<TYPE>& p) const
    {
        return p.GetNext();
    }

    const ARG_TYPE GetNext(CONST_POSITION<TYPE>& p) const
    {
        return p.GetNext();
    }

    ARG_TYPE GetPrev(REVERSE_POSITION<TYPE>& p) const
    {
        return p.GetPrev();
    }

    ARG_TYPE GetPrev(CONST_REVERSE_POSITION<TYPE>& p) const
    {
        return p.GetPrev();
    }

    // ---------------------- Operations -----------------------

    void RemoveAll()
    {
        m_data.clear();
    }

    void RemoveTail()
    {
        m_data.erase(--m_data.end());
    }

    void AddTail(TYPE item)
    {
        m_data.push_back(item);
    }

    void AddHead(TYPE item)
    {
        m_data.push_front(item);
    }

    void InsertBefore(POSITION<TYPE>& pos, TYPE item)
    {
        pos.InsertBefore(item);
    }

private:
    std::list<TYPE> m_data;
};
}

#endif // !NOVAC_PPPLIB_CLIST_H
