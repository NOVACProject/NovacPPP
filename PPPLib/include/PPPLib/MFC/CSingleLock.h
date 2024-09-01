#ifndef NOVAC_PPPLIB_CSINGLE_LOCK_H
#define NOVAC_PPPLIB_CSINGLE_LOCK_H

#include "CCriticalSection.h"
#include <memory>

namespace novac
{
class CSingleLock
{
public:

    // ---------------------- Construction -----------------------
    CSingleLock()
        : m_section(nullptr)
    {
    }

    CSingleLock(CCriticalSection* section)
        : m_section(section)
    {
    }

    ~CSingleLock()
    {
        m_lock.reset();
        m_section = nullptr;
    }

    void Lock()
    {
        m_lock = std::make_unique<std::lock_guard<std::mutex>>(m_section->m_mutex);
    }

    void Unlock()
    {
        m_lock.reset();
    }

    bool IsLocked() const
    {
        if (m_lock) { return true; }
        else { return false; }
    }

private:
    std::unique_ptr<std::lock_guard<std::mutex>> m_lock;
    CCriticalSection* m_section;
};
}

#endif // !NOVAC_PPPLIB_CSINGLE_LOCK_H
