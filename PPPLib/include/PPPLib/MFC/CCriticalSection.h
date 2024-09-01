#ifndef NOVAC_PPPLIB_CRITICALSECTION_H
#define NOVAC_PPPLIB_CRITICALSECTION_H

#include <mutex>

namespace novac
{
// Emulating a MFC CriticalSection
class CCriticalSection
{
public:

    CCriticalSection()
    {

    }

    ~CCriticalSection()
    {
    }

    std::mutex m_mutex;
};
}

#endif // !NOVAC_PPPLIB_CRITICALSECTION_H
