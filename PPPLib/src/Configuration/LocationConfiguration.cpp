#include <PPPLib/Configuration/LocationConfiguration.h>
#include <sstream>

namespace Configuration
{
    std::ostream& operator << (std::ostream& out, const CInstrumentLocation& location)
    {
        out << "[" << location.m_locationName << "]: " << location.m_validFrom << " to " << location.m_validTo;
        return out;
    }

    void CLocationConfiguration::Clear()
    {
        m_locationNum = 0;
        for (int k = 0; k < MAX_N_LOCATIONS; ++k)
        {
            this->m_location[k].Clear();
        }
    }

    void CLocationConfiguration::InsertLocation(const CInstrumentLocation& loc)
    {
        if (m_locationNum == MAX_N_LOCATIONS)
        {
            return;
        }

        // insert the location into the array
        m_location[m_locationNum] = loc;

        // increase the counter for the number of locationss
        ++m_locationNum;

        return;
    }

    int CLocationConfiguration::GetLocation(int index, CInstrumentLocation& loc) const
    {
        if (index < 0 || index >= m_locationNum)
            return 1;

        // copy the data to the requested parameter
        loc = m_location[index];

        return 0;
    }

    unsigned long CLocationConfiguration::GetLocationNum() const
    {
        return m_locationNum;
    }

    void CLocationConfiguration::CheckSettings() const
    {
        if (m_locationNum == 0)
        {
            throw std::invalid_argument("No location defined");
        }

        // Check the time ranges
        for (int k = 0; k < m_locationNum; ++k)
        {
            // check that the time range is valid
            if (m_location[k].m_validFrom >= m_location[k].m_validTo) 
            {
                std::stringstream msg;
                msg << "Time range number " << k << " is invalid. Valid from (" << m_location[k].m_validFrom << ") is is not before Valid to (" << m_location[k].m_validTo << ")";
                throw std::invalid_argument(msg.str());
            }
            
            // check if this time range overlaps some other 
            for (int j = k + 1; j < m_locationNum; ++j) 
            {
                if ((m_location[k].m_validFrom < m_location[j].m_validFrom) && (m_location[k].m_validTo > m_location[j].m_validFrom)) 
                {
                    std::stringstream msg;
                    msg << "Time ranges: (" << m_location[k] << ") and (" << m_location[j] << ") are overlapping.";
                    throw std::invalid_argument(msg.str());
                }
                else if ((m_location[j].m_validFrom < m_location[k].m_validFrom) && (m_location[j].m_validTo > m_location[k].m_validFrom)) 
                {
                    std::stringstream msg;
                    msg << "Time ranges: (" << m_location[k] << ") and (" << m_location[j] << ") are overlapping.";
                    throw std::invalid_argument(msg.str());
                }
            }
        }
    }
}
