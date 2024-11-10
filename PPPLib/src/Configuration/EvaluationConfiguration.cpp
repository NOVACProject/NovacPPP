#include <PPPLib/Configuration/EvaluationConfiguration.h>
#include <PPPLib/Definitions.h>
#include <SpectralEvaluation/StringUtils.h>
#include <sstream>

namespace Configuration
{

void CEvaluationConfiguration::Clear()
{
    m_windows.clear();
}

void CEvaluationConfiguration::InsertFitWindow(const novac::CFitWindow& window, const novac::CDateTime& validFrom, const novac::CDateTime& validTo)
{
    FitWindowWithTime timedWindow;
    timedWindow.window = window;
    timedWindow.validFrom = validFrom;
    timedWindow.validTo = validTo;

    m_windows.push_back(timedWindow);
}

int CEvaluationConfiguration::SetFitWindow(size_t index, const novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo)
{
    if (index > MAX_FIT_WINDOWS)
    {
        throw std::invalid_argument("Invalid fit window index, cannot create more than MAX_FIT_WINDOWS windows.");
    }

    if (index >= m_windows.size())
    {
        m_windows.resize(1 + index);
    }

    FitWindowWithTime timedWindow;
    timedWindow.window = window;
    timedWindow.validFrom = validFrom;
    timedWindow.validTo = validTo;
    m_windows[index] = std::move(timedWindow);

    return 0;
}

int CEvaluationConfiguration::GetFitWindow(size_t index, novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo) const
{
    if (index >= m_windows.size())
    {
        return 1;
    }

    window = m_windows[index].window;
    validFrom = m_windows[index].validFrom;
    validTo = m_windows[index].validTo;

    return 0;
}

std::ostream& operator << (std::ostream& out, const FitWindowWithTime& eval)
{
    out << "[" << eval.window.name << "]: " << eval.validFrom << " to " << eval.validTo;
    return out;
}

void  CEvaluationConfiguration::CheckSettings() const
{

    // make sure that at least one fit-window is defined
    const size_t nWindows = m_windows.size();
    if (nWindows == 0)
    {
        throw std::invalid_argument("No fit window defined");
    }

    // Check the time ranges
    for (size_t k = 0; k < nWindows; ++k)
    {
        if (m_windows[k].validFrom >= m_windows[k].validTo)
        {
            std::stringstream msg;
            msg << "Evaluation configuration for " << this->m_serial << " has invalid time range: " << m_windows[k];
            throw std::invalid_argument(msg.str());
        }
        if (m_windows[k].window.offsetRemovalRange.from > m_windows[k].window.offsetRemovalRange.to)
        {
            std::stringstream msg;
            msg << "Fit window for " << this->m_serial << " has invalid offset removal range: " << m_windows[k];
            throw std::invalid_argument(msg.str());
        }

        // check if this time range overlaps some other 
        for (size_t j = k + 1; j < nWindows; ++j)
        {
            if (m_windows[j].window.channel != m_windows[k].window.channel)
            {
                continue; // no use to compare master and slave...
            }
            else if (!EqualsIgnoringCase(m_windows[j].window.name, m_windows[k].window.name))
            {
                continue; // if the windows have different names, then don't compare...
            }
            else
            {
                if ((m_windows[k].validFrom < m_windows[j].validFrom) && (m_windows[k].validTo > m_windows[j].validFrom))
                {
                    std::stringstream msg;
                    msg << "Evaluation configuration for " << this->m_serial << " has overlapping time ranges: (" << m_windows[k] << ") and (" << m_windows[j] << ")";
                    throw std::invalid_argument(msg.str());
                }
                else if ((m_windows[j].validFrom < m_windows[k].validFrom) && (m_windows[j].validTo > m_windows[k].validFrom))
                {
                    std::stringstream msg;
                    msg << "Evaluation configuration for " << this->m_serial << " has overlapping time ranges: (" << m_windows[k] << ") and (" << m_windows[j] << ")";
                    throw std::invalid_argument(msg.str());
                }
            }
        }
    }

    return;
}

}  // namespace Configuration