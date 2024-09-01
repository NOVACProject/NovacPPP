#include <PPPLib/Configuration/EvaluationConfiguration.h>
#include <SpectralEvaluation/StringUtils.h>
#include <sstream>

using namespace Configuration;

CEvaluationConfiguration::CEvaluationConfiguration(void)
{
}

CEvaluationConfiguration::~CEvaluationConfiguration(void)
{
    Clear();
}

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

int CEvaluationConfiguration::SetFitWindow(int index, const novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo)
{
    if (index < 0)
    {
        return 1;
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

int CEvaluationConfiguration::GetFitWindow(int index, novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo) const
{
    if (index < 0 || index >= m_windows.size())
        return 1;

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
    int nWindows = static_cast<int>(m_windows.size());
    if (nWindows == 0)
    {
        throw std::invalid_argument("No fit window defined");
    }

    // Check the time ranges
    for (int k = 0; k < nWindows; ++k)
    {
        if (m_windows[k].validFrom >= m_windows[k].validTo)
        {
            std::stringstream msg;
            msg << "Evaluation configuration for " << this->m_serial << " has invalid time range: " << m_windows[k];
            throw std::invalid_argument(msg.str());
        }

        // check if this time range overlaps some other 
        for (int j = k + 1; j < nWindows; ++j)
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