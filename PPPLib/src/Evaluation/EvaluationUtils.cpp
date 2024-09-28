#include <PPPLib/Evaluation/EvaluationUtils.h>
#include <PPPLib/Logging.h>
#include <sstream>

bool Evaluation::IsGoodEnoughToEvaluate(
    const novac::CScanFileHandler& scan,
    const novac::CFitWindow& fitWindow,
    const novac::SpectrometerModel& model,
    const Configuration::CInstrumentLocation& instrLocation,
    const Configuration::CUserConfiguration& userSettings,
    ReasonForScanRejection& reason,
    std::string& reasonMessage)
{
    novac::CSpectrum skySpectrum;

    // Check that the sky-spectrum is ok
    scan.GetSky(skySpectrum);
    if (skySpectrum.IsDark())
    {
        reasonMessage = "Sky spectrum is dark";
        reason = ReasonForScanRejection::SkySpectrumDark;
        return false;
    }

    if ((instrLocation.m_instrumentType == INSTRUMENT_TYPE::INSTR_GOTHENBURG && skySpectrum.ExposureTime() > userSettings.m_maxExposureTime_got) ||
        (instrLocation.m_instrumentType == INSTRUMENT_TYPE::INSTR_HEIDELBERG && skySpectrum.ExposureTime() > userSettings.m_maxExposureTime_hei))
    {
        std::stringstream msg;
        msg << "Sky spectrum has too long exposure time (" << skySpectrum.ExposureTime() << " ms)";
        reasonMessage = msg.str();
        reason = ReasonForScanRejection::SkySpectrumTooLongExposureTime;
        return false;
    }

        const double dynamicRange = skySpectrum.NumSpectra() * model.maximumIntensityForSingleReadout;

    if (skySpectrum.MaxValue(fitWindow.fitLow, fitWindow.fitHigh) >= dynamicRange)
    {
        reasonMessage = "Sky spectrum is saturated in fit region";
        reason = ReasonForScanRejection::SkySpectrumSaturated;
        return false;
    }

    return true;
}