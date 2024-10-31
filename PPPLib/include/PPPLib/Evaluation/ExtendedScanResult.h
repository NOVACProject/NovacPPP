#pragma once

#include <PPPLib/Definitions.h>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/NovacEnums.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>

namespace Evaluation
{
/** The struct <b>CExtendedScanResult</b> is a container used to
    store the result of the evaluation of one .pak - file using several
    fit-windows.

    This is used in the PostProcessing, each scan can be evaluated in different
        wavelength regions and using different settings for the fit. The results
        in these different wavelength regions can then be compared for other purposes
        such as studying the radiative transfer...

    To save memory does	this not contain the column datas, only the names of the .txt files
        where they can be found.

    To reduce the number of times we need to read data from the evaluation-log files
        this also contains the calculated properties of the plume in the scan, such as the
        completeness, the angle of centre of mass of the plume etc.
*/
struct CExtendedScanResult
{
public:

    CExtendedScanResult()
    {
        m_evalLogFile.reserve(MAX_FIT_WINDOWS);
        m_fitWindowName.reserve(MAX_FIT_WINDOWS);
    }

    CExtendedScanResult(std::string serial, novac::CDateTime scanStartTime, novac::MeasurementMode mode)
        : m_instrumentSerial(serial), m_measurementMode(mode), m_startTime(scanStartTime)
    {
        m_evalLogFile.reserve(MAX_FIT_WINDOWS);
        m_fitWindowName.reserve(MAX_FIT_WINDOWS);
    }

    std::string m_instrumentSerial = "";

    novac::MeasurementMode m_measurementMode = novac::MeasurementMode::Flux;

    /** The full path and file-name to the .pak-file containing the spectra
        from which this result was computed. */
    std::string m_pakFile = "";

    /** The full path and file-name of the evaluation log file that was generated.
        The file itself contains the result of the evaluation */
    std::vector<std::string> m_evalLogFile;

    /** The full name of the fit-window that was used to generate each evaluation
        result, this is typically a name such as 'SO2', 'SO2_low' or 'O4'
         m_fitWindowName[i] is the name of the fit-window behind the result in m_evalLogFile[i]. */
    std::vector<std::string> m_fitWindowName;

    /** The date and time that the scan was generated. */
    novac::CDateTime m_startTime;

    /** The properties of this scan. This is only evaluated in the main-fit window */
    novac::CPlumeInScanProperty m_scanProperties;

};
}
