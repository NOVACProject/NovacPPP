#pragma once

#include <PPPLib/MFC/CString.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <PPPLib/Definitions.h>
#include <SpectralEvaluation/DateTime.h>

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
    /** The full path and file-name to the .pak-file containing the spectra
        from which this result was computed. */
    std::string m_pakFile = "";

    /** The full path and file-name of the evaluation log file that was generated.
        The file itself contains the result of the evaluation */
    std::string m_evalLogFile[MAX_FIT_WINDOWS];

    /** The full name of the fit-window that was used to generate each evaluation
        result, this is typically a name such as 'SO2', 'SO2_low' or 'O4'
         m_fitWindowName[i] is the name of the fit-window behind the result in m_evalLogFile[i]. */
    std::string m_fitWindowName[MAX_FIT_WINDOWS];

    /** The date and time that the scan was generated. In UTC, taken from the .pak-file  */
    novac::CDateTime m_startTime;

    /** The properties of this scan. This is only evaluated in the main-fit window */
    novac::CPlumeInScanProperty m_scanProperties;

};
}
