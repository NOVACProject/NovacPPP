#ifndef NOVAC_PPPLIB_CFILE_UTILS_H
#define NOVAC_PPPLIB_CFILE_UTILS_H

#include <PPPLib/MFC/CString.h>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/NovacEnums.h>

namespace novac
{
class CFileUtils
{
public:
    /** Take out the file name from a long path
        @param fileName path of the file	*/
    static void GetFileName(novac::CString& fileName);

    /** Take out the directory from a long path name.
        @param fileName - the complete path of the file */
    static void GetDirectory(novac::CString& fileName);

    /** Takes the filename of an evaluation log or pak-file and extracts the
        Serial-number of the spectrometer, the date the scan was performed
        and the start-time of the scan from the filename.
        NOTICE: This only looks at the file name and works if the file name is on the standard of the NovacProgram
        which is "D2J2134_170129_0317_1.pak".
        @return true If the parsing is successful. */
    static bool GetInfoFromFileName(const novac::CString fileName, CDateTime& start, novac::CString& serial, int& channel, MeasurementMode& mode);

    /** Judges if the provided .pak file is a complete file from the file name only.
        @return true if the file is from an incomplete scan */
    static bool IsIncompleteFile(const novac::CString& fileName);
};
}  // namespace novac

#endif  // NOVAC_PPPLIB_CFILE_UTILS_H