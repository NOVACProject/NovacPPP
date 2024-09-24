#pragma once

#include <PPPLib/PPPLib.h>
#include <PPPLib/Evaluation/ScanResult.h>
#include <SpectralEvaluation/Definitions.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/MFC/CArray.h>
#include <SpectralEvaluation/Definitions.h>

namespace FileHandler
{

class CEvaluationLogFileHandler
{
public:
    CEvaluationLogFileHandler();

    /** The evaluation log */
    novac::CString m_evaluationLog;

    // ------------------- PUBLIC METHODS -------------------------

    /** Reads the conents of the provided evaluation log and fills in all the members of this class. */
    RETURN_CODE ReadEvaluationLog();

    /** Writes the contents of the array 'm_scan' to a new evaluation-log file */
    // TODO:refactor the software version here, this isn't that pretty.
    RETURN_CODE WriteEvaluationLog(const novac::CString fileName, int softwareMajorNumber, int softwareMinorNumber);

    /** Returns true if the scan number 'scanNo' in the most recently read
            evaluation log file is a wind speed measurement. */
    bool IsWindSpeedMeasurement(int scanNo);

    /** Returns true if the scan number 'scanNo' in the most recently read
            evaluation log file is a wind speed measurement of heidelberg type. */
    bool IsWindSpeedMeasurement_Heidelberg(int scanNo);

    /** Appends the evaluation result of one spectrum to the given string.
            @param info - the information about the spectrum
            @param result - the evaluation result, can be NULL
            @param string - will on return be filled with the output line to be written to the evaluation-log.
            @return SUCCESS - always */
    static RETURN_CODE FormatEvaluationResult(const novac::CSpectrumInfo* info, const novac::CEvaluationResult* result, INSTRUMENT_TYPE iType, double maxIntensity, int nSpecies, novac::CString& string);

    // ------------------- PUBLIC DATA -------------------------

    /** Information from the evaluated scans */
    std::vector<Evaluation::CScanResult> m_scan;

    /** Information of the wind field used to calculate the flux of each scan */
    std::vector<Meteorology::CWindField> m_windField;

    /** The species that were found in this evaluation log */
    novac::CString m_specie[20];

    /** The number of species found in the evaluation log */
    long m_specieNum;

    /** The currently selected specie */
    long m_curSpecie;

    /** The instrument-type for the instrument that generated the results */
    INSTRUMENT_TYPE m_instrumentType;

    /** The additional spectrum information of one spectrum. */
    novac::CSpectrumInfo m_specInfo;

protected:

    typedef struct LogColumns
    {
        int column[MAX_N_REFERENCES];
        int columnError[MAX_N_REFERENCES];
        int shift[MAX_N_REFERENCES];
        int shiftError[MAX_N_REFERENCES];
        int squeeze[MAX_N_REFERENCES];
        int squeezeError[MAX_N_REFERENCES];
        int intensity;
        int fitIntensity;
        int peakSaturation;
        int fitSaturation;
        int offset;
        int delta;
        int chiSquare;
        int nSpec;
        int expTime;
        int position;
        int position2;
        int nSpecies;
        int starttime;
        int stoptime;
        int name;
    }LogColumns;

    /** Data structure to remember what column corresponds to which value in the evaluation log */
    LogColumns m_col;

    /** The result from the evaluation of one spectrum. */
    novac::CEvaluationResult m_evResult;

    /** Reads the header line for the scan information and retrieves which
        column represents which value. */
    void ParseScanHeader(const char szLine[8192]);

    /** Reads and parses the XML-shaped 'scanInfo' header before the scan */
    void ParseScanInformation(novac::CSpectrumInfo& scanInfo, double& flux, FILE* f);

    /** Reads and parses the XML-shaped 'fluxInfo' header before the scan */
    void ParseFluxInformation(Meteorology::CWindField& windField, double& flux, FILE* f);

    /** Resets the information about which column data is stored in */
    void ResetColumns();

    /** Resets the old scan information */
    void ResetScanInformation();

    /** Makes a quick scan through the evaluation-log
        to count the number of scans in it */
    long CountScansInFile();

    /** Sorts the scans in order of collection */
    void SortScans();

    /** Returns true if the scans are already ordered */
    bool IsSorted();

    /** Sorts the CScanResult-objects in the given array.
            Algorithm based on MergeSort (~O(NlogN)) */
    static void SortScans(novac::CArray<Evaluation::CScanResult, Evaluation::CScanResult&>& array, bool ascending = true);

    /** Sorts the CScanResult-objects in the given array.
            Algorithm based on BubbleSort (~O(N2))
            Quite efficient for small arrays since the elements does not have to be copied
                and thus uses very little memory */
    static void BubbleSortScans(novac::CArray<Evaluation::CScanResult, Evaluation::CScanResult&>& array, bool ascending = true);

    /** Merges the two arrays in a sorted way and stores the
            result in the output-array 'result' */
    static void MergeArrays(novac::CArray<Evaluation::CScanResult, Evaluation::CScanResult&>& array1, novac::CArray<Evaluation::CScanResult, Evaluation::CScanResult&>& array2, novac::CArray<Evaluation::CScanResult, Evaluation::CScanResult&>& result, bool ascending = true);
};

}