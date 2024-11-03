#pragma once

#include <SpectralEvaluation/Definitions.h>
#include <SpectralEvaluation/Log.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

#include <PPPLib/PPPLib.h>
#include <PPPLib/Evaluation/ScanResult.h>
#include <PPPLib/Configuration/UserConfiguration.h>

namespace FileHandler
{

class CEvaluationLogFileHandler
{
public:
    CEvaluationLogFileHandler(
        novac::ILogger& log,
        std::string evaluationLog,
        novac::Molecule molecule,
        novac::SpectrometerModel* spectrometerModel = nullptr);

    /** The evaluation log */
    const std::string m_evaluationLog;

    /** The model of the spectrometer which created the scans in this file.
        If this is specified in the constructor, then this will be used when reading the data.
        If this is NOT specified in the constructor then this will be guessed from the data when reading. */
    novac::SpectrometerModel m_spectrometerModel;

    // ------------------- PUBLIC METHODS -------------------------

    /** Reads the conents of the provided evaluation log and fills in all the members of this class. */
    RETURN_CODE ReadEvaluationLog();

    /** Writes the contents of the array 'm_scan' to a new evaluation-log file */
    // TODO:refactor the software version here, this isn't that pretty.
    RETURN_CODE WriteEvaluationLog(const std::string& fileName, novac::SpectrometerModel spectrometerModel, int softwareMajorNumber, int softwareMinorNumber);

    /** Appends the evaluation result of one spectrum to the given string.
            @param info - the information about the spectrum
            @param result - the evaluation result, can be NULL
            @param string - will on return be filled with the output line to be written to the evaluation-log.
            @return SUCCESS - always */
    static RETURN_CODE FormatEvaluationResult(const novac::CSpectrumInfo* info, const novac::CEvaluationResult* result, novac::NovacInstrumentType iType, double maxIntensity, size_t nSpecies, novac::CString& string);

    // ------------------- PUBLIC DATA -------------------------

    /** Information from the evaluated scans */
    std::vector<Evaluation::CScanResult> m_scan;

    /** Information of the wind field used to calculate the flux of each scan */
    std::vector<Meteorology::WindField> m_windField;

    /** The names of the species that were found in this evaluation log */
    std::vector<std::string> m_specieName;

    novac::Molecule m_molecule;

    /** The instrument-type for the instrument that generated the results */
    novac::NovacInstrumentType m_instrumentType;

    /** The additional spectrum information of one spectrum. */
    novac::CSpectrumInfo m_specInfo;

private:

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

    novac::ILogger& m_log;

    /** Reads the header line for the scan information and retrieves which
        column represents which value. */
    void ParseScanHeader(const char szLine[8192]);

    /** Reads and parses the XML-shaped 'scanInfo' header before the scan */
    void ParseScanInformation(novac::CSpectrumInfo& scanInfo, double& flux, FILE* f);

    /** Reads and parses the XML-shaped 'fluxInfo' header before the scan */
    void ParseFluxInformation(Meteorology::WindField& windField, double& flux, FILE* f);

    /** Resets the information about which column data is stored in */
    void ResetColumns();

    /** Resets the old scan information */
    void ResetScanInformation();

    /** Makes a quick scan through the evaluation-log
        to count the number of scans in it */
    size_t CountScansInFile();

    /** Sorts the scans in order of collection */
    void SortScans();

    /** Returns true if the scans are already ordered */
    bool IsSorted();
};

}