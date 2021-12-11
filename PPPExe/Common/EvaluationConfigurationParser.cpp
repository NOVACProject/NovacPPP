#include "EvaluationConfigurationParser.h"
#include "Common.h"
#include <memory>
#include <cstring>
#include <SpectralEvaluation/Configuration/DarkSettings.h>

using namespace FileHandler;
using namespace novac;


int CEvaluationConfigurationParser::ReadConfigurationFile(const novac::CString& fileName,
            Configuration::CEvaluationConfiguration& settings,
            Configuration::CDarkCorrectionConfiguration& darkSettings,
            Configuration::CInstrumentCalibrationConfiguration& calibrationSettings)
{

    // 1. Open the file
    if (!Open(fileName))
    {
        return FAIL;
    }

    // parse the file
    while (szToken = NextToken()) {
        // no use to parse empty lines
        if (strlen(szToken) < 3)
            continue;

        if (novac::Equals(szToken, "serial", strlen("serial"))) {
            this->Parse_StringItem("/serial", settings.m_serial);
            continue;
        }

        if (novac::Equals(szToken, "fitWindow", strlen("fitWindow"))) {
            novac::CFitWindow tmpWindow;
            novac::CDateTime validFrom, validTo;

            Parse_FitWindow(tmpWindow, validFrom, validTo);

            settings.InsertFitWindow(tmpWindow, &validFrom, &validTo);
        }

        if (novac::Equals(szToken, "DarkCorrection", strlen("DarkCorrection"))) {
            Configuration::CDarkSettings dSettings;
            novac::CDateTime validFrom, validTo;

            Parse_DarkCorrection(dSettings, validFrom, validTo);

            darkSettings.InsertDarkCurrentCorrectionSettings(dSettings, &validFrom, &validTo);
        }

        if (novac::Equals(szToken, "Calibration", strlen("Calibration"))) {
            Parse_CalibrationSettings(calibrationSettings);
        }
    }
    Close();

    return 0;
}

int CEvaluationConfigurationParser::WriteConfigurationFile(
    const novac::CString& fileName,
    const Configuration::CEvaluationConfiguration& settings,
    const Configuration::CDarkCorrectionConfiguration& darkSettings,
    const Configuration::CInstrumentCalibrationConfiguration& calibrationSettings)
{
    novac::CString indent, str;
    novac::CFitWindow window;
    Configuration::CDarkSettings dSettings;
    novac::CDateTime from, to;

    // open the file
    FILE* f = fopen(fileName, "w");
    if (f == NULL)
        return 1;

    // write the header lines and the start of the file
    fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(f, "<!-- This is the configuration file for the evaluation of spectra in the NOVAC Post Processing Program -->\n\n");
    fprintf(f, "<EvaluationConfiguration>\n");
    indent.Format("\t");

    // Write the serial-number of the spectrometer for which this configuration is valid
    fprintf(f, "\t<serial>%s</serial>\n", (const char*)settings.m_serial);

    // ------ loop through each of the fit windows and write them to file --------
    unsigned long nWindows = settings.GetFitWindowNum();
    for (unsigned int k = 0; k < nWindows; ++k) {
        settings.GetFitWindow(k, window, from, to);

        fprintf(f, "\t<fitWindow>\n");

        // the channel of the spectrometer
        fprintf(f, "\t\t<channel>%d</channel>\n", window.channel);

        // the name of the fit-window
        fprintf(f, "\t\t<name>%s</name>\n", window.name.c_str());

        // the time-range when the fit-window is valid
        fprintf(f, "\t\t<validFrom>%04d.%02d.%02d</validFrom>\n", from.year, from.month, from.day);
        fprintf(f, "\t\t<validTo>%04d.%02d.%02d</validTo>\n", to.year, to.month, to.day);

        // The size of the spectra and the interlace-steps
        fprintf(f, "\t\t<specLength>%d</specLength>\n", window.specLength);
        fprintf(f, "\t\t<interlaceStep>%d</interlaceStep>\n", window.interlaceStep);

        // the option for the polynomial to use
        fprintf(f, "\t\t<polyOrder>%d</polyOrder>\n", window.polyOrder);

        // the type of fit to use
        fprintf(f, "\t\t<fitType>%d</fitType>\n", window.fitType);

        // the boundaries of the fit (in pixels)
        fprintf(f, "\t\t<fitLow>%d</fitLow>\n", window.fitLow);
        fprintf(f, "\t\t<fitHigh>%d</fitHigh>\n", window.fitHigh);

        // If we should use a pre-calibrated solar-spectrum to calibrate
        //  the shift & squeeze of the spectra
        if (window.fraunhoferRef.m_path.size() > 3) {
            fprintf(f, "\t\t<wavelengthCalibration>\n");
            fprintf(f, "\t\t<fraunhoferSpec>%s</fraunhoferSpec>\n", window.fraunhoferRef.m_path.c_str());
            fprintf(f, "\t\t</wavelengthCalibration>\n");
        }

        // Each of the references...
        for (int j = 0; j < window.nRef; ++j) {
            fprintf(f, "\t\t<Reference>\n");
            fprintf(f, "\t\t\t<name>%s</name>\n", window.ref[j].m_specieName.c_str());
            fprintf(f, "\t\t\t<path>%s</path>\n", window.ref[j].m_path.c_str());

            // The value for the shift
            fprintf(f, "\t\t\t<shiftOption>%d</shiftOption>\n", window.ref[j].m_shiftOption);
            if (window.ref[j].m_shiftOption != novac::SHIFT_TYPE::SHIFT_FREE)
                fprintf(f, "\t\t\t<shiftValue>%lf</shiftValue>\n", window.ref[j].m_shiftValue);

            // The value for the squeeze
            fprintf(f, "\t\t\t<squeezeOption>%d</squeezeOption>\n", window.ref[j].m_shiftOption);
            if (window.ref[j].m_shiftOption != novac::SHIFT_TYPE::SHIFT_FREE)
                fprintf(f, "\t\t\t<squeezeValue>%lf</squeezeValue>\n", window.ref[j].m_shiftValue);

            // The value for the column
            fprintf(f, "\t\t\t<columnOption>%d</columnOption>\n", window.ref[j].m_columnOption);
            if (window.ref[j].m_columnOption != novac::SHIFT_TYPE::SHIFT_FREE)
                fprintf(f, "\t\t\t<columnValue>%lf</columnValue>\n", window.ref[j].m_columnValue);

            fprintf(f, "\t\t</Reference>\n");
        }

        fprintf(f, "\t</fitWindow>\n");
    }

    // ------ loop through each of the dark-current settings and write them to file --------
    for (unsigned int k = 0; k < nWindows; ++k) {
        darkSettings.GetDarkSettings(k, dSettings, from, to);

        fprintf(f, "\t<DarkCorrection>\n");

        // the time-range when the dark-current settings is valid
        fprintf(f, "\t\t<validFrom>%04d.%02d.%02d</validFrom>\n", from.year, from.month, from.day);
        fprintf(f, "\t\t<validTo>%04d.%02d.%02d</validTo>\n", to.year, to.month, to.day);

        if (dSettings.m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MEASURED_IN_SCAN) {
            // only use a dark-spectrum with the same exp.-time
            fprintf(f, "\t\t<dark>SCAN</dark>\n");
        }
        else if (dSettings.m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MODEL_ALWAYS) {
            // always model the dark-spectrum
            fprintf(f, "\t\t<dark>MODEL</dark>\n");

            // dark-current
            if (dSettings.m_darkCurrentOption == Configuration::DARK_MODEL_OPTION::MEASURED_IN_SCAN) {
                fprintf(f, "\t\t<darkCurrent>SCAN</darkCurrent>\n");
            }
            else {
                fprintf(f, "\t\t<darkCurrent>USER</darkCurrent>\n");
                fprintf(f, "\t\t<darkCurrentSpec>%s</darkCurrentSpec>\n", dSettings.m_darkCurrentSpec.c_str());
            }

            // offset
            if (dSettings.m_offsetOption == Configuration::DARK_MODEL_OPTION::MEASURED_IN_SCAN) {
                fprintf(f, "\t\t<offset>SCAN</offset>\n");
            }
            else {
                fprintf(f, "\t\t<offset>USER</offset>\n");
                fprintf(f, "\t\t<offsetSpec>%s</offsetSpec>\n", dSettings.m_offsetSpec.c_str());
            }

        }
        else if (dSettings.m_darkSpecOption == Configuration::DARK_SPEC_OPTION::USER_SUPPLIED) {
            fprintf(f, "\t\t<dark>USER</dark>\n");
            fprintf(f, "\t\t<darkCurrentSpec>%s</darkCurrentSpec>\n", dSettings.m_darkCurrentSpec.c_str());
            fprintf(f, "\t\t<offsetSpec>%s</offsetSpec>\n", dSettings.m_offsetSpec.c_str());
        }

        fprintf(f, "\t</DarkCorrection>\n");
    }

    // The instrument calibration settings
    if (calibrationSettings.m_initialCalibrationFile.size() > 0 || calibrationSettings.m_instrumentLineshapeFile.size() > 0)
    {
        fprintf(f, "\t<Calibration>\n");
        fprintf(f, "\t\t<initialCalibrationFile>%s</initialCalibrationFile>\n", calibrationSettings.m_initialCalibrationFile.c_str());
        fprintf(f, "\t\t<initialInstrumentLineshapeFile>%s</initialInstrumentLineshapeFile>\n", calibrationSettings.m_instrumentLineshapeFile.c_str());
        fprintf(f, "\t</Calibration>\n");
    }

    fprintf(f, "</EvaluationConfiguration>\n");

    // remember to close the file when we're done
    fclose(f);

    return 0;
}

void SaveSlitFunctionAndWavelengthCalibration(novac::CFitWindow& window, novac::CString& slitfunctionFile, novac::CString& wavelengthCalibFile)
{
    if (slitfunctionFile.GetLength() > 0)
    {
        for (int ii = 0; ii < window.nRef; ++ii)
        {
            window.ref[ii].m_slitFunctionFile = slitfunctionFile.ToStdString();
        }
    }
    if (wavelengthCalibFile.GetLength() > 0)
    {
        for (int ii = 0; ii < window.nRef; ++ii)
        {
            window.ref[ii].m_wavelengthCalibrationFile = wavelengthCalibFile.ToStdString();
        }
    }
}

int CEvaluationConfigurationParser::Parse_FitWindow(novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo) {
    window.Clear();
    novac::CString slitfunctionFile, wavelengthCalibFile;

    // parse the file
    while (szToken = NextToken()) {
        // no use to parse empty lines
        if (strlen(szToken) < 2)
            continue;

        // ignore comments
        if (Equals(szToken, "!--", 3)) {
            continue;
        }

        // end of fit-window section
        if (Equals(szToken, "/fitWindow")) {
            SaveSlitFunctionAndWavelengthCalibration(window, slitfunctionFile, wavelengthCalibFile);
            return 0;
        }

        if (Equals(szToken, "fitWindow"))
        {
            novac::CFitWindow child;
            novac::CDateTime childValidFrom, childValidTo;
            Parse_FitWindow(child, childValidFrom, childValidTo);
            window.child.push_back(child);
            continue;
        }


        if (Equals(szToken, "name")) {
            Parse_StringItem("/name", window.name);
            continue;
        }

        if (Equals(szToken, "validFrom")) {
            Parse_Date("/validFrom", validFrom);
            continue;
        }

        if (Equals(szToken, "validTo")) {
            Parse_Date("/validTo", validTo);
            continue;
        }

        if (Equals(szToken, "fitLow")) {
            Parse_IntItem("/fitLow", window.fitLow);
            continue;
        }

        if (Equals(szToken, "fitHigh")) {
            Parse_IntItem("/fitHigh", window.fitHigh);
            continue;
        }

        if (Equals(szToken, "polyOrder")) {
            Parse_IntItem("/polyOrder", window.polyOrder);
            continue;
        }

        if (Equals(szToken, "fitType")) {
            Parse_IntItem("/fitType", (int&)window.fitType); // TODO: Will this be ok????
            continue;
        }

        if (Equals(szToken, "channel")) {
            Parse_IntItem("/channel", window.channel);
            continue;
        }

        if (Equals(szToken, "specLength")) {
            Parse_IntItem("/specLength", window.specLength);
            continue;
        }

        if (Equals(szToken, "fOptShift")) {
            int flagToParse = 0;
            Parse_IntItem("/fOptShift", flagToParse);
            window.findOptimalShift = (flagToParse > 0);
            continue;
        }

        if (Equals(szToken, "shiftSky")) {
            int flagToParse = 0;
            Parse_IntItem("/shiftSky", flagToParse);
            window.shiftSky = (flagToParse > 0);
            continue;
        }

        if (Equals(szToken, "interlaceStep")) {
            Parse_IntItem("/interlaceStep", window.interlaceStep);
            continue;
        }

        if (Equals(szToken, "interlaced")) {
            Parse_IntItem("/interlaced", window.interlaceStep);
            window.interlaceStep += 1;
            continue;
        }

        if (Equals(szToken, "fraunhoferSpec", 14)) {
            // Parse the settings for the wavelength calibration
            this->Parse_PathItem("/fraunhoferSpec", window.fraunhoferRef.m_path);
            continue;
        }

        if (Equals(szToken, "slitFunction")) {
            // This is the path to a reference which needs to be convolved before we can continue.
            Parse_PathItem("/slitFunction", slitfunctionFile);
            continue;
        }

        if (Equals(szToken, "wavlengthCalibration")) {
            // This is the path to a reference which needs to be convolved before we can continue.
            Parse_PathItem("/wavlengthCalibration", wavelengthCalibFile);
            continue;
        }

        if (Equals(szToken, "Reference", 9)) {
            Parse_Reference(window);
            continue;
        }
    }

    return 1;
}

int CEvaluationConfigurationParser::Parse_CalibrationSettings(Configuration::CInstrumentCalibrationConfiguration& calibrationSettings)
{
    while (szToken = NextToken()) {
        // no use to parse empty lines
        if (strlen(szToken) < 2)
            continue;

        // ignore comments
        if (Equals(szToken, "!--", 3)) {
            continue;
        }

        // end of dark-correction section
        if (Equals(szToken, "/Calibration")) {
            return 0;
        }

        if (Equals(szToken, "initialCalibrationFile")) {
            Parse_StringItem("/initialCalibrationFile", calibrationSettings.m_initialCalibrationFile);
            continue;
        }

        if (Equals(szToken, "initialInstrumentLineshapeFile")) {
            Parse_StringItem("/initialInstrumentLineshapeFile", calibrationSettings.m_instrumentLineshapeFile);
            continue;
        }
    }

    return 1;
}


int CEvaluationConfigurationParser::Parse_Reference(novac::CFitWindow& window) {
    int nRef = window.nRef;

    // the actual reading loop
    while (szToken = NextToken()) {

        // no use to parse empty lines
        if (strlen(szToken) < 3)
            continue;

        // ignore comments
        if (Equals(szToken, "!--", 3)) {
            continue;
        }

        if (Equals(szToken, "/Reference")) {
            ++window.nRef;
            return 0;
        }

        if (Equals(szToken, "name")) {
            Parse_StringItem("/name", window.ref[nRef].m_specieName);
            continue;
        }

        if (Equals(szToken, "filtered")) {
            novac::CString str;
            Parse_StringItem("/filtered", str);
            if (Equals(str, "HP500")) {
                window.ref[nRef].m_isFiltered = true;
            }
            continue;
        }

        if (Equals(szToken, "path")) {
            // This is the path to a pre-convolved reference. Just read the path and read the reference from there.
            Parse_PathItem("/path", window.ref[nRef].m_path);
            continue;
        }

        if (Equals(szToken, "crossSection")) {
            // This is the path to a reference which needs to be convolved before we can continue.
            Parse_PathItem("/crossSection", window.ref[nRef].m_crossSectionFile);
            continue;
        }

        if (Equals(szToken, "shiftOption")) {
            int tmpInt;
            Parse_IntItem("/shiftOption", tmpInt);
            switch (tmpInt) {
            case 0: window.ref[nRef].m_shiftOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: window.ref[nRef].m_shiftOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: window.ref[nRef].m_shiftOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: window.ref[nRef].m_shiftOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            }
            continue;
        }

        if (Equals(szToken, "shiftValue")) {
            Parse_FloatItem("/shiftValue", window.ref[nRef].m_shiftValue);
            continue;
        }

        if (Equals(szToken, "squeezeOption")) {
            int tmpInt;
            Parse_IntItem("/squeezeOption", tmpInt);
            switch (tmpInt) {
            case 0: window.ref[nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: window.ref[nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: window.ref[nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: window.ref[nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            }
            continue;
        }

        if (Equals(szToken, "squeezeValue")) {
            Parse_FloatItem("/squeezeValue", window.ref[nRef].m_squeezeValue);
            continue;
        }

        if (Equals(szToken, "columnOption")) {
            int tmpInt;
            Parse_IntItem("/columnOption", tmpInt);
            switch (tmpInt) {
            case 0: window.ref[nRef].m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: window.ref[nRef].m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: window.ref[nRef].m_columnOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: window.ref[nRef].m_columnOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            }
            continue;
        }

        if (Equals(szToken, "columnValue")) {
            Parse_FloatItem("/columnValue", window.ref[nRef].m_columnValue);
            continue;
        }
    }

    return FAIL;
}

int CEvaluationConfigurationParser::Parse_DarkCorrection(Configuration::CDarkSettings& dSettings, novac::CDateTime& validFrom, novac::CDateTime& validTo) {
    dSettings.Clear();
    novac::CString str;

    while (szToken = NextToken()) {
        // no use to parse empty lines
        if (strlen(szToken) < 2)
            continue;

        // ignore comments
        if (Equals(szToken, "!--", 3)) {
            continue;
        }

        // end of dark-correction section
        if (Equals(szToken, "/DarkCorrection")) {
            return 0;
        }

        // valid interval
        if (Equals(szToken, "validFrom")) {
            Parse_Date("/validFrom", validFrom);
            continue;
        }
        if (Equals(szToken, "validTo")) {
            Parse_Date("/validTo", validTo);
            continue;
        }

        // the option for the dark
        if (Equals(szToken, "dark")) {
            Parse_StringItem("/dark", str);

            if (Equals(str, "MODEL")) {
                dSettings.m_darkSpecOption = Configuration::DARK_SPEC_OPTION::MODEL_ALWAYS;
            }
            else if (Equals(str, "USER")) {
                dSettings.m_darkSpecOption = Configuration::DARK_SPEC_OPTION::USER_SUPPLIED;
            }
            else {
                dSettings.m_darkSpecOption = Configuration::DARK_SPEC_OPTION::MEASURED_IN_SCAN;
            }
            continue;
        }

        if (Equals(szToken, "darkCurrentSpec")) {
            Parse_StringItem("/darkCurrentSpec", dSettings.m_darkCurrentSpec);
            continue;
        }

        if (Equals(szToken, "darkCurrent")) {
            Parse_StringItem("/darkCurrent", str);

            if (Equals(str, "USER")) {
                dSettings.m_darkCurrentOption = Configuration::DARK_MODEL_OPTION::USER_SUPPLIED;
            }
            else {
                dSettings.m_darkCurrentOption = Configuration::DARK_MODEL_OPTION::MEASURED_IN_SCAN;
            }
            continue;
        }

        if (Equals(szToken, "offsetSpec")) {
            Parse_StringItem("/offsetSpec", dSettings.m_offsetSpec);
            continue;
        }

        if (Equals(szToken, "offset")) {
            Parse_StringItem("/offset", str);

            if (Equals(str, "USER")) {
                dSettings.m_offsetOption = Configuration::DARK_MODEL_OPTION::USER_SUPPLIED;
            }
            else {
                dSettings.m_offsetOption = Configuration::DARK_MODEL_OPTION::MEASURED_IN_SCAN;
            }
            continue;
        }
    }

    return 1;
}

