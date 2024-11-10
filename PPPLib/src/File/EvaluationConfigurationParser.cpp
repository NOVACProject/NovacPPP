#include <PPPLib/File/EvaluationConfigurationParser.h>
#include <PPPLib/Logging.h>
#include <PPPLib/PPPLib.h>

#include <memory>
#include <cstring>
#include <SpectralEvaluation/Configuration/DarkSettings.h>

using namespace novac;

namespace FileHandler
{

static char start = 's';

RETURN_CODE CEvaluationConfigurationParser::ReadConfigurationFile(const novac::CString& fileName,
    Configuration::CEvaluationConfiguration& settings,
    Configuration::CDarkCorrectionConfiguration& darkSettings,
    Configuration::CInstrumentCalibrationConfiguration& calibrationSettings)
{
    // 1. Open the file
    if (!Open(fileName))
    {
        return RETURN_CODE::FAIL;
    }

    // parse the file, one line at a time.
    szToken = &start;
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        if (novac::Equals(szToken, "serial", strlen("serial")))
        {
            this->Parse_StringItem("/serial", settings.m_serial);
            continue;
        }

        if (novac::Equals(szToken, "fitWindow", strlen("fitWindow")))
        {
            novac::CFitWindow tmpWindow;
            novac::CDateTime validFrom, validTo;

            if (RETURN_CODE::FAIL == Parse_FitWindow(tmpWindow, validFrom, validTo))
            {
                return RETURN_CODE::FAIL;
            }

            settings.InsertFitWindow(tmpWindow, validFrom, validTo);
        }

        if (novac::Equals(szToken, "DarkCorrection", strlen("DarkCorrection")))
        {
            Configuration::CDarkSettings dSettings;
            novac::CDateTime validFrom, validTo;

            if (RETURN_CODE::FAIL == Parse_DarkCorrection(dSettings, validFrom, validTo))
            {
                return RETURN_CODE::FAIL;
            }

            darkSettings.InsertDarkCurrentCorrectionSettings(dSettings, validFrom, validTo);
        }

        if (novac::Equals(szToken, "Calibration", strlen("Calibration")))
        {
            if (RETURN_CODE::FAIL == Parse_CalibrationSettings(calibrationSettings))
            {
                return RETURN_CODE::FAIL;
            }
        }
    }
    Close();

    return RETURN_CODE::SUCCESS;
}

RETURN_CODE CEvaluationConfigurationParser::WriteConfigurationFile(
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
        return RETURN_CODE::FAIL;

    // write the header lines and the start of the file
    fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(f, "<!-- This is the configuration file for the evaluation of spectra in the NOVAC Post Processing Program -->\n\n");
    fprintf(f, "<EvaluationConfiguration>\n");
    indent.Format("\t");

    // Write the serial-number of the spectrometer for which this configuration is valid
    fprintf(f, "\t<serial>%s</serial>\n", settings.m_serial.c_str());

    // ------ loop through each of the fit windows and write them to file --------
    const size_t nWindows = settings.NumberOfFitWindows();
    for (size_t k = 0; k < nWindows; ++k)
    {
        settings.GetFitWindow(k, window, from, to);

        fprintf(f, "\t<fitWindow>\n");

        // the channel of the spectrometer
        fprintf(f, "\t\t<channel>%d</channel>\n", window.channel);

        // the name of the fit-window
        fprintf(f, "\t\t<name>%s</name>\n", window.name.c_str());

        // the time-range when the fit-window is valid
        fprintf(f, "\t\t<validFrom>%04d.%02d.%02dT%02d:%02d:%02d</validFrom>\n", from.year, from.month, from.day, from.hour, from.minute, from.second);
        fprintf(f, "\t\t<validTo>%04d.%02d.%02dT%02d:%02d:%02d</validTo>\n", to.year, to.month, to.day, to.hour, to.minute, to.second);

        // The size of the spectra and the interlace-steps
        fprintf(f, "\t\t<specLength>%d</specLength>\n", window.specLength);
        fprintf(f, "\t\t<interlaceStep>%d</interlaceStep>\n", window.interlaceStep);

        // the option for the polynomial to use
        fprintf(f, "\t\t<polyOrder>%d</polyOrder>\n", window.polyOrder);

        // the type of fit to use
        fprintf(f, "\t\t<fitType>%d</fitType>\n", (int)window.fitType);

        // the boundaries of the fit (in pixels)
        fprintf(f, "\t\t<fitLow>%d</fitLow>\n", window.fitLow);
        fprintf(f, "\t\t<fitHigh>%d</fitHigh>\n", window.fitHigh);

        // If we should use a pre-calibrated solar-spectrum to calibrate
        //  the shift & squeeze of the spectra
        if (window.fraunhoferRef.m_path.size() > 3)
        {
            fprintf(f, "\t\t<wavelengthCalibration>\n");
            fprintf(f, "\t\t\t<fraunhoferSpec>%s</fraunhoferSpec>\n", window.fraunhoferRef.m_path.c_str());
            fprintf(f, "\t\t</wavelengthCalibration>\n");
        }

        // Each of the references...
        for (const auto& reference : window.reference)
        {
            fprintf(f, "\t\t<Reference>\n");
            fprintf(f, "\t\t\t<name>%s</name>\n", reference.m_specieName.c_str());
            fprintf(f, "\t\t\t<path>%s</path>\n", reference.m_path.c_str());

            // The value for the shift
            fprintf(f, "\t\t\t<shiftOption>%d</shiftOption>\n", (int)reference.m_shiftOption);
            if (reference.m_shiftOption != novac::SHIFT_TYPE::SHIFT_FREE)
            {
                fprintf(f, "\t\t\t<shiftValue>%lf</shiftValue>\n", reference.m_shiftValue);
            }

            // The value for the squeeze
            fprintf(f, "\t\t\t<squeezeOption>%d</squeezeOption>\n", (int)reference.m_squeezeOption);
            if (reference.m_squeezeOption != novac::SHIFT_TYPE::SHIFT_FREE)
            {
                fprintf(f, "\t\t\t<squeezeValue>%lf</squeezeValue>\n", reference.m_squeezeValue);
            }

            // The value for the column
            fprintf(f, "\t\t\t<columnOption>%d</columnOption>\n", (int)reference.m_columnOption);
            if (reference.m_columnOption != novac::SHIFT_TYPE::SHIFT_FREE)
            {
                fprintf(f, "\t\t\t<columnValue>%lf</columnValue>\n", reference.m_columnValue);
            }
            fprintf(f, "\t\t</Reference>\n");
        }

        fprintf(f, "\t</fitWindow>\n");
    }

    // ------ loop through each of the dark-current settings and write them to file --------
    const size_t numberOfDarkSettings = darkSettings.GetSettingsNum();
    for (size_t k = 0; k < numberOfDarkSettings; ++k)
    {
        darkSettings.GetDarkSettings(k, dSettings, from, to);

        fprintf(f, "\t<DarkCorrection>\n");

        // the time-range when the dark-current settings is valid
        fprintf(f, "\t\t<validFrom>%04d.%02d.%02dT%02d:%02d:%02d</validFrom>\n", from.year, from.month, from.day, from.hour, from.minute, from.second);
        fprintf(f, "\t\t<validTo>%04d.%02d.%02dT%02d:%02d:%02d</validTo>\n", to.year, to.month, to.day, to.hour, to.minute, to.second);

        if (dSettings.m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MEASURED_IN_SCAN)
        {
            // only use a dark-spectrum with the same exp.-time
            fprintf(f, "\t\t<dark>SCAN</dark>\n");
        }
        else if (dSettings.m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MODEL_ALWAYS)
        {
            // always model the dark-spectrum
            fprintf(f, "\t\t<dark>MODEL</dark>\n");

            // dark-current
            if (dSettings.m_darkCurrentOption == Configuration::DARK_MODEL_OPTION::MEASURED_IN_SCAN)
            {
                fprintf(f, "\t\t<darkCurrent>SCAN</darkCurrent>\n");
            }
            else
            {
                fprintf(f, "\t\t<darkCurrent>USER</darkCurrent>\n");
                fprintf(f, "\t\t<darkCurrentSpec>%s</darkCurrentSpec>\n", dSettings.m_darkCurrentSpec.c_str());
            }

            // offset
            if (dSettings.m_offsetOption == Configuration::DARK_MODEL_OPTION::MEASURED_IN_SCAN)
            {
                fprintf(f, "\t\t<offset>SCAN</offset>\n");
            }
            else
            {
                fprintf(f, "\t\t<offset>USER</offset>\n");
                fprintf(f, "\t\t<offsetSpec>%s</offsetSpec>\n", dSettings.m_offsetSpec.c_str());
            }

        }
        else if (dSettings.m_darkSpecOption == Configuration::DARK_SPEC_OPTION::USER_SUPPLIED)
        {
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
        if (calibrationSettings.m_instrumentLineshapeFile.size() > 0)
        {
            fprintf(f, "\t\t<initialInstrumentLineshapeFile>%s</initialInstrumentLineshapeFile>\n", calibrationSettings.m_instrumentLineshapeFile.c_str());
        }
        fprintf(f, "\t</Calibration>\n");
    }

    fprintf(f, "</EvaluationConfiguration>\n");

    // remember to close the file when we're done
    fclose(f);

    return RETURN_CODE::SUCCESS;
}

static void SaveSlitFunctionAndWavelengthCalibration(novac::CFitWindow& window, const std::string& slitfunctionFile, const std::string& wavelengthCalibFile)
{
    if (slitfunctionFile.size() > 0)
    {
        for (auto& reference : window.reference)
        {
            reference.m_slitFunctionFile = slitfunctionFile;
        }
    }
    if (wavelengthCalibFile.size() > 0)
    {
        for (auto& reference : window.reference)
        {
            reference.m_wavelengthCalibrationFile = wavelengthCalibFile;
        }
    }
}

RETURN_CODE CEvaluationConfigurationParser::Parse_FitWindow(novac::CFitWindow& window, novac::CDateTime& validFrom, novac::CDateTime& validTo)
{
    window.Clear();
    std::string slitfunctionFile, wavelengthCalibFile;

    // parse the file, one line at a time.
    szToken = &start;
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        // ignore comments
        if (Equals(szToken, "!--", 3))
        {
            continue;
        }

        // end of fit-window section
        if (Equals(szToken, "/fitWindow"))
        {
            SaveSlitFunctionAndWavelengthCalibration(window, slitfunctionFile, wavelengthCalibFile);
            return RETURN_CODE::SUCCESS;
        }

        if (Equals(szToken, "fitWindow"))
        {
            novac::CFitWindow child;
            novac::CDateTime childValidFrom, childValidTo;
            Parse_FitWindow(child, childValidFrom, childValidTo);
            window.child.push_back(child);
            continue;
        }

        if (Equals(szToken, "offsetFrom"))
        {
            Parse_SizeItem("/offsetFrom", window.offsetRemovalRange.from);
            continue;
        }

        if (Equals(szToken, "offsetTo"))
        {
            Parse_SizeItem("/offsetTo", window.offsetRemovalRange.to);
            continue;
        }

        if (Equals(szToken, "name"))
        {
            Parse_StringItem("/name", window.name);
            continue;
        }

        if (Equals(szToken, "validFrom"))
        {
            Parse_Date("/validFrom", validFrom);
            continue;
        }

        if (Equals(szToken, "validTo"))
        {
            Parse_Date("/validTo", validTo);
            continue;
        }

        if (Equals(szToken, "fitLow"))
        {
            Parse_IntItem("/fitLow", window.fitLow);
            continue;
        }

        if (Equals(szToken, "fitHigh"))
        {
            Parse_IntItem("/fitHigh", window.fitHigh);
            continue;
        }

        if (Equals(szToken, "polyOrder"))
        {
            Parse_IntItem("/polyOrder", window.polyOrder);
            continue;
        }

        if (Equals(szToken, "fitType"))
        {
            int fitType = 0;
            std::string valueParsed;
            Parse_IntItem("/fitType", fitType, valueParsed);
            if (fitType < 0 || fitType >(int)FIT_TYPE::FIT_POLY)
            {
                throw EvaluationConfigurationException(m_filename, "FitType does not match any known fit type.", valueParsed);
            }
            window.fitType = (FIT_TYPE)fitType;
            continue;
        }

        if (Equals(szToken, "channel"))
        {
            Parse_IntItem("/channel", window.channel);
            continue;
        }

        if (Equals(szToken, "specLength"))
        {
            Parse_IntItem("/specLength", window.specLength);
            continue;
        }

        if (Equals(szToken, "fOptShift"))
        {
            int flagToParse = 0;
            Parse_IntItem("/fOptShift", flagToParse);
            window.findOptimalShift = (flagToParse > 0);
            continue;
        }

        if (Equals(szToken, "shiftSky"))
        {
            int flagToParse = 0;
            Parse_IntItem("/shiftSky", flagToParse);
            window.shiftSky = (flagToParse > 0);
            continue;
        }

        if (Equals(szToken, "interlaceStep"))
        {
            Parse_IntItem("/interlaceStep", window.interlaceStep);
            continue;
        }

        if (Equals(szToken, "interlaced"))
        {
            Parse_IntItem("/interlaced", window.interlaceStep);
            window.interlaceStep += 1;
            continue;
        }

        if (Equals(szToken, "fraunhoferSpec", 14))
        {
            // Parse the settings for the wavelength calibration
            this->Parse_PathItem("/fraunhoferSpec", window.fraunhoferRef.m_path);
            continue;
        }

        if (Equals(szToken, "slitFunction"))
        {
            // This is the path to a reference which needs to be convolved before we can continue.
            Parse_PathItem("/slitFunction", slitfunctionFile);
            continue;
        }

        if (Equals(szToken, "wavlengthCalibration"))
        {
            // This is the path to a reference which needs to be convolved before we can continue.
            Parse_PathItem("/wavlengthCalibration", wavelengthCalibFile);
            continue;
        }

        if (Equals(szToken, "Reference", 9))
        {
            if (RETURN_CODE::FAIL == Parse_Reference(window))
            {
                return RETURN_CODE::FAIL;
            }
            continue;
        }
    }

    return RETURN_CODE::FAIL;
}

RETURN_CODE CEvaluationConfigurationParser::Parse_CalibrationSettings(Configuration::CInstrumentCalibrationConfiguration& calibrationSettings)
{
    // parse the file, one line at a time.
    szToken = &start;
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        // ignore comments
        if (Equals(szToken, "!--", 3))
        {
            continue;
        }

        // end of dark-correction section
        if (Equals(szToken, "/Calibration"))
        {
            return RETURN_CODE::SUCCESS;
        }

        if (Equals(szToken, "initialCalibrationFile"))
        {
            Parse_PathItem("/initialCalibrationFile", calibrationSettings.m_initialCalibrationFile);
            continue;
        }

        if (Equals(szToken, "initialInstrumentLineshapeFile"))
        {
            Parse_PathItem("/initialInstrumentLineshapeFile", calibrationSettings.m_instrumentLineshapeFile);
            continue;
        }
    }

    return RETURN_CODE::FAIL;
}


RETURN_CODE CEvaluationConfigurationParser::Parse_Reference(novac::CFitWindow& window)
{
    novac::CReferenceFile newReference;

    // parse the file, one line at a time.
    szToken = &start;
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        // ignore comments
        if (Equals(szToken, "!--", 3))
        {
            continue;
        }

        if (Equals(szToken, "/Reference"))
        {
            window.reference.push_back(newReference);
            return RETURN_CODE::SUCCESS;
        }

        if (Equals(szToken, "name"))
        {
            Parse_StringItem("/name", newReference.m_specieName);
            continue;
        }

        if (Equals(szToken, "filtered"))
        {
            novac::CString str;
            Parse_StringItem("/filtered", str);
            if (Equals(str, "HP500"))
            {
                newReference.m_isFiltered = true;
            }
            continue;
        }

        if (Equals(szToken, "path"))
        {
            // This is the path to a pre-convolved reference. Just read the path and read the reference from there.
            Parse_PathItem("/path", newReference.m_path);
            continue;
        }

        if (Equals(szToken, "crossSection"))
        {
            // This is the path to a reference which needs to be convolved before we can continue.
            Parse_PathItem("/crossSection", newReference.m_crossSectionFile);
            continue;
        }

        if (Equals(szToken, "shiftOption"))
        {
            int tmpInt;
            std::string valueParsed;
            Parse_IntItem("/shiftOption", tmpInt, valueParsed);
            switch (tmpInt)
            {
            case 0: newReference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: newReference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: newReference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: newReference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            default: throw EvaluationConfigurationException(m_filename, "Failed to parse shift option type.", valueParsed);
            }
            continue;
        }

        if (Equals(szToken, "shiftValue"))
        {
            Parse_FloatItem("/shiftValue", newReference.m_shiftValue);
            continue;
        }

        if (Equals(szToken, "squeezeOption"))
        {
            int tmpInt;
            std::string valueParsed;
            Parse_IntItem("/squeezeOption", tmpInt, valueParsed);
            switch (tmpInt)
            {
            case 0: newReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: newReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: newReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: newReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            default: throw EvaluationConfigurationException(m_filename, "Failed to parse squeeze option type.", valueParsed);
            }
            continue;
        }

        if (Equals(szToken, "squeezeValue"))
        {
            Parse_FloatItem("/squeezeValue", newReference.m_squeezeValue);
            continue;
        }

        if (Equals(szToken, "columnOption"))
        {
            int tmpInt;
            std::string valueParsed;
            Parse_IntItem("/columnOption", tmpInt, valueParsed);
            switch (tmpInt)
            {
            case 0: newReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE; break;
            case 1: newReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX; break;
            case 2: newReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_LINK; break;
            case 3: newReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_LIMIT; break;
            default: throw EvaluationConfigurationException(m_filename, "Failed to parse column option type.", valueParsed);
            }
            continue;
        }

        if (Equals(szToken, "columnValue"))
        {
            Parse_FloatItem("/columnValue", newReference.m_columnValue);
            continue;
        }
    }

    return RETURN_CODE::FAIL;
}

RETURN_CODE CEvaluationConfigurationParser::Parse_DarkCorrection(Configuration::CDarkSettings& dSettings, novac::CDateTime& validFrom, novac::CDateTime& validTo)
{
    dSettings.Clear();
    novac::CString str;

    // parse the file, one line at a time.
    szToken = &start;
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        // ignore comments
        if (Equals(szToken, "!--", 3))
        {
            continue;
        }

        // end of dark-correction section
        if (Equals(szToken, "/DarkCorrection"))
        {
            return RETURN_CODE::SUCCESS;
        }

        // valid interval
        if (Equals(szToken, "validFrom"))
        {
            Parse_Date("/validFrom", validFrom);
            continue;
        }
        if (Equals(szToken, "validTo"))
        {
            Parse_Date("/validTo", validTo);
            continue;
        }

        // the option for the dark
        if (Equals(szToken, "dark"))
        {
            Parse_StringItem("/dark", str);

            if (Equals(str, "MODEL"))
            {
                dSettings.m_darkSpecOption = Configuration::DARK_SPEC_OPTION::MODEL_ALWAYS;
            }
            else if (Equals(str, "USER"))
            {
                dSettings.m_darkSpecOption = Configuration::DARK_SPEC_OPTION::USER_SUPPLIED;
            }
            else
            {
                dSettings.m_darkSpecOption = Configuration::DARK_SPEC_OPTION::MEASURED_IN_SCAN;
            }
            continue;
        }

        if (Equals(szToken, "darkCurrentSpec"))
        {
            Parse_StringItem("/darkCurrentSpec", dSettings.m_darkCurrentSpec);
            continue;
        }

        if (Equals(szToken, "darkCurrent"))
        {
            Parse_StringItem("/darkCurrent", str);

            if (Equals(str, "USER"))
            {
                dSettings.m_darkCurrentOption = Configuration::DARK_MODEL_OPTION::USER_SUPPLIED;
            }
            else
            {
                dSettings.m_darkCurrentOption = Configuration::DARK_MODEL_OPTION::MEASURED_IN_SCAN;
            }
            continue;
        }

        if (Equals(szToken, "offsetSpec"))
        {
            Parse_StringItem("/offsetSpec", dSettings.m_offsetSpec);
            continue;
        }

        if (Equals(szToken, "offset"))
        {
            Parse_StringItem("/offset", str);

            if (Equals(str, "USER"))
            {
                dSettings.m_offsetOption = Configuration::DARK_MODEL_OPTION::USER_SUPPLIED;
            }
            else
            {
                dSettings.m_offsetOption = Configuration::DARK_MODEL_OPTION::MEASURED_IN_SCAN;
            }
            continue;
        }
    }

    return RETURN_CODE::FAIL;
}

}
