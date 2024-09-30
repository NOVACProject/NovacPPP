#include <PPPLib/File/SetupFileReader.h>
#include <PPPLib/Exceptions.h>
#include <PPPLib/Logging.h>
#include <cstring>

using namespace FileHandler;
using namespace novac;

void CSetupFileReader::ReadSetupFile(const novac::CString& filename, Configuration::CNovacPPPConfiguration& setup)
{
    // 1. Open the file
    if (!Open(filename))
    {
        std::string message = std::string("Failed to open setup file for reading: '") + filename.std_str() + std::string("'");
        m_log.Error(message);
        throw PPPLib::FileIoException(message);
    }

    // parse the file, one line at a time.
    szToken = "start";
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        //* Look for the xml tag 'instrument' and use Parse_Instrument and Parse_Location to read serial number and location to object 'setup' */
        if (novac::Equals(szToken, "instrument", 9))
        {
            Configuration::CInstrumentConfiguration instrument;
            Parse_Instrument(instrument);
            setup.m_instrument.push_back(instrument);
            continue;
        }


    }//end while

    Close();
}

void CSetupFileReader::Parse_Instrument(Configuration::CInstrumentConfiguration& instr)
{

    // parse the file, one line at a time.
    szToken = "start";
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        if (novac::Equals(szToken, "/instrument", 11))
        {
            return;
        }

        if (novac::Equals(szToken, "serial", 6))
        {
            Parse_StringItem("/serial", instr.m_serial);
            continue;
        }
        if (novac::Equals(szToken, "location", 8))
        {
            Parse_Location(instr.m_location);
        }
    }
}

// Parse for location area and store in the LocationConfiguration object
void CSetupFileReader::Parse_Location(Configuration::CLocationConfiguration& loc)
{

    Configuration::CInstrumentLocation location;

    // parse the file, one line at a time.
    szToken = "start";
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        if (novac::Equals(szToken, "/location", 9))
        {
            // insert the information that we have into the array and return.
            loc.InsertLocation(location);
            return;
        }

        if (novac::Equals(szToken, "name", 4))
        {
            Parse_StringItem("/name", location.m_locationName);
            continue;
        }
        if (novac::Equals(szToken, "latitude", 8))
        {
            Parse_FloatItem("/latitude", location.m_latitude);
            continue;
        }
        if (Equals(szToken, "longitude", 9))
        {
            Parse_FloatItem("/longitude", location.m_longitude);
            continue;
        }
        if (Equals(szToken, "altitude", 8))
        {
            Parse_IntItem("/altitude", location.m_altitude);
            continue;
        }
        if (Equals(szToken, "compass", 7))
        {
            Parse_FloatItem("/compass", location.m_compass);
            continue;
        }
        if (Equals(szToken, "coneangle", 9))
        {
            Parse_FloatItem("/coneangle", location.m_coneangle);
            continue;
        }
        if (Equals(szToken, "tilt", 4))
        {
            Parse_FloatItem("/tilt", location.m_tilt);
            continue;
        }
        if (Equals(szToken, "type", 4))
        {
            Parse_IntItem("/type", (int&)location.m_instrumentType);
            continue;
        }
        if (Equals(szToken, "spectrometerModel", 17))
        {
            SpectrometerModel newSpectrometer;
            Parse_CustomSpectrometer(newSpectrometer);
            if (!newSpectrometer.IsUnknown())
            {
                CSpectrometerDatabase::GetInstance().AddModel(newSpectrometer);
            }
            continue;
        }
        if (Equals(szToken, "spec", 4))
        {
            Parse_StringItem("/spec", location.m_spectrometerModel);
            continue;
        }
        if (Equals(szToken, "volcano", 7))
        {
            Parse_StringItem("/volcano", location.m_volcano);
            continue;
        }
        if (Equals(szToken, "time_from", 9))
        {
            Parse_Date("/time_from", location.m_validFrom);
            continue;
        }
        if (Equals(szToken, "time_to", 7))
        {
            Parse_Date("/time_to", location.m_validTo);
            continue;
        }
    }
}

/** This takes care of writing the contents of a setup data-structure to file
    Only the part regarding the instrument's location will be written to the file */
RETURN_CODE CSetupFileReader::WriteSetupFile(const novac::CString& fileName, const Configuration::CNovacPPPConfiguration& setup)
{
    Configuration::CInstrumentLocation instrLocation;

    // Open the file
    FILE* f = fopen(fileName, "w");
    if (f == nullptr)
        return RETURN_CODE::FAIL;

    // Write the header
    fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
    fprintf(f, "<!-- This is the global configuration file for the evaluation of spectra in the NOVAC Post Processing Program -->\n");

    // the main item
    fprintf(f, "<NovacPPPConfiguration>\n");

    // loop through each instrument
    for (int k = 0; k < setup.NumberOfInstruments(); ++k)
    {
        const Configuration::CInstrumentConfiguration& instr = setup.m_instrument[k];

        fprintf(f, "\t<instrument>\n");

        // The serial-number
        fprintf(f, "\t\t<serial>%s</serial>\n", (const char*)instr.m_serial);

        // loop through the locations
        for (unsigned int j = 0; j < instr.m_location.GetLocationNum(); ++j)
        {
            instr.m_location.GetLocation(j, instrLocation);

            // The location
            fprintf(f, "\t\t<location>\n");

            fprintf(f, "\t\t\t<name>%s</name>\n", (const char*)instrLocation.m_locationName);
            fprintf(f, "\t\t\t<volcano>%s</volcano>\n", (const char*)instrLocation.m_volcano);
            fprintf(f, "\t\t\t<time_from>%04d.%02d.%02d</time_from>\n", instrLocation.m_validFrom.year, instrLocation.m_validFrom.month, instrLocation.m_validFrom.day);
            fprintf(f, "\t\t\t<time_to>%04d.%02d.%02d</time_to>\n", instrLocation.m_validTo.year, instrLocation.m_validTo.month, instrLocation.m_validTo.day);
            fprintf(f, "\t\t\t<latitude>%.5lf</latitude>\n", instrLocation.m_latitude);
            fprintf(f, "\t\t\t<longitude>%.5lf</longitude>\n", instrLocation.m_longitude);
            fprintf(f, "\t\t\t<altitude>%d</altitude>\n", instrLocation.m_altitude);
            fprintf(f, "\t\t\t<compass>%.1lf</compass>\n", instrLocation.m_compass);
            fprintf(f, "\t\t\t<coneangle>%.0lf</coneangle>\n", instrLocation.m_coneangle);
            fprintf(f, "\t\t<tilt>%.0lf</tilt>\n", instrLocation.m_tilt);
            fprintf(f, "\t\t<type>%d</type>\n", (int)instrLocation.m_instrumentType);

            const SpectrometerModel currentSpectrometer = CSpectrometerDatabase::GetInstance().GetModel(instrLocation.m_spectrometerModel);
            {
                fprintf(f, "\t\t<spectrometerModel>\n");
                fprintf(f, "\t\t\t<name>%s</name>\n", currentSpectrometer.modelName.c_str());
                fprintf(f, "\t\t\t<maxIntensity>%lf</maxIntensity>\n", currentSpectrometer.maximumIntensityForSingleReadout);
                fprintf(f, "\t\t\t<numberOfPixels>%d</numberOfPixels>\n", currentSpectrometer.numberOfPixels);
                fprintf(f, "\t\t</spectrometerModel>\n");
            }

            fprintf(f, "\t\t</location>\n");
        }

        fprintf(f, "\t</instrument>\n");
    }

    fprintf(f, "</NovacPPPConfiguration>\n");

    // remember to close the file
    fclose(f);

    return RETURN_CODE::SUCCESS;
}

void CSetupFileReader::Parse_CustomSpectrometer(SpectrometerModel& model)
{
    // parse the file, one line at a time.
    szToken = "start";
    while (szToken != nullptr)
    {
        szToken = NextToken();

        // no use to parse empty lines
        if (szToken == nullptr || strlen(szToken) < 3)
            continue;

        if (novac::Equals(szToken, "/spectrometerModel", 8))
        {
            return;
        }

        if (novac::Equals(szToken, "name", 4))
        {
            Parse_StringItem("/name", model.modelName);
            continue;
        }

        if (novac::Equals(szToken, "maxIntensity", 12))
        {
            Parse_FloatItem("/maxIntensity", model.maximumIntensityForSingleReadout);
            continue;
        }

        if (novac::Equals(szToken, "numberOfPixels", 14))
        {
            Parse_IntItem("/numberOfPixels", model.numberOfPixels);
            continue;
        }
    }
}