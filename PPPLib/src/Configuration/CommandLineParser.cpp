#include <PPPLib/Configuration/CommandLineParser.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/VolcanoInfo.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/MFC/CStringTokenizer.h>
#include <SpectralEvaluation/StringUtils.h>

#include <cstring>
#include <sstream>

namespace Configuration
{
void CommandLineParser::ParseCommandLineOptions(
    const std::vector<std::string>& arguments,
    Configuration::CUserConfiguration& userSettings,
    novac::CVolcanoInfo& volcanoes,
    std::string& exePath,
    novac::ILogger& log)
{
    char seps[] = " \t";
    std::vector<char> buffer(16384, 0);
    novac::CString parameter;
    novac::CString errorMessage;

    // a local copy of the command line parameters
    novac::CString commandLine;
    for (size_t arg = 0; arg < arguments.size(); ++arg)
    {
        commandLine = commandLine + " " + novac::CString(arguments[arg]);
    }

    novac::LogContext context;

    // Go through the received input parameters and set the approprate option
    novac::CStringTokenizer tokenizer(commandLine.c_str(), seps);
    const char* token = tokenizer.NextToken();

    while (nullptr != token)
    {
        std::string currentToken{ token };
        Trim(currentToken);

        // The first date which we should analyze data from
        if (novac::Equals(currentToken, FLAG(str_fromDate), strlen(FLAG(str_fromDate))))
        {
            parameter.Format(token + strlen(FLAG(str_fromDate)) + 1);

            novac::CDateTime parsedDate;
            if (!novac::CDateTime::ParseDate(parameter, parsedDate))
            {
                std::stringstream msg;
                msg << "Could not parse From date: " << parameter.std_str();
                throw std::invalid_argument(msg.str());
            }
            else
            {
                log.Information(context.With("cmd", str_fromDate), "Set From date: " + parameter.std_str());
                userSettings.m_fromDate = parsedDate;
            }

            token = tokenizer.NextToken();
            continue;
        }

        // The last date which we should analyze data from
        if (novac::Equals(currentToken, FLAG(str_toDate), strlen(FLAG(str_toDate))))
        {
            parameter.Format(token + strlen(FLAG(str_toDate)) + 1);

            novac::CDateTime parsedDate;
            if (!novac::CDateTime::ParseDate(parameter, parsedDate))
            {
                std::stringstream msg;
                msg << "Could not parse To date: " << parameter.std_str();
                throw std::invalid_argument(msg.str());
            }
            else
            {
                log.Information(context.With("cmd", str_toDate), "Set To date: " + parameter.std_str());
                userSettings.m_toDate = parsedDate;
            }

            token = tokenizer.NextToken();
            continue;
        }

        // the volcano to process data from 
        if (novac::Equals(currentToken, FLAG(str_volcano), strlen(FLAG(str_volcano))))
        {
            parameter.Format(token + strlen(FLAG(str_volcano)) + 1);

            int volcano = volcanoes.GetVolcanoIndex(parameter);
            if (volcano >= 0)
            {
                log.Information(context.With("cmd", str_volcano), "Set volcano: " + parameter.std_str());
                userSettings.m_volcano = volcano;
            }
            else
            {
                std::stringstream msg;
                msg << "Could not find volcano: " << parameter.std_str();
                throw std::invalid_argument(msg.str());
            }
            token = tokenizer.NextToken();
            continue;
        }

        // the working directory, used to override the location of the configurations
        if (novac::Equals(currentToken, FLAG(str_workingDirectory), strlen(FLAG(str_workingDirectory))))
        {
            if (sscanf(currentToken.c_str() + strlen(FLAG(str_workingDirectory)), "%s", buffer.data()))
            {
                exePath = std::string(buffer.data());
                Trim(exePath, " \t\"");
                exePath = Filesystem::AppendPathSeparator(exePath);

                log.Information(context.With("cmd", str_workingDirectory), "Set working directory: " + exePath);
            }
            token = tokenizer.NextToken();
            continue;
        }

        // the maximum number of threads
        if (novac::Equals(currentToken, FLAG(str_maxThreadNum), strlen(FLAG(str_maxThreadNum))))
        {
            if (1 == sscanf(currentToken.c_str() + strlen(FLAG(str_maxThreadNum)), "%ld", &userSettings.m_maxThreadNum))
            {
                log.Information(context.With("cmd", str_maxThreadNum), "Set max number of threads");
                userSettings.m_maxThreadNum = std::max(userSettings.m_maxThreadNum, (unsigned long)1);
            }
            token = tokenizer.NextToken();
            continue;
        }

        // The options for the local directory
        if (novac::Equals(currentToken, FLAG(str_includeSubDirectories_Local), strlen(FLAG(str_includeSubDirectories_Local))))
        {
            if (1 == sscanf(currentToken.c_str() + strlen(FLAG(str_includeSubDirectories_Local)), "%d", &userSettings.m_includeSubDirectories_Local))
            {
                log.Information(context.With("cmd", str_includeSubDirectories_Local), "Updated includeSubDirectories_Local");
            }
            token = tokenizer.NextToken();
            continue;
        }
        if (novac::Equals(currentToken, FLAG(str_LocalDirectory), strlen(FLAG(str_LocalDirectory))))
        {
            if (sscanf(currentToken.c_str() + strlen(FLAG(str_LocalDirectory)), "%s", buffer.data()))
            {
                userSettings.m_LocalDirectory.Format("%s", buffer.data());
                userSettings.m_LocalDirectory = Filesystem::AppendPathSeparator(userSettings.m_LocalDirectory);

                log.Information(context.With("cmd", str_LocalDirectory), "Set local directory: " + userSettings.m_LocalDirectory.std_str());
            }
            else
            {
                userSettings.m_LocalDirectory = "";
            }
            token = tokenizer.NextToken();
            continue;
        }

        // The options for the FTP directory
        if (novac::Equals(currentToken, FLAG(str_includeSubDirectories_FTP), strlen(FLAG(str_includeSubDirectories_FTP))))
        {
            if (1 == sscanf(currentToken.c_str() + strlen(FLAG(str_includeSubDirectories_FTP)), "%d", &userSettings.m_includeSubDirectories_FTP))
            {
                log.Information(context.With("cmd", str_includeSubDirectories_FTP), "Updated include FTP sub directories");
            }
            token = tokenizer.NextToken();
            continue;
        }

        if (novac::Equals(currentToken, FLAG(str_FTPDirectory), strlen(FLAG(str_FTPDirectory))))
        {
            if (sscanf(currentToken.c_str() + strlen(FLAG(str_FTPDirectory)), "%s", buffer.data()))
            {
                userSettings.m_FTPDirectory.Format("%s", buffer.data());
                log.Information(context.With("cmd", str_FTPDirectory), "Updated FTP directory: " + userSettings.m_FTPDirectory.std_str());
            }
            else
            {
                userSettings.m_FTPDirectory.Format("");
            }
            token = tokenizer.NextToken();
            continue;
        }

        if (novac::Equals(currentToken, FLAG(str_FTPUsername), strlen(FLAG(str_FTPUsername))))
        {
            if (sscanf(currentToken.c_str() + strlen(FLAG(str_FTPUsername)), "%s", buffer.data()))
            {
                log.Information(context.With("cmd", str_FTPUsername), "Updated FTP username");
                userSettings.m_FTPUsername.Format("%s", buffer.data());
            }
            token = tokenizer.NextToken();
            continue;
        }
        if (novac::Equals(currentToken, FLAG(str_FTPPassword), strlen(FLAG(str_FTPPassword))))
        {
            if (sscanf(currentToken.c_str() + strlen(FLAG(str_FTPPassword)), "%s", buffer.data()))
            {
                log.Information(context.With("cmd", str_FTPPassword), "Updated FTP password");
                userSettings.m_FTPPassword.Format("%s", buffer.data());
            }
            token = tokenizer.NextToken();
            continue;
        }

        // If we should upload the results to the NovacFTP server at the end...
        if (novac::Equals(currentToken, FLAG(str_uploadResults), strlen(FLAG(str_uploadResults))))
        {
            if (1 == sscanf(currentToken.c_str() + strlen(FLAG(str_uploadResults)), "%d", &userSettings.m_uploadResults))
            {
                log.Information(context.With("cmd", str_uploadResults), "Updated upload results");
            }
            token = tokenizer.NextToken();
            continue;
        }

        // The output directory
        if (novac::Equals(currentToken, FLAG(str_outputDirectory), strlen(FLAG(str_outputDirectory))))
        {
            if (sscanf(currentToken.c_str() + strlen(FLAG(str_outputDirectory)), "%s", buffer.data()))
            {
                userSettings.m_outputDirectory.Format("%s", buffer.data());
                userSettings.m_outputDirectory = Filesystem::AppendPathSeparator(userSettings.m_outputDirectory);

                log.Information(context.With("cmd", str_outputDirectory), "Updated output directory: " + userSettings.m_outputDirectory.std_str());
            }
            token = tokenizer.NextToken();
            continue;
        }

        // The temporary directory
        if (novac::Equals(currentToken, FLAG(str_tempDirectory), strlen(FLAG(str_tempDirectory))))
        {
            if (sscanf(currentToken.c_str() + strlen(FLAG(str_tempDirectory)), "%s", buffer.data()))
            {
                userSettings.m_tempDirectory.Format("%s", buffer.data());
                userSettings.m_tempDirectory = Filesystem::AppendPathSeparator(userSettings.m_tempDirectory);

                log.Information(context.With("cmd", str_tempDirectory), "Updated temp directory: " + userSettings.m_tempDirectory.std_str());
            }
            token = tokenizer.NextToken();
            continue;
        }

        // The windField file
        int N = (int)strlen(FLAG(str_windFieldFile));
        if (novac::Equals(currentToken, FLAG(str_windFieldFile), N))
        {
            if (sscanf(currentToken.c_str() + N, "%s", buffer.data()))
            {
                userSettings.m_windFieldFile.Format("%s", buffer.data());

                log.Information(context.With("cmd", str_windFieldFile), "Updated wind field file: " + userSettings.m_windFieldFile.std_str());
            }
            token = tokenizer.NextToken();
            continue;
        }

        // The processing mode
        if (novac::Equals(currentToken, FLAG(str_processingMode), strlen(FLAG(str_processingMode))))
        {
            if (1 == sscanf(currentToken.c_str() + strlen(FLAG(str_processingMode)), "%d", (int*)&userSettings.m_processingMode))
            {
                log.Information(context.With("cmd", str_processingMode), "Updated processing mode");
            }
            token = tokenizer.NextToken();
            continue;
        }

        // the molecule
        if (novac::Equals(currentToken, FLAG(str_molecule), strlen(FLAG(str_molecule))))
        {
            if (sscanf(currentToken.c_str() + strlen(FLAG(str_molecule)), "%s", buffer.data()))
            {
                const std::string moleculeName = std::string(buffer.data());
                if (novac::Equals(moleculeName, "BrO"))
                {
                    userSettings.m_molecule = StandardMolecule::BrO;
                }
                else if (novac::Equals(moleculeName, "NO2"))
                {
                    userSettings.m_molecule = StandardMolecule::NO2;
                }
                else if (novac::Equals(moleculeName, "O3"))
                {
                    userSettings.m_molecule = StandardMolecule::O3;
                }
                else
                {
                    userSettings.m_molecule = StandardMolecule::SO2;
                }
                log.Information(context.With("cmd", str_molecule), "Updated molecule");
            }
            token = tokenizer.NextToken();
            continue;
        }

        // get the next token
        token = tokenizer.NextToken();
    }
}
}