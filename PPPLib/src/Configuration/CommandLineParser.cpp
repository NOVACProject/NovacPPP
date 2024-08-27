#include <PPPLib/Configuration/CommandLineParser.h>
#include <PPPLib/Configuration/UserConfiguration.h>
#include <PPPLib/File/Filesystem.h>
#include <PPPLib/VolcanoInfo.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/MFC/CStringTokenizer.h>
#include <SpectralEvaluation/StringUtils.h>

namespace Configuration
{
    void CommandLineParser::ParseCommandLineOptions(
        const std::vector<std::string>& arguments,
        Configuration::CUserConfiguration& userSettings,
        novac::CVolcanoInfo& volcanoes,
        std::string& exePath)
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
                    errorMessage.Format("Could not parse date: %s", (const char*)parameter);
                    ShowMessage(errorMessage);
                }
                else
                {
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
                    errorMessage.Format("Could not parse date: %s", (const char*)parameter);
                    ShowMessage(errorMessage);
                }
                else
                {
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
                    userSettings.m_volcano = volcano;
                }
                else
                {
                    errorMessage.Format("Could not find volcano: %s", (const char*)parameter);
                    ShowMessage(errorMessage);
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
                }
                token = tokenizer.NextToken();
                continue;
            }

            // the maximum number of threads
            if (novac::Equals(currentToken, FLAG(str_maxThreadNum), strlen(FLAG(str_maxThreadNum))))
            {
                sscanf(currentToken.c_str() + strlen(FLAG(str_maxThreadNum)), "%ld", &userSettings.m_maxThreadNum);
                userSettings.m_maxThreadNum = std::max(userSettings.m_maxThreadNum, (unsigned long)1);
                token = tokenizer.NextToken();
                continue;
            }

            // The options for the local directory
            if (novac::Equals(currentToken, FLAG(str_includeSubDirectories_Local), strlen(FLAG(str_includeSubDirectories_Local))))
            {
                sscanf(currentToken.c_str() + strlen(FLAG(str_includeSubDirectories_Local)), "%d", &userSettings.m_includeSubDirectories_Local);
                token = tokenizer.NextToken();
                continue;
            }
            if (novac::Equals(currentToken, FLAG(str_LocalDirectory), strlen(FLAG(str_LocalDirectory))))
            {
                if (sscanf(currentToken.c_str() + strlen(FLAG(str_LocalDirectory)), "%s", buffer.data()))
                {
                    userSettings.m_LocalDirectory.Format("%s", buffer.data());
                    userSettings.m_LocalDirectory = Filesystem::AppendPathSeparator(userSettings.m_LocalDirectory);
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
                sscanf(currentToken.c_str() + strlen(FLAG(str_includeSubDirectories_FTP)), "%d", &userSettings.m_includeSubDirectories_FTP);
                token = tokenizer.NextToken();
                continue;
            }

            if (novac::Equals(currentToken, FLAG(str_FTPDirectory), strlen(FLAG(str_FTPDirectory))))
            {
                if (sscanf(currentToken.c_str() + strlen(FLAG(str_FTPDirectory)), "%s", buffer.data()))
                {
                    userSettings.m_FTPDirectory.Format("%s", buffer.data());
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
                    userSettings.m_FTPUsername.Format("%s", buffer.data());
                }
                token = tokenizer.NextToken();
                continue;
            }
            if (novac::Equals(currentToken, FLAG(str_FTPPassword), strlen(FLAG(str_FTPPassword))))
            {
                if (sscanf(currentToken.c_str() + strlen(FLAG(str_FTPPassword)), "%s", buffer.data()))
                {
                    userSettings.m_FTPPassword.Format("%s", buffer.data());
                }
                token = tokenizer.NextToken();
                continue;
            }

            // If we should upload the results to the NovacFTP server at the end...
            if (novac::Equals(currentToken, FLAG(str_uploadResults), strlen(FLAG(str_uploadResults))))
            {
                sscanf(currentToken.c_str() + strlen(FLAG(str_uploadResults)), "%d", &userSettings.m_uploadResults);
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
                }
                token = tokenizer.NextToken();
                continue;
            }

            // The processing mode
            if (novac::Equals(currentToken, FLAG(str_processingMode), strlen(FLAG(str_processingMode))))
            {
                sscanf(currentToken.c_str() + strlen(FLAG(str_processingMode)), "%d", (int*)&userSettings.m_processingMode);
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
                        userSettings.m_molecule = MOLEC_BRO;
                    }
                    else if (novac::Equals(moleculeName, "NO2"))
                    {
                        userSettings.m_molecule = MOLEC_NO2;
                    }
                    else if (novac::Equals(moleculeName, "O3"))
                    {
                        userSettings.m_molecule = MOLEC_O3;
                    }
                    else {
                        userSettings.m_molecule = MOLEC_SO2;
                    }
                }
                token = tokenizer.NextToken();
                continue;
            }

            // get the next token
            token = tokenizer.NextToken();
        }
    }
}