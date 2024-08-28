#pragma once

#include <vector>
#include <PPPLib/Logging.h>

namespace novac
{
    class CVolcanoInfo;
}

namespace Configuration
{
    class CUserConfiguration;

    /* The CommandLineParser is used to parse whatever options the user has given to the command line */
    class CommandLineParser
    {
    public:
        /**
         * @brief Parses the given command line options
         * @param arguments The options given from the command line.
         * @param userSettings The user configuration which will be filled with data from the arguments.
         * @param volcanoes The global list of volcanoes in the program, this is used to look up the name / number of the volcano.
         * @param exePath An optionally set path to the executable. This will be updated only if the parameter is set in the command line.
         */
        static void ParseCommandLineOptions(
            const std::vector<std::string>& arguments,
            Configuration::CUserConfiguration& userSettings,
            novac::CVolcanoInfo& volcanoes,
            std::string& exePath,
            ILogger& log);
    };
}