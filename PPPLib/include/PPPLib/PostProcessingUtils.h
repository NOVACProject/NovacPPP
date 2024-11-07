#pragma once

#include <string>

namespace Configuration
{
class CNovacPPPConfiguration;
}

namespace novac
{
class ILogger;
class LogContext;
class CFitWindow;

struct directorySetup
{
    std::string tempDirectory;
    std::string executableDirectory;
};

/** Prepares for the evaluation of the spectra by reading in all the reference files that are needed.
    This will modify the contents of the setup to contain the references.
    @throws std::invalid_argument if the references files could not be found or not be read. */
void PrepareEvaluation(novac::ILogger& logger, std::string tempDirectory, Configuration::CNovacPPPConfiguration& setup);

/** Prepares for the evaluation by reading in all the reference files in the fit window and filtering them when needed.
    @throws novac::InvalidReferenceException if any of the references files could not be found or not be read. */
void PrepareFitWindow(novac::ILogger& logger, novac::LogContext& instrumentContext, const std::string& instrumentSerial, novac::CFitWindow& window, const directorySetup& setup);

}