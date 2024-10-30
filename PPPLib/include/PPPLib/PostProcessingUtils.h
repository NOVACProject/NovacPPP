#pragma once

#include <string>

namespace Configuration
{
class CNovacPPPConfiguration;
}

namespace novac
{
class ILogger;

/** Prepares for the evaluation of the spectra by reading in all the reference files that are needed.
    This will modify the contents of the setup to contain the references.
    @throws std::invalid_argument if the references files could not be found or not be read. */
void PrepareEvaluation(novac::ILogger& logger, std::string tempDirectory, Configuration::CNovacPPPConfiguration& setup);

}