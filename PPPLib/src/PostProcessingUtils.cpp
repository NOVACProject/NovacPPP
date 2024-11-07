#include <PPPLib/PostProcessingUtils.h>
#include <SpectralEvaluation/Log.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <PPPLib/Configuration/NovacPPPConfiguration.h>

#include <PPPLib/File/Filesystem.h>

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace novac
{

static void ConvolveReference(novac::ILogger& logger, novac::LogContext context, std::string exePath, std::string tempFilePath, novac::CReferenceFile& ref)
{
    // Make sure the high-res section do exist.
    if (!IsExistingFile(ref.m_crossSectionFile))
    {
        const std::string fullPath = Filesystem::GetAbsolutePathFromRelative(ref.m_crossSectionFile, exePath);
        if (IsExistingFile(fullPath))
        {
            ref.m_crossSectionFile = fullPath;
        }
        else
        {
            throw InvalidReferenceException("Cannot create reference file for '" + ref.Name() + "'. Could not find given cross section file: " + ref.m_crossSectionFile);
        }
    }

    // Make sure the slit-function do exist.
    if (!IsExistingFile(ref.m_slitFunctionFile))
    {
        const std::string fullPath = Filesystem::GetAbsolutePathFromRelative(ref.m_slitFunctionFile, exePath);
        if (IsExistingFile(fullPath))
        {
            ref.m_slitFunctionFile = fullPath;
        }
        else
        {
            throw InvalidReferenceException("Cannot create reference file for '" + ref.Name() + "'. Could not find given slit function file: " + ref.m_slitFunctionFile);
        }
    }

    // Make sure the wavelength calibration do exist.
    if (!IsExistingFile(ref.m_wavelengthCalibrationFile))
    {
        const std::string fullPath = Filesystem::GetAbsolutePathFromRelative(ref.m_wavelengthCalibrationFile, exePath);
        if (IsExistingFile(fullPath))
        {
            ref.m_wavelengthCalibrationFile = fullPath;
        }
        else
        {
            throw InvalidReferenceException("Cannot create reference file for '" + ref.Name() + "'. Could not find given wavelength calibration  file: " + ref.m_wavelengthCalibrationFile);
        }
    }

    // Now do the convolution
    logger.Information(context, "Convolving reference.");
    if (ref.ConvolveReference())
    {
        throw InvalidReferenceException("Cannot create reference file for '" + ref.Name() + "'. Convolution failed.");
    }

    // Save the resulting reference, for reference...
    SaveCrossSectionFile(tempFilePath, *ref.m_data);
}


void PrepareEvaluation(novac::ILogger& logger, std::string tempDirectory, Configuration::CNovacPPPConfiguration& setup)
{
    logger.Information("--- Reading References --- ");

    if (setup.m_instrument.size() == 0)
    {
        throw std::invalid_argument("No instruments were configured.");
    }

    novac::LogContext context;

    directorySetup directories;
    directories.tempDirectory = tempDirectory;
    directories.executableDirectory = setup.m_executableDirectory;

    // this is true if we failed to prepare the evaluation...
    bool failure = false;

    // Loop through each of the configured instruments
    for (size_t instrumentIndex = 0; instrumentIndex < setup.m_instrument.size(); ++instrumentIndex)
    {
        auto instrumentContext = context.With(novac::LogContext::Device, setup.m_instrument[instrumentIndex].m_serial.std_str());

        // For each instrument, loop through the fit-windows that are defined
        const size_t numberOfFitWindows = setup.m_instrument[instrumentIndex].m_eval.NumberOfFitWindows();
        for (size_t fitWindowIndex = 0; fitWindowIndex < numberOfFitWindows; ++fitWindowIndex)
        {
            novac::CFitWindow window;
            CDateTime fromTime, toTime; //  these are not used but must be passed onto GetFitWindow...
            if (setup.m_instrument[instrumentIndex].m_eval.GetFitWindow(fitWindowIndex, window, fromTime, toTime))
            {
                logger.Error(instrumentContext, "Failed to get fit window from configuration.");
                failure = true;
                continue;
            }

            PrepareFitWindow(logger, instrumentContext, setup.m_instrument[instrumentIndex].m_serial.std_str(), window, directories);

            // If we've made it this far, then we've managed to read in all the references.
            // Now store the data in setup
            setup.m_instrument[instrumentIndex].m_eval.SetFitWindow(fitWindowIndex, window, fromTime, toTime);
        }
    }

    if (failure)
    {
        throw std::invalid_argument("failed to setup evaluation");
    }
}

void PrepareFitWindow(novac::ILogger& logger, novac::LogContext& instrumentContext, const std::string& instrumentSerial, novac::CFitWindow& window, const directorySetup& setup)
{
    auto windowContext = instrumentContext.With(novac::LogContext::FitWindow, window.name);

    // For each reference in the fit-window, read it in and make sure that it exists...
    for (size_t referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex)
    {
        auto referenceContext = windowContext.With(novac::LogContext::FileName, window.ref[referenceIndex].m_path);

        if (window.ref[referenceIndex].m_path.empty())
        {
            // The reference file was not given in the configuration. Try to generate a configuration
            //  from the cross section, slit-function and wavelength calibration. These three must then 
            //  exist or the evaluation fails.
            novac::CString tempFile;
            tempFile.Format("%s%s_%s.xs", setup.tempDirectory.c_str(), instrumentSerial.c_str(), window.ref[referenceIndex].m_specieName.c_str());

            ConvolveReference(logger, referenceContext, setup.executableDirectory, tempFile.std_str(), window.ref[referenceIndex]);
        }
        else
        {
            if (!IsExistingFile(window.ref[referenceIndex].m_path))
            {
                // the file does not exist, try to change it to include the path of the configuration-directory...
                const std::string fileName = Filesystem::GetAbsolutePathFromRelative(window.ref[referenceIndex].m_path, setup.executableDirectory);
                if (IsExistingFile(fileName))
                {
                    window.ref[referenceIndex].m_path = fileName;
                }
                else
                {
                    throw InvalidReferenceException("Cannot find reference file: '" + fileName + "' defined for fit window :'" + window.name + "' for instrument: " + instrumentSerial);
                }
            }

            // Read in the cross section
            window.ref[referenceIndex].ReadCrossSectionDataFromFile();

            // Verify that the reference has values in the given range. Throws InvalidReferenceException if it doesn't.
            window.ref[referenceIndex].VerifyReferenceValues(window.fitLow, window.fitHigh);

            // Make a check of the data range as well, in order to show the user.
            {
                std::pair<size_t, size_t> indices;
                const auto minMaxValues = MinMax(
                    window.ref[referenceIndex].m_data->m_crossSection.begin() + window.fitLow,
                    window.ref[referenceIndex].m_data->m_crossSection.begin() + window.fitHigh,
                    indices);

                std::stringstream msg;
                msg << "Reference has values in range [" << minMaxValues.first << ", " << minMaxValues.second << "]";
                if (std::abs(minMaxValues.first) > 1e-6 || std::abs(minMaxValues.second) > 1e-6)
                {
                    msg << ". This seems large. Please verify that the reference is scaled to molecules/cm2.";
                }
                logger.Information(referenceContext, msg.str());
            }

            // If we are supposed to high-pass filter the spectra then
            // we should also high-pass filter the cross-sections
            if (window.fitType == novac::FIT_TYPE::FIT_HP_DIV || window.fitType == novac::FIT_TYPE::FIT_HP_SUB)
            {
                if (window.ref[referenceIndex].m_isFiltered == false)
                {
                    logger.Information(referenceContext, "High pass filtering reference.");
                    if (novac::Equals(window.ref[referenceIndex].m_specieName, "ring"))
                    {
                        HighPassFilter_Ring(*window.ref[referenceIndex].m_data);
                    }
                    else
                    {
                        HighPassFilter(*window.ref[referenceIndex].m_data, CrossSectionUnit::cm2_molecule);
                    }
                }
                else
                {
                    throw InvalidReferenceException("Reference file is filtered. This is not supported in the NovacPPP. Reference: '" + window.ref[referenceIndex].Name() + "' defined for fit window :'" + window.name + "' for instrument: " + instrumentSerial);
                }
            }
        }// endif
    }

    // If the window also contains a fraunhofer-reference then read it too.
    if (window.fraunhoferRef.m_path.size() > 4)
    {
        auto referenceContext = instrumentContext.With(novac::LogContext::FileName, window.fraunhoferRef.m_path);

        if (!IsExistingFile(window.fraunhoferRef.m_path))
        {
            // the file does not exist, try to change it to include the path of the configuration-directory...
            const std::string fileName = Filesystem::GetAbsolutePathFromRelative(window.fraunhoferRef.m_path, setup.executableDirectory);

            if (IsExistingFile(fileName))
            {
                window.fraunhoferRef.m_path = fileName;
            }
            else
            {
                throw InvalidReferenceException("Cannot find Fraunhofer reference file defined for fit window :'" + window.name + "' for instrument: " + instrumentSerial);
            }
        }

        window.fraunhoferRef.ReadCrossSectionDataFromFile();

        if (window.fitType == novac::FIT_TYPE::FIT_HP_DIV || window.fitType == novac::FIT_TYPE::FIT_HP_SUB)
        {
            logger.Information(referenceContext, "High pass filtering Fraunhofer reference.");
            HighPassFilter_Ring(*window.fraunhoferRef.m_data);
        }
        else
        {
            logger.Information(referenceContext, "Running log on Fraunhofer reference.");
            Log(*window.fraunhoferRef.m_data);
        }
    }
}

}