#pragma once

/** The class CInstrumentCalibrationConfiguration defines per-instrument settings
    on how to perform calibrations. This will only be used in calibration mode. */
namespace Configuration
{
class CInstrumentCalibrationConfiguration
{
public:

    /** Path to the intial calibration file (either .std, .clb or .xs).
        If this is a file in the extended std format then it may also contain the instrument line shape
        (and hence make the instrumentLineshapeFile unnecessary). */
    std::string m_initialCalibrationFile;

    /** Path to the initial instrument line shape file (.slf) if any is provided.
        Ususally not set if m_initialCalibrationFile is .std. */
    std::string m_instrumentLineshapeFile;

};
}
