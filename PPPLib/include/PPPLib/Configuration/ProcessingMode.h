#pragma once

/** The options for how to do the post-processing.
    This also defines what types of operations will be performed */
enum class PROCESSING_MODE
{
    // Default mode, fliuxes are calculated from each successfully evaluated scan.
    PROCESSING_MODE_FLUX,

    // Calculation of molecular ratios in the plume, e.g. BrO/SO2. Not fully implemented.
    PROCESSING_MODE_COMPOSITION,


    PROCESSING_MODE_STRATOSPHERE,

    // Calculation of plume heights and plume directions.
    PROCESSING_MODE_GEOMETRY,

    // Performing instrument calibrations only, useful for prepairing for a later evaluation run.
    PROCESSING_MODE_INSTRUMENT_CALIBRATION
};
