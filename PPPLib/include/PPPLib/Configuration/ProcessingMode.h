#pragma once

/** The options for how to do the post-processing.
    This also defines what types of operations will be performed */
enum class ProcessingMode
{
    // Default mode, fliuxes are calculated from each successfully evaluated scan.
    Flux,

    // Calculation of molecular ratios in the plume, e.g. BrO/SO2. Not fully implemented.
    Composition,

    // Calculation of stratospheric composition. Not fully implemented
    Stratosphere,

    // Calculation of plume heights and plume directions.
    Geometry,

    // Performing instrument calibrations only, useful for prepairing for a later evaluation run.
    InstrumentCalibration
};
