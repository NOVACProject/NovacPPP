#ifndef METSOURCE_H
#define METSOURCE_H

namespace Meteorology
{

/**
    The data-type MeteorologySource is used to identify where a piece of meteorological
    information comes from. These sources have different kind of 'quality' i.e.
    they are judged to not be equally good, e.g. is a wind-direction determined
    by a geometrical calculation determined to be more accurate than a wind-direction
    from a met. model. To retrieve the quality of a MeteorologySource call the global
    function 'GetSourceQuality' defined in WindField.h */

enum class MeteorologySource
{
    None,               /** no source - used to say that this contains no information */
    Default,            /** default value */
    User,                /** this data was given by the user somehow */
    EcmwfForecast,        /** this data originates from a forecast by the ECMWF */
    EcmwfAnalysis,        /** this data originates from analysis data from the ECMWF */
    GeometryCalculationTwoInstruments,    /** this data was obtained by a geometrical calculation involving _two_ scanning instruments */
    GeometryCalculationSingleInstrument,/** this data was obtained by a geometrical calculation involving _one_ scanning instrument */
    DualBeamMeasurement,    /** this is a dual-beam measurement */
    ModelledWrf,                /** this data originates from a run of the WRF met. model */
    NoaaGdas,                /** this data originates from the GDAS model from NOAA */
    NoaaFnl,                /** this data originates from the FNL model from NOAA */
    NumberOfSources    /** This is not a source, but is instead used to determine how many possible sources we have */
};
}

#endif