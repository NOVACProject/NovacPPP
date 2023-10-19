#ifndef METSOURCE_H
#define METSOURCE_H

namespace Meteorology {
    /**
        The data-type MET_SOURCE is used to identify where a piece of meteorological
        information comes from. These sources have different kind of 'quality' i.e.
        they are judged to not be equally good, e.g. is a wind-direction determined
        by a geometrical calculation determined to be more accurate than a wind-direction
        from a met. model. To retrieve the quality of a MET_SOURCE call the global
        function 'GetSourceQuality' defined in WindField.h
    */

    enum MET_SOURCE {
        MET_NONE,               /** no source - used to say that this contains no information */
        MET_DEFAULT,            /** default value */
        MET_USER,                /** this data was given by the user somehow */
        MET_ECMWF_FORECAST,        /** this data originates from a forecast by the ECMWF */
        MET_ECMWF_ANALYSIS,        /** this data originates from analysis data from the ECMWF */
        MET_GEOMETRY_CALCULATION,    /** this data was obtained by a geometrical calculation involving _two_ scanning instruments */
        MET_GEOMETRY_CALCULATION_SINGLE_INSTR,/** this data was obtained by a geometrical calculation involving _one_ scanning instrument */
        MET_DUAL_BEAM_MEASUREMENT,    /** this is a dual-beam measurement */
        MET_MODEL_WRF,                /** this data originates from a run of the WRF met. model */
        MET_NOAA_GDAS,                /** this data originates from the GDAS model from NOAA */
        MET_NOAA_FNL,                /** this data originates from the FNL model from NOAA */
        MET_NUMBER_OF_DEFINED_SOURCES    /** This is not a source, but is instead used to determine how many possible sources we have */
    };
}

#endif