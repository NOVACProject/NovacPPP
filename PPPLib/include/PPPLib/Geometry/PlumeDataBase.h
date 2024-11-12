#pragma once

#include "PlumeHeight.h"
#include "GeometryResult.h"
#include <SpectralEvaluation/DateTime.h>
#include <PPPLib/MFC/CString.h>
#include <PPPLib/Configuration/UserConfiguration.h>

// include the list-template from the C++ standard library
#include <list>

namespace Geometry
{
/** An instance of the class <b>CPlumeDataBase</b> can be used to
    keep track of the height of the gas plume over time.

    Potentially, this class can also be used to keep track of other
        properties of the plume, such as dispersion, rise, etc...
*/
class CPlumeDataBase
{
public:
    CPlumeDataBase(const Configuration::CUserConfiguration& userSettings);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The name of the volcano for which the database is valid. */
    novac::CString m_volcano;

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Retrieves the plume height at a given time.
        @param time - the time for which the wind field should be retrieved
        @param plumeHeight - will on successful return be filled with the information
            about the plume height at the requested time.
        @return true if the wind field could be retrieved, otherwise false.
     */
    bool GetPlumeHeight(const novac::CDateTime& time, PlumeHeight& plumeHeight) const;

    /** Inserts a plume height into the database */
    void InsertPlumeHeight(const PlumeHeight& plumeHeight);

    /** Inserts a calculated plume height into the database */
    void InsertPlumeHeight(const CGeometryResult& geomResult);

    /** Writes the contents of this database to file.
        @return 0 on success. */
    int WriteToFile(const novac::CString& fileName) const;

private:

    const Configuration::CUserConfiguration& m_userSettings;

    // the structure PlumeData is used to hold the information about the plume for a
    //  single point in time
    class PlumeData
    {
    public:
        PlumeData();
        PlumeData(const PlumeData& p);
        ~PlumeData();
        PlumeData& operator=(const PlumeData& p);

        /** The data needs to be labelled with the time when it is valid */
        novac::CDateTime validFrom;
        novac::CDateTime validTo;

        // The plume altitude (meters above sea level)
        double altitude;
        double altitudeError;
        Meteorology::MeteorologySource altitudeSource;
    };

    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    /** This is the database of plume information. Each item in the list
        holds the information of the plume for a single time frame.

        Each PlumeData object in the list MUST have an unique time frame.
        */
    std::list <PlumeData> m_dataBase;


    // ----------------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------------
    // ----------------------------------------------------------------------


    // Calculates the average and error of the plume heights in the given list
    void CalculateAverageHeight(const std::list <PlumeData>& plumeList, double& averageAltitude, double& altitudeError) const;

};
}