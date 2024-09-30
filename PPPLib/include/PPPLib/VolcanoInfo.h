#ifndef PPP_LIB_VOLCANO_INFO_H
#define PPP_LIB_VOLCANO_INFO_H

#include <vector>
#include <PPPLib/MFC/CString.h>

/** The <b>CVolcanoInfo</b>-class is a class that stores known information
        about a set of volcanoes. This information can then later be used in the program
        for various purposes.*/
namespace novac
{
class CVolcanoInfo
{
public:
    CVolcanoInfo();

    ~CVolcanoInfo() = default;

    // -----------------------------------------------------------
    // ---------------------- PUBLIC DATA ------------------------
    // -----------------------------------------------------------

    /** The number of volcanoes that are configured */
    unsigned int m_volcanoNum;

    /** This is the number of volcanoes that are automatically
        configured. Any additional volcanoes are configured
        by the user... */
    unsigned int m_preConfiguredVolcanoNum;

    // -------------------------------------------------------------
    // --------------------- PUBLIC METHODS ------------------------
    // -------------------------------------------------------------

    /** Adds a new volcano to the list */
    void AddVolcano(const novac::CString& name, const novac::CString& number, const novac::CString& country, double latitude, double longitude, double altitude, double hoursToGMT = 0.0, int observatory = 1);
    void UpdateVolcano(unsigned int index, const novac::CString& name, const novac::CString& number, const novac::CString& country, double latitude, double longitude, double altitude, double hoursToGMT = 0.0, int observatory = 1);

    /** Retrieves the name of the volcano with the given index */
    void GetVolcanoName(unsigned int index, novac::CString& name);

    /** Retrieves the code of the volcano with the given index */
    void GetVolcanoCode(unsigned int index, novac::CString& code);
    const novac::CString GetVolcanoCode(unsigned int index);

    /** Retrieves the location of the volcano */
    void GetVolcanoLocation(unsigned int index, novac::CString& location) const;
    novac::CString GetVolcanoLocation(unsigned int index) const;

    /** Retrieves the simplified name of the volcano with the given index */
    void GetSimpleVolcanoName(unsigned int index, novac::CString& name) const;
    novac::CString GetSimpleVolcanoName(unsigned int index) const;

    /** Retrieves the volcano index from a given name (or code) */
    int GetVolcanoIndex(const novac::CString& name);

    /** Retrieves the volcano position from the given index */
    double GetPeakLatitude(int index);
    double GetPeakLatitude(const novac::CString& name) { return GetPeakLatitude(GetVolcanoIndex(name)); }
    double GetPeakLongitude(int index);
    double GetPeakLongitude(const novac::CString& name) { return GetPeakLongitude(GetVolcanoIndex(name)); }
    double GetPeakAltitude(int index);
    double GetPeakAltitude(const novac::CString& name) { return GetPeakAltitude(GetVolcanoIndex(name)); }

    /** Retrieves the time-zone this volcano is in */
    double GetHoursToGMT(int index);
    double GetHoursToGMT(const novac::CString& name) { return GetHoursToGMT(GetVolcanoIndex(name)); }

    /** Retrieves the observatory that monitors this volcano */
    int GetObservatoryIndex(int index);
    int GetObservatoryIndex(const novac::CString& name) { return GetObservatoryIndex(GetVolcanoIndex(name)); }

private:
    class CVolcano
    {
    public:
        CVolcano();
        CVolcano(const novac::CString& name, const novac::CString& number, const novac::CString& country, double latitude, double longitude, double altitude, double hoursToGMT = 0.0, int observatory = 1);
        CVolcano(const novac::CString& name, const novac::CString& simpleName, const novac::CString& number, const novac::CString& country, double latitude, double longitude, double altitude, double hoursToGMT = 0.0, int observatory = 1);
        ~CVolcano() = default;

        /** The name of the volcano */
        novac::CString m_name;

        /** The simplified name of the volcano */
        novac::CString m_simpleName;

        /** The number of the volcano. This is from the
            Smithsonian's inventory of the worlds volcanoes
            http://www.volcano.si.edu */
        novac::CString m_number;

        /** The country where this volcano is located */
        novac::CString m_country;

        /** The latitude of the peak(s) */
        double m_peakLatitude;

        /** The longitude of the peak(s) */
        double m_peakLongitude;

        /** The altitude of the peak(s) (masl) */
        double m_peakHeight;

        /** The number of hours to GMT, used to calculate the local-time from the GPS-time */
        double m_hoursToGMT;

        /** The observatory in charge of this volcano */
        int m_observatory;
    };

    // ----------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------
    // ----------------------------------------------------------------

    /** The list of volcanoes that belongs to this CVolcanoInfo object */
    std::vector<CVolcano> m_volcanoes;

    // ----------------------------------------------------------------
    // --------------------- PRIVATE METHODS --------------------------
    // ----------------------------------------------------------------

    /** Fills in all the fields into the database */
    void InitializeDatabase();

    void InitializeDatabase_01();
    void InitializeDatabase_02();
    void InitializeDatabase_03();
    void InitializeDatabase_04();
    void InitializeDatabase_05();
    void InitializeDatabase_06();
    void InitializeDatabase_07();
    void InitializeDatabase_08();
    void InitializeDatabase_09();
    void InitializeDatabase_10();
    void InitializeDatabase_11();
    void InitializeDatabase_12();
    void InitializeDatabase_13();
    void InitializeDatabase_14();
    void InitializeDatabase_15();
    void InitializeDatabase_16();
    void InitializeDatabase_17();
    void InitializeDatabase_18();
    void InitializeDatabase_19();
};
}

#endif  // PPP_LIB_VOLCANO_INFO_H