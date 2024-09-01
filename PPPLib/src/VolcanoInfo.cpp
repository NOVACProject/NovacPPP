#include <PPPLib/VolcanoInfo.h>

namespace novac
{

CVolcanoInfo::CVolcanoInfo()
{
    InitializeDatabase();
}

CVolcanoInfo::CVolcano::CVolcano()
{
    m_name = "";
    m_simpleName = "";
    m_number = "";
    m_country = "";

    m_peakLatitude = 0.0;
    m_peakLongitude = 0.0;
    m_peakHeight = 0.0;

    m_hoursToGMT = 0.0;

    m_observatory = 0;
}

CVolcanoInfo::CVolcano::CVolcano(const novac::CString& name, const novac::CString& number, const novac::CString& country, double latitude, double longitude, double altitude, double hoursToGMT, int observatory)
{
    m_name.Format(name);
    SimplifyString(name, m_simpleName);
    m_number.Format(number);
    m_country.Format(country);

    m_peakLatitude = latitude;
    m_peakLongitude = longitude;
    m_peakHeight = altitude;

    m_hoursToGMT = hoursToGMT;

    m_observatory = observatory;
}

CVolcanoInfo::CVolcano::CVolcano(const novac::CString& name, const novac::CString& simpleName, const novac::CString& number, const novac::CString& country, double latitude, double longitude, double altitude, double hoursToGMT, int observatory)
{
    m_name = name;

    if (simpleName.GetLength() == 0)
    {
        SimplifyString(name, m_simpleName);
    }
    else
    {
        m_simpleName = simpleName;
    }

    m_number = number;
    m_country = country;

    m_peakLatitude = latitude;
    m_peakLongitude = longitude;
    m_peakHeight = altitude;

    m_hoursToGMT = hoursToGMT;

    m_observatory = observatory;
}

/** Adds a new volcano to the list */
void CVolcanoInfo::AddVolcano(const novac::CString& name, const novac::CString& number, const novac::CString& country, double latitude, double longitude, double altitude, double hoursToGMT, int observatory)
{
    novac::CString simpleName;
    SimplifyString(name, simpleName);

    m_volcanoes.push_back(CVolcano(name, simpleName, number, country, latitude, longitude, altitude, hoursToGMT, observatory));

    ++m_volcanoNum;
}

void CVolcanoInfo::UpdateVolcano(unsigned int index, const novac::CString& name, const novac::CString& number, const novac::CString& country, double latitude, double longitude, double altitude, double hoursToGMT, int observatory)
{
    if (index >= this->m_volcanoNum)
    {
        AddVolcano(name, number, country, latitude, longitude, altitude, hoursToGMT, observatory);
    }
    else
    {
        CVolcano& volcano = m_volcanoes.at(index);
        volcano.m_name.Format(name);
        SimplifyString(name, volcano.m_simpleName);
        volcano.m_number.Format(number);
        volcano.m_country.Format(country);
        volcano.m_peakHeight = altitude;
        volcano.m_peakLatitude = latitude;
        volcano.m_peakLongitude = longitude;
        volcano.m_observatory = observatory;
        volcano.m_hoursToGMT = hoursToGMT;
    }
}

int CVolcanoInfo::GetVolcanoIndex(const novac::CString& name)
{
    static unsigned int lastIndex = 0;

    // first try with the same volcano as the last time
    //	this function was called...
    CVolcano& vol = m_volcanoes.at(lastIndex);
    if (Equals(vol.m_name, name) || Equals(vol.m_simpleName, name) || Equals(vol.m_number, name))
    {
        return (int)lastIndex;
    }

    // this was not the same volcano as last time, search for it...
    lastIndex = 0;
    for (std::vector<CVolcano>::iterator iter = m_volcanoes.begin(); iter != m_volcanoes.end(); ++iter)
    {
        if (Equals(name, iter->m_name) || Equals(name, iter->m_simpleName) || Equals(name, iter->m_number))
        {
            return (int)lastIndex;
        }
        ++lastIndex;
    }

    return -1; // no volcano found
}

void CVolcanoInfo::GetVolcanoName(unsigned int index, novac::CString& name)
{
    if (index >= m_volcanoNum)
    {
        name = "";
    }
    else
    {
        CVolcano& vol = m_volcanoes.at(index);
        name.Format(vol.m_name);
    }
}

novac::CString CVolcanoInfo::GetSimpleVolcanoName(unsigned int index) const
{
    novac::CString name;
    this->GetSimpleVolcanoName(index, name);
    return name;
}

/** Retrieves the location of the volcano */
void CVolcanoInfo::GetVolcanoLocation(unsigned int index, novac::CString& location) const
{
    if (index >= m_volcanoNum)
    {
        location = "";
    }
    else
    {
        const CVolcano& vol = m_volcanoes.at(index);
        location.Format(vol.m_country);
    }
}

novac::CString CVolcanoInfo::GetVolcanoLocation(unsigned int index) const
{
    novac::CString location;

    this->GetVolcanoLocation(index, location);

    return location;
}

void CVolcanoInfo::GetVolcanoCode(unsigned int index, novac::CString& code)
{
    if (index >= m_volcanoNum)
    {
        code = "";
    }
    else
    {
        CVolcano& vol = m_volcanoes.at(index);
        code.Format(vol.m_number);
    }
}

const novac::CString CVolcanoInfo::GetVolcanoCode(unsigned int index)
{
    novac::CString name;
    this->GetVolcanoCode(index, name);
    return name;
}

void CVolcanoInfo::GetSimpleVolcanoName(unsigned int index, novac::CString& name) const
{
    if (index >= m_volcanoNum)
    {
        name = "";
    }
    else
    {
        const CVolcano& vol = m_volcanoes.at(index);
        name.Format(vol.m_simpleName);
    }
}

double CVolcanoInfo::GetPeakLatitude(unsigned int index)
{
    if (index >= m_volcanoNum)
    {
        return 0.0;
    }
    else
    {
        CVolcano& vol = m_volcanoes.at(index);
        return vol.m_peakLatitude;
    }
}
double CVolcanoInfo::GetPeakLongitude(unsigned int index)
{
    if (index >= m_volcanoNum)
    {
        return 0.0;
    }
    else
    {
        CVolcano& vol = m_volcanoes.at(index);
        return vol.m_peakLongitude;
    }
}
double CVolcanoInfo::GetPeakAltitude(unsigned int index)
{
    if (index >= m_volcanoNum)
    {
        return 0.0;
    }
    else
    {
        CVolcano& vol = m_volcanoes.at(index);
        return vol.m_peakHeight;
    }
}

/** Retrieves the time-zone this volcano is in */
double CVolcanoInfo::GetHoursToGMT(unsigned int index)
{
    if (index >= m_volcanoNum)
    {
        return 0.0;
    }
    else
    {
        CVolcano& vol = m_volcanoes.at(index);
        return vol.m_hoursToGMT;
    }
}

/** Retrieves the observatory that monitors this volcano */
int CVolcanoInfo::GetObservatoryIndex(unsigned int index)
{
    if (index >= m_volcanoNum)
    {
        return 1;
    }
    else
    {
        CVolcano& vol = m_volcanoes.at(index);
        return vol.m_observatory;
    }
}


void CVolcanoInfo::InitializeDatabase()
{
    m_volcanoes.push_back(CVolcano("Arenal", "arenal", "1405-033", "Costa Rica", 10.46, -84.7, 1670, -6, 11));
    m_volcanoes.push_back(CVolcano("Poás", "poas", "1405-02=", "Costa Rica", 10.20, -84.23, 2708, -6, 11));
    m_volcanoes.push_back(CVolcano("Turrialba", "turrialba", "1405-07=", "Costa Rica", 10.025, -83.767, 3340, -6, 11));
    m_volcanoes.push_back(CVolcano("Santa Ana", "santa_ana", "1403-02=", "El Salvador", 13.85, -89.63, 2381, -6, 13));
    m_volcanoes.push_back(CVolcano("San Miguel", "san_miguel", "1403-10=", "El Salvador", 13.43, -88.27, 2130, -6, 13));
    m_volcanoes.push_back(CVolcano("Popocatépetl", "popocatepetl", "1401-09=", "México", 19.02, -98.62, 5426, -6, 17));
    m_volcanoes.push_back(CVolcano("Fuego de Colima", "fuego_de_colima", "1401-04=", "México", 19.51, -103.62, 3850, -6, 17));
    m_volcanoes.push_back(CVolcano("San Cristóbal", "san_cristobal", "1404-02=", "Nicaragua", 12.70, -87.0, 1745, -6, 4));
    m_volcanoes.push_back(CVolcano("Masaya", "masaya", "1404-10=", "Nicaragua", 11.98, -86.14, 635, -6, 4));
    m_volcanoes.push_back(CVolcano("Galeras", "galeras", "1501-08=", "Colombia", 1.22, -77.37, 4276, -5, 5));
    m_volcanoes.push_back(CVolcano("Nevado del Ruiz", "nevado_del_ruiz", "1501-02=", "Colombia", 4.89, -75.32, 5321, -5, 5));
    m_volcanoes.push_back(CVolcano("Nevado del Huila", "nevado_del_huila", "1501-05=", "Colombia", 2.93, -76.03, 5364, -5, 5));
    m_volcanoes.push_back(CVolcano("Nyiragongo", "nyiragongo", "0203-03=", "Democratic Republic of Congo", -1.516732, 29.24668, 3470, +2, 11));
    m_volcanoes.push_back(CVolcano("Nyamuragira", "nyamuragira", "0203-02=", "Democratic Republic of Congo", -1.41, 29.20, 3058, +2, 11));
    m_volcanoes.push_back(CVolcano("Etna", "etna", "0101-06=", "Italy", 37.752, 14.995, 3300, +1, 6));
    m_volcanoes.push_back(CVolcano("La Soufrière", "la_soufriere", "1600-06=", "France", 16.05, -61.67, 1467, 0, 9));
    m_volcanoes.push_back(CVolcano("Tungurahua", "tungurahua", "1502-08=", "Ecuador", -1.467, -78.442, 5023, -5, 3));
    m_volcanoes.push_back(CVolcano("Cotopaxi", "cotopaxi", "1502-05=", "Ecuador", -0.677, -78.436, 5911, -5, 3));
    m_volcanoes.push_back(CVolcano("Pacaya", "pacaya", "1402-11=", "Guatemala", 14.381, -90.601, 2552, -6, 8));
    m_volcanoes.push_back(CVolcano("Piton de la Fournaise", "piton_de_la_fournaise", "0303-02=", "France", -21.23, 55.71, 2632, +4, 9));
    m_volcanoes.push_back(CVolcano("Santiaguito", "santiaguito", "1402-03=", "Guatemala", 14.756, -91.552, 3772, -6, 8));
    m_volcanoes.push_back(CVolcano("Fuego (Guatemala)", "fuego_guatemala", "1402-09=", "Guatemala", 14.473, -90.880, 3763, -6, 8));
    m_volcanoes.push_back(CVolcano("Vulcano", "vulcano", "0101-05=", "Italy", 38.404, 14.986, 500, +1, 7));
    m_volcanoes.push_back(CVolcano("Stromboli", "stromboli", "0101-04=", "Italy", 38.789, 15.213, 924, +1, 7));
    m_volcanoes.push_back(CVolcano("Yasur", "yasur", "0507-10=", "Indonesia", -19.53, 169.442, 361, +11, 2));
    m_volcanoes.push_back(CVolcano("Villarrica", "villarrica", "1507-12=", "Chile", -39.42, -71.93, 2847, -7, 2));
    m_volcanoes.push_back(CVolcano("Llaima", "llaima", "1507-11=", "Chile", -38.692, -71.729, 3125, -7, 2));
    m_volcanoes.push_back(CVolcano("Antisana", "antisana", "1502-03=", "Ecuador", -0.481, -78.141, 5753, -5, 3));
    m_volcanoes.push_back(CVolcano("El Reventador", "el_reventador", "1502-01=", "Ecuador", -0.077, -77.656, 3562, -5, 3));
    m_volcanoes.push_back(CVolcano("Sangay", "sangay", "1502-09=", "Ecuador", -2.002, -78.341, 5230, -5, 3));
    m_volcanoes.push_back(CVolcano("Cayambe", "cayambe", "1502-006", "Ecuador", 0.029, -77.986, 5790, -5, 3));
    m_volcanoes.push_back(CVolcano("Guagua Pichincha", "guagua_pichincha", "1502-02=", "Ecuador", -0.171, -78.598, 4784, -5, 3));
    m_volcanoes.push_back(CVolcano("Merapi", "merapi", "0603-25=", "Indonesia", -7.542, 110.442, 2968, +11, 1));
    m_volcanoes.push_back(CVolcano("Kilauea", "kilauea", "1302-01-", "USA", 19.421, -155.287, 1222, -10, 18));
    m_volcanoes.push_back(CVolcano("Chalmers", "chalmers", "0000-000", "Sweden", 0.0, 0.0, 0, +1, 1));
    m_volcanoes.push_back(CVolcano("Harestua", "harestua", "0000-000", "Norway", 60.21, 10.75, 600, +1, 1));

    InitializeDatabase_01();
    InitializeDatabase_02();
    InitializeDatabase_03();
    InitializeDatabase_04();
    InitializeDatabase_05();
    InitializeDatabase_06();
    InitializeDatabase_07();
    InitializeDatabase_08();
    InitializeDatabase_09();
    InitializeDatabase_10();
    InitializeDatabase_11();
    InitializeDatabase_12();
    InitializeDatabase_13();
    InitializeDatabase_14();
    InitializeDatabase_15();
    InitializeDatabase_16();
    InitializeDatabase_17();
    InitializeDatabase_18();
    InitializeDatabase_19();


    m_volcanoNum = (unsigned int)m_volcanoes.size();
    m_preConfiguredVolcanoNum = m_volcanoNum;
}


void CVolcanoInfo::InitializeDatabase_01()
{
    m_volcanoes.push_back(CVolcano("West Eifel Volc Field", "", "0100-01-", "Germany", -50.17, 6.85, 600, 1, 0));
    m_volcanoes.push_back(CVolcano("Chaine des Puys", "", "0100-02-", "France", -45.775, 2.97, 1464, 1, 0));
    m_volcanoes.push_back(CVolcano("Olot Volc Field", "", "0100-03-", "Spain", -42.17, 2.53, 893, 1, 0));
    m_volcanoes.push_back(CVolcano("Calatrava Volc Field", "", "0100-04-", "Spain", -38.87, -4.02, 1117, 1, 0));
    m_volcanoes.push_back(CVolcano("Larderello", "", "0101-001", "Italy", -43.25, 10.87, 500, 1, 0));
    m_volcanoes.push_back(CVolcano("Vulsini", "", "0101-003", "Italy", -42.6, 11.93, 800, 1, 0));
    m_volcanoes.push_back(CVolcano("Alban Hills", "", "0101-004", "Italy", -41.73, 12.7, 949, 1, 0));
    m_volcanoes.push_back(CVolcano("Campi Flegrei", "", "0101-01=", "Italy", -40.827, 14.139, 458, 1, 0));
    m_volcanoes.push_back(CVolcano("Vesuvius", "", "0101-02=", "Italy", -40.821, 14.426, 1281, 1, 0));
    m_volcanoes.push_back(CVolcano("Ischia", "", "0101-03=", "Italy", -40.73, 13.897, 789, 1, 0));
    m_volcanoes.push_back(CVolcano("Panarea", "", "0101-041", "Italy", -38.63, 15.07, 421, 1, 0));
    m_volcanoes.push_back(CVolcano("Lipari", "", "0101-042", "Italy", -38.48, 14.95, 602, 1, 0));
    m_volcanoes.push_back(CVolcano("Pantelleria", "", "0101-071", "Italy", -36.77, 12.02, 836, 1, 0));
    m_volcanoes.push_back(CVolcano("Campi Flegrei Mar Sicilia", "", "0101-07=", "Italy", -37.1, 12.7, -8, 1, 0));
    m_volcanoes.push_back(CVolcano("Methana", "", "0102-02=", "Greece", -37.615, 23.336, 760, 1, 0));
    m_volcanoes.push_back(CVolcano("Mílos", "", "0102-03=", "Greece", -36.699, 24.439, 751, 1, 0));
    m_volcanoes.push_back(CVolcano("Santorini", "", "0102-04=", "Greece", -36.404, 25.396, 367, 1, 0));
    m_volcanoes.push_back(CVolcano("Yali", "", "0102-051", "Greece", -36.671, 27.14, 180, 1, 0));
    m_volcanoes.push_back(CVolcano("Nisyros", "", "0102-05=", "Greece", -36.586, 27.16, 698, 1, 0));
    m_volcanoes.push_back(CVolcano("Kula", "", "0103-00-", "Turkey", -38.58, 28.52, 750, 1, 0));
    m_volcanoes.push_back(CVolcano("Karapinar Field", "", "0103-001", "Turkey", -37.67, 33.65, 1302, 1, 0));
    m_volcanoes.push_back(CVolcano("Hasan Dagi", "", "0103-002", "Turkey", -38.13, 34.17, 3253, 1, 0));
    m_volcanoes.push_back(CVolcano("Göllü Dag", "", "0103-003", "Turkey", -38.25, 34.57, 2143, 1, 0));
    m_volcanoes.push_back(CVolcano("Acigöl-Nevsehir", "", "0103-004", "Turkey", -38.57, 34.52, 1689, 1, 0));
    m_volcanoes.push_back(CVolcano("Karaca Dag", "", "0103-011", "Turkey", -37.67, 39.83, 1957, 1, 0));
    m_volcanoes.push_back(CVolcano("Erciyes Dagi", "", "0103-01=", "Turkey", -38.52, 35.48, 3916, 1, 0));
    m_volcanoes.push_back(CVolcano("Süphan Dagi", "", "0103-021", "Turkey", -38.92, 42.82, 4158, 1, 0));
    m_volcanoes.push_back(CVolcano("Girekol", "", "0103-022", "Turkey", -39.17, 43.33, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Nemrut Dagi", "", "0103-02=", "Turkey", -38.65, 42.23, 2948, 1, 0));
    m_volcanoes.push_back(CVolcano("Tendürek Dagi", "", "0103-03=", "Turkey", -39.37, 43.87, 3584, 1, 0));
    m_volcanoes.push_back(CVolcano("Ararat", "", "0103-04-", "Turkey", -39.7, 44.3, 5165, 1, 0));
    m_volcanoes.push_back(CVolcano("Kars Plateau", "", "0103-05-", "Turkey", -40.75, 42.9, 3000, 1, 0));
    m_volcanoes.push_back(CVolcano("Elbrus", "", "0104-01-", "Russia", -43.33, 42.45, 5633, 1, 0));
    m_volcanoes.push_back(CVolcano("Kasbek", "", "0104-02-", "Georgia", -42.7, 44.5, 5050, 1, 0));
    m_volcanoes.push_back(CVolcano("Kabargin Oth Group", "", "0104-03-", "Georgia", -42.55, 44, 3650, 1, 0));
    m_volcanoes.push_back(CVolcano("Aragats", "", "0104-06-", "Armenia", -40.53, 44.2, 4095, 1, 0));
    m_volcanoes.push_back(CVolcano("Ghegam Ridge", "", "0104-07-", "Armenia", -40.275, 44.75, 3597, 1, 0));
    m_volcanoes.push_back(CVolcano("Dar-Alages", "", "0104-08-", "Armenia", -39.7, 45.542, 3329, 1, 0));
    m_volcanoes.push_back(CVolcano("Porak", "", "0104-09-", "Armenia", -40.02, 45.78, 2800, 1, 0));
    m_volcanoes.push_back(CVolcano("Tskhouk-Karckar", "", "0104-10-", "Armenia", -39.73, 46.02, 3000, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_02()
{
    m_volcanoes.push_back(CVolcano("Jebel at Tair", "", "0201-01=", "Red Sea", -15.55, 41.83, 244, 1, 0));
    m_volcanoes.push_back(CVolcano("Zukur", "", "0201-021", "Red Sea", -14.02, 42.75, 624, 1, 0));
    m_volcanoes.push_back(CVolcano("Hanish", "", "0201-022", "Red Sea", -13.72, 42.73, 422, 1, 0));
    m_volcanoes.push_back(CVolcano("Jebel Zubair", "", "0201-02=", "Red Sea", -15.05, 42.18, 191, 1, 0));
    m_volcanoes.push_back(CVolcano("Jalua", "", "0201-03=", "Ethiopia", -15.042, 39.82, 713, 1, 0));
    m_volcanoes.push_back(CVolcano("Dallol", "", "0201-041", "Ethiopia", -14.242, 40.3, -48, 1, 0));
    m_volcanoes.push_back(CVolcano("Alid", "", "0201-04=", "Ethiopia", -14.88, 39.92, 904, 1, 0));
    m_volcanoes.push_back(CVolcano("Gada Ale", "", "0201-05=", "Ethiopia", -13.975, 40.408, 287, 1, 0));
    m_volcanoes.push_back(CVolcano("Alu", "", "0201-06=", "Ethiopia", -13.825, 40.508, 429, 1, 0));
    m_volcanoes.push_back(CVolcano("Borale Ale", "", "0201-071", "Ethiopia", -13.725, 40.6, 668, 1, 0));
    m_volcanoes.push_back(CVolcano("Dalaffilla", "", "0201-07=", "Ethiopia", -13.792, 40.55, 613, 1, 0));
    m_volcanoes.push_back(CVolcano("Erta Ale", "", "0201-08=", "Ethiopia", -13.6, 40.67, 613, 1, 0));
    m_volcanoes.push_back(CVolcano("Hayli Gubbi", "", "0201-091", "Ethiopia", -13.5, 40.72, 521, 1, 0));
    m_volcanoes.push_back(CVolcano("Ale Bagu", "", "0201-09=", "Ethiopia", -13.52, 40.63, 1031, 1, 0));
    m_volcanoes.push_back(CVolcano("Nabro", "", "0201-101", "Ethiopia", -13.37, 41.7, 2218, 1, 0));
    m_volcanoes.push_back(CVolcano("Mallahle", "", "0201-102", "Ethiopia", -13.27, 41.65, 1875, 1, 0));
    m_volcanoes.push_back(CVolcano("Sork Ale", "", "0201-103", "Ethiopia", -13.18, 41.725, 1611, 1, 0));
    m_volcanoes.push_back(CVolcano("Asavyo", "", "0201-104", "Ethiopia", -13.07, 41.6, 1200, 1, 0));
    m_volcanoes.push_back(CVolcano("Mat Ala", "", "0201-105", "Ethiopia", -13.1, 41.15, 523, 1, 0));
    m_volcanoes.push_back(CVolcano("Tat Ali", "", "0201-106", "Ethiopia", -13.28, 41.07, 700, 1, 0));
    m_volcanoes.push_back(CVolcano("Borawli", "", "0201-107", "Ethiopia", -13.3, 40.98, 812, 1, 0));
    m_volcanoes.push_back(CVolcano("Dubbi", "", "0201-10=", "Ethiopia", -13.58, 41.808, 1625, 1, 0));
    m_volcanoes.push_back(CVolcano("Ma Alalta", "", "0201-111", "Ethiopia", -13.02, 40.2, 1815, 1, 0));
    m_volcanoes.push_back(CVolcano("Alayta", "", "0201-112", "Ethiopia", -12.88, 40.57, 1501, 1, 0));
    m_volcanoes.push_back(CVolcano("Dabbahu", "", "0201-113", "Ethiopia", -12.6, 40.48, 1442, 1, 0));
    m_volcanoes.push_back(CVolcano("Dabbayra", "", "0201-114", "Ethiopia", -12.38, 40.07, 1302, 1, 0));
    m_volcanoes.push_back(CVolcano("Manda Hararo", "", "0201-115", "Ethiopia", -12.17, 40.82, 600, 1, 0));
    m_volcanoes.push_back(CVolcano("Groppo", "", "0201-116", "Ethiopia", -11.73, 40.25, 930, 1, 0));
    m_volcanoes.push_back(CVolcano("Afderà", "", "0201-11=", "Ethiopia", -13.08, 40.85, 1295, 1, 0));
    m_volcanoes.push_back(CVolcano("Borawli", "", "0201-121", "Ethiopia", -11.63, 41.45, 875, 1, 0));
    m_volcanoes.push_back(CVolcano("Manda-Inakir", "", "0201-122", "Ethiopia", -12.38, 42.2, 600, 1, 0));
    m_volcanoes.push_back(CVolcano("Mousa Alli", "", "0201-123", "Ethiopia", -12.47, 42.4, 2028, 1, 0));
    m_volcanoes.push_back(CVolcano("Gufa", "", "0201-124", "Ethiopia", -12.55, 42.53, 600, 1, 0));
    m_volcanoes.push_back(CVolcano("Assab Volc Field", "", "0201-125", "Ethiopia", -12.95, 42.43, 987, 1, 0));
    m_volcanoes.push_back(CVolcano("Ardoukôba", "", "0201-126", "Djibouti", -11.58, 42.47, 298, 1, 0));
    m_volcanoes.push_back(CVolcano("Kurub", "", "0201-12=", "Ethiopia", -11.88, 41.208, 625, 1, 0));
    m_volcanoes.push_back(CVolcano("Dama Ali", "", "0201-141", "Ethiopia", -11.28, 41.63, 1068, 1, 0));
    m_volcanoes.push_back(CVolcano("Yangudi", "", "0201-151", "Ethiopia", -10.58, 41.042, 1383, 1, 0));
    m_volcanoes.push_back(CVolcano("Gabillema", "", "0201-15=", "Ethiopia", -11.08, 41.27, 1459, 1, 0));
    m_volcanoes.push_back(CVolcano("Ayelu", "", "0201-16=", "Ethiopia", -10.082, 40.702, 2145, 1, 0));
    m_volcanoes.push_back(CVolcano("Hertali", "", "0201-171", "Ethiopia", -9.78, 40.33, 900, 1, 0));
    m_volcanoes.push_back(CVolcano("Liado Hayk", "", "0201-172", "Ethiopia", -9.57, 40.28, 878, 1, 0));
    m_volcanoes.push_back(CVolcano("Adwa", "", "0201-17=", "Ethiopia", -10.07, 40.84, 1733, 1, 0));
    m_volcanoes.push_back(CVolcano("Dofen", "", "0201-18=", "Ethiopia", -9.35, 40.13, 1151, 1, 0));
    m_volcanoes.push_back(CVolcano("Beru", "", "0201-191", "Ethiopia", -8.95, 39.75, 1100, 1, 0));
    m_volcanoes.push_back(CVolcano("Fentale", "", "0201-19=", "Ethiopia", -8.975, 39.93, 2007, 1, 0));
    m_volcanoes.push_back(CVolcano("Kone", "", "0201-20-", "Ethiopia", -8.8, 39.692, 1619, 1, 0));
    m_volcanoes.push_back(CVolcano("Boset-Bericha", "", "0201-21-", "Ethiopia", -8.558, 39.475, 2447, 1, 0));
    m_volcanoes.push_back(CVolcano("Bishoftu Volc Field", "", "0201-22-", "Ethiopia", -8.78, 38.98, 1850, 1, 0));
    m_volcanoes.push_back(CVolcano("Sodore", "", "0201-222", "Ethiopia", -8.43, 39.35, 1765, 1, 0));
    m_volcanoes.push_back(CVolcano("Gedamsa Caldera", "", "0201-23-", "Ethiopia", -8.35, 39.18, 1984, 1, 0));
    m_volcanoes.push_back(CVolcano("Bora-Bericcio", "", "0201-24-", "Ethiopia", -8.27, 39.03, 2285, 1, 0));
    m_volcanoes.push_back(CVolcano("Tullu Moje", "", "0201-25-", "Ethiopia", -8.158, 39.13, 2349, 1, 0));
    m_volcanoes.push_back(CVolcano("East Zway", "", "0201-252", "Ethiopia", -7.95, 38.93, 1889, 1, 0));
    m_volcanoes.push_back(CVolcano("Butajiri-Silti Field", "", "0201-26-", "Ethiopia", -8.05, 38.35, 2281, 1, 0));
    m_volcanoes.push_back(CVolcano("Alutu", "", "0201-27-", "Ethiopia", -7.77, 38.78, 2335, 1, 0));
    m_volcanoes.push_back(CVolcano("O'a Caldera", "", "0201-28-", "Ethiopia", -7.47, 38.58, 2075, 1, 0));
    m_volcanoes.push_back(CVolcano("Corbetti Caldera", "", "0201-29-", "Ethiopia", -7.18, 38.43, 2320, 1, 0));
    m_volcanoes.push_back(CVolcano("Bilate River Field", "", "0201-291", "Ethiopia", -7.07, 38.1, 1700, 1, 0));
    m_volcanoes.push_back(CVolcano("Tepi", "", "0201-292", "Ethiopia", -7.42, 35.43, 2728, 1, 0));
    m_volcanoes.push_back(CVolcano("Hobicha Caldera", "", "0201-293", "Ethiopia", -6.78, 37.83, 1800, 1, 0));
    m_volcanoes.push_back(CVolcano("Chiracha", "", "0201-30-", "Ethiopia", -6.65, 38.12, 1650, 1, 0));
    m_volcanoes.push_back(CVolcano("Tosa Sucha", "", "0201-31-", "Ethiopia", -5.93, 37.57, 1650, 1, 0));
    m_volcanoes.push_back(CVolcano("Korath Range", "", "0201-32-", "Ethiopia", -5.1, 35.88, 912, 1, 0));
    m_volcanoes.push_back(CVolcano("Mega Basalt Field", "", "0201-33-", "Ethiopia", -4.08, 37.42, 1067, 1, 0));
    m_volcanoes.push_back(CVolcano("North Island", "", "0202-001", "Africa-E", -4.07, 36.05, 520, 1, 0));
    m_volcanoes.push_back(CVolcano("Central Island", "", "0202-01=", "Africa-E", -3.5, 36.042, 550, 1, 0));
    m_volcanoes.push_back(CVolcano("Marsabit", "", "0202-021", "Africa-E", -2.32, 37.97, 1707, 1, 0));
    m_volcanoes.push_back(CVolcano("South Island", "", "0202-02=", "Africa-E", -2.63, 36.6, 800, 1, 0));
    m_volcanoes.push_back(CVolcano("The Barrier", "", "0202-03=", "Africa-E", -2.32, 36.57, 1032, 1, 0));
    m_volcanoes.push_back(CVolcano("Namarunu", "", "0202-04-", "Africa-E", -1.98, 36.43, 817, 1, 0));
    m_volcanoes.push_back(CVolcano("Segererua Plateau", "", "0202-05-", "Africa-E", -1.57, 37.9, 699, 1, 0));
    m_volcanoes.push_back(CVolcano("Emuruangogolak", "", "0202-051", "Africa-E", -1.5, 36.33, 1328, 1, 0));
    m_volcanoes.push_back(CVolcano("Silali", "", "0202-052", "Africa-E", -1.15, 36.23, 1528, 1, 0));
    m_volcanoes.push_back(CVolcano("Paka", "", "0202-053", "Africa-E", -0.92, 36.18, 1697, 1, 0));
    m_volcanoes.push_back(CVolcano("Korosi", "", "0202-054", "Africa-E", -0.77, 36.12, 1446, 1, 0));
    m_volcanoes.push_back(CVolcano("Ol Kokwe", "", "0202-055", "Africa-E", -0.62, 36.075, 1130, 1, 0));
    m_volcanoes.push_back(CVolcano("Nyambeni Hills", "", "0202-056", "Africa-E", -0.23, 37.87, 750, 1, 0));
    m_volcanoes.push_back(CVolcano("Menengai", "", "0202-06=", "Africa-E", 0.2, 36.07, 2278, 1, 0));
    m_volcanoes.push_back(CVolcano("Elmenteita Badlands", "", "0202-071", "Africa-E", 0.52, 36.27, 2126, 1, 0));
    m_volcanoes.push_back(CVolcano("Homa Mountain", "", "0202-07=", "Africa-E", 0.38, 34.5, 1751, 1, 0));
    m_volcanoes.push_back(CVolcano("Ol Doinyo Eburru", "", "0202-08=", "Africa-E", 0.65, 36.22, 2856, 1, 0));
    m_volcanoes.push_back(CVolcano("Olkaria", "", "0202-09=", "Africa-E", 0.904, 36.292, 2434, 1, 0));
    m_volcanoes.push_back(CVolcano("Longonot", "", "0202-10=", "Africa-E", 0.914, 36.446, 2776, 1, 0));
    m_volcanoes.push_back(CVolcano("Suswa", "", "0202-11=", "Africa-E", 1.175, 36.35, 2356, 1, 0));
    m_volcanoes.push_back(CVolcano("Ol Doinyo Lengai", "", "0202-12=", "Africa-E", 2.764, 35.914, 2962, 1, 0));
    m_volcanoes.push_back(CVolcano("Chyulu Hills", "", "0202-13=", "Africa-E", 2.68, 37.88, 2188, 1, 0));
    m_volcanoes.push_back(CVolcano("Kilimanjaro", "", "0202-15=", "Africa-E", 3.07, 37.35, 5895, 1, 0));
    m_volcanoes.push_back(CVolcano("Igwisi Hills", "", "0202-161", "Africa-E", 4.87, 31.92, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("SW Usangu Basin", "", "0202-163", "Africa-E", 8.75, 33.8, 2179, 1, 0));
    m_volcanoes.push_back(CVolcano("Ngozi", "", "0202-164", "Africa-E", 8.97, 33.57, 2622, 1, 0));
    m_volcanoes.push_back(CVolcano("Izumbwe-Mpoli", "", "0202-165", "Africa-E", 8.93, 33.4, 1568, 1, 0));
    m_volcanoes.push_back(CVolcano("Rungwe", "", "0202-166", "Africa-E", 9.13, 33.67, 2961, 1, 0));
    m_volcanoes.push_back(CVolcano("Meru", "", "0202-16=", "Africa-E", 3.25, 36.75, 4565, 1, 0));
    m_volcanoes.push_back(CVolcano("Kieyo", "", "0202-17=", "Africa-E", 9.23, 33.78, 2175, 1, 0));
    m_volcanoes.push_back(CVolcano("Fort Portal", "", "0203-001", "Africa-C", -0.7, 30.25, 1615, 1, 0));
    m_volcanoes.push_back(CVolcano("Kyatwa", "", "0203-002", "Africa-C", -0.45, 30.25, 1430, 1, 0));
    m_volcanoes.push_back(CVolcano("Katwe-Kikorongo", "", "0203-003", "Africa-C", 0.08, 29.92, 1067, 1, 0));
    m_volcanoes.push_back(CVolcano("Bunyaruguru", "", "0203-004", "Africa-C", 0.2, 30.08, 1554, 1, 0));
    m_volcanoes.push_back(CVolcano("Katunga", "", "0203-005", "Africa-C", 0.471, 30.191, 1707, 1, 0));
    m_volcanoes.push_back(CVolcano("May-ya-moto", "", "0203-01=", "Africa-C", 0.93, 29.33, 950, 1, 0));
    m_volcanoes.push_back(CVolcano("Karisimbi", "", "0203-04-", "Africa-C", 1.5, 29.45, 4507, 1, 0));
    m_volcanoes.push_back(CVolcano("Visoke", "", "0203-05-", "Africa-C", 1.47, 29.492, 3711, 1, 0));
    m_volcanoes.push_back(CVolcano("Muhavura", "", "0203-06-", "Africa-C", 1.38, 29.67, 4127, 1, 0));
    m_volcanoes.push_back(CVolcano("Bufumbira", "", "0203-07-", "Africa-C", 1.23, 29.72, 2440, 1, 0));
    m_volcanoes.push_back(CVolcano("Tshibinda", "", "0203-08-", "Africa-C", 2.32, 28.75, 1460, 1, 0));
    m_volcanoes.push_back(CVolcano("Sao Tome", "", "0204-001", "Africa-W", -0.2, 6.58, 2024, 1, 0));
    m_volcanoes.push_back(CVolcano("San Carlos", "", "0204-002", "Africa-W", -3.35, 8.52, 2260, 1, 0));
    m_volcanoes.push_back(CVolcano("San Joaquin", "", "0204-003", "Africa-W", -3.35, 8.63, 2009, 1, 0));
    m_volcanoes.push_back(CVolcano("Santa Isabel", "", "0204-004", "Africa-W", -3.58, 8.75, 3007, 1, 0));
    m_volcanoes.push_back(CVolcano("Tombel Graben", "", "0204-011", "Africa-W", -4.75, 9.67, 500, 1, 0));
    m_volcanoes.push_back(CVolcano("Cameroon", "", "0204-01=", "Africa-W", -4.203, 9.17, 4095, 1, 0));
    m_volcanoes.push_back(CVolcano("Manengouba", "", "0204-02-", "Africa-W", -5.03, 9.83, 2411, 1, 0));
    m_volcanoes.push_back(CVolcano("Oku Volc Field", "", "0204-03-", "Africa-W", -6.25, 10.5, 3011, 1, 0));
    m_volcanoes.push_back(CVolcano("Ngaoundere Plateau", "", "0204-04-", "Africa-W", -7.25, 13.67, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Biu Plateau", "", "0204-05-", "Africa-W", -10.75, 12, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Todra Volc Field", "", "0205-001", "Africa-N", -17.68, 8.5, 1780, 1, 0));
    m_volcanoes.push_back(CVolcano("Tin Zaouatene Volc Field", "", "0205-002", "Africa-N", -19.83, 2.83, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("In Ezzane Volc Field", "", "0205-003", "Africa-N", -23, 10.83, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Tahalra Volc Field", "", "0205-004", "Africa-N", -22.67, 5, 1467, 1, 0));
    m_volcanoes.push_back(CVolcano("Atakor Volc Field", "", "0205-005", "Africa-N", -23.33, 5.83, 2918, 1, 0));
    m_volcanoes.push_back(CVolcano("Manzaz Volc Field", "", "0205-006", "Africa-N", -23.92, 5.83, 1672, 1, 0));
    m_volcanoes.push_back(CVolcano("Haruj", "", "0205-007", "Africa-N", -27.25, 17.5, 1200, 1, 0));
    m_volcanoes.push_back(CVolcano("Wau-en-Namus", "", "0205-008", "Africa-N", -25.05, 17.55, 547, 1, 0));
    m_volcanoes.push_back(CVolcano("Tarso Tôh", "", "0205-009", "Africa-N", -21.33, 16.33, 2000, 1, 0));
    m_volcanoes.push_back(CVolcano("Tarso Toussidé", "", "0205-01=", "Africa-N", -21.03, 16.45, 3265, 1, 0));
    m_volcanoes.push_back(CVolcano("Emi Koussi", "", "0205-021", "Africa-N", -19.8, 18.53, 3415, 1, 0));
    m_volcanoes.push_back(CVolcano("Tarso Voon", "", "0205-02=", "Africa-N", -20.92, 17.28, 3100, 1, 0));
    m_volcanoes.push_back(CVolcano("Jebel Marra", "", "0205-03-", "Africa-N", -12.95, 24.27, 3042, 1, 0));
    m_volcanoes.push_back(CVolcano("Kutum Volc Field", "", "0205-04-", "Africa-N", -14.57, 25.85, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Meidob Volc Field", "", "0205-05-", "Africa-N", -15.32, 26.47, 2000, 1, 0));
    m_volcanoes.push_back(CVolcano("Bayuda Volc Field", "", "0205-06-", "Africa-N", -18.33, 32.75, 670, 1, 0));
    m_volcanoes.push_back(CVolcano("Jebel Umm Arafieb", "", "0205-07-", "Africa-N", -18.17, 33.83, 0, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_03()
{
    m_volcanoes.push_back(CVolcano("Sharat Kovakab", "", "0300-01-", "Syria", -36.53, 40.85, 534, 1, 0));
    m_volcanoes.push_back(CVolcano("Golan Heights", "", "0300-03-", "Syria", -33.1, 35.97, 1197, 1, 0));
    m_volcanoes.push_back(CVolcano("Es Safa", "", "0300-05-", "Syria", -33.25, 37.07, 979, 1, 0));
    m_volcanoes.push_back(CVolcano("Jabal ad Druze", "", "0300-06-", "Syria", -32.658, 36.425, 1803, 1, 0));
    m_volcanoes.push_back(CVolcano("Al Harrah", "", "0301-001", "Arabia-W", -31.08, 38.42, 1100, 1, 0));
    m_volcanoes.push_back(CVolcano("Harrat ar Rahah", "", "0301-01=", "Arabia-W", -27.8, 36.17, 1950, 1, 0));
    m_volcanoes.push_back(CVolcano("Harrat Uwayrid", "", "0301-02=", "Arabia-W", -27.08, 37.25, 1920, 1, 0));
    m_volcanoes.push_back(CVolcano("Harrat Lunayyir", "", "0301-04-", "Arabia-W", -25.17, 37.75, 1370, 1, 0));
    m_volcanoes.push_back(CVolcano("Harrat Ithnayn", "", "0301-05=", "Arabia-W", -26.58, 40.2, 1625, 1, 0));
    m_volcanoes.push_back(CVolcano("Harrat Khaybar", "", "0301-06=", "Arabia-W", -25, 39.92, 2093, 1, 0));
    m_volcanoes.push_back(CVolcano("Harrat Kishb", "", "0301-071", "Arabia-W", -22.8, 41.38, 1475, 1, 0));
    m_volcanoes.push_back(CVolcano("Harrat al Birk", "", "0301-072", "Arabia-W", -18.37, 41.63, 381, 1, 0));
    m_volcanoes.push_back(CVolcano("Harrat Rahat", "", "0301-07=", "Arabia-W", -23.08, 39.78, 1744, 1, 0));
    m_volcanoes.push_back(CVolcano("Jabal Yar", "", "0301-08-", "Arabia-W", -17.05, 42.83, 305, 1, 0));
    m_volcanoes.push_back(CVolcano("Harra of Arhab", "", "0301-09-", "Arabia-S", -15.63, 44.08, 3100, 1, 0));
    m_volcanoes.push_back(CVolcano("Jabal el- Marha", "", "0301-10-", "Arabia-S", -15.245, 44.236, 2506, 1, 0));
    m_volcanoes.push_back(CVolcano("Jabal Haylan", "", "0301-11-", "Arabia-S", -15.43, 44.78, 1550, 1, 0));
    m_volcanoes.push_back(CVolcano("Harras of Dhamar", "", "0301-12-", "Arabia-S", -14.57, 44.67, 3500, 1, 0));
    m_volcanoes.push_back(CVolcano("Harra es- Sawâd", "", "0301-16-", "Arabia-S", -13.58, 46.12, 1737, 1, 0));
    m_volcanoes.push_back(CVolcano("Harra of Bal Haf", "", "0301-17-", "Arabia-S", -14.05, 48.33, 233, 1, 0));
    m_volcanoes.push_back(CVolcano("Bir Borhut", "", "0301-18-", "Arabia-S", -15.55, 50.63, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Sahand", "", "0302-001", "Iran", -37.75, 46.43, 3707, 1, 0));
    m_volcanoes.push_back(CVolcano("Sabalan", "", "0302-002", "Iran", -38.25, 47.92, 4811, 1, 0));
    m_volcanoes.push_back(CVolcano("Damavand", "", "0302-01-", "Iran", -35.951, 52.109, 5670, 1, 0));
    m_volcanoes.push_back(CVolcano("Qal'eh Hasan Ali", "", "0302-02-", "Iran", -29.4, 57.57, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Bazman", "", "0302-03-", "Iran", -28.07, 60, 3490, 1, 0));
    m_volcanoes.push_back(CVolcano("Taftan", "", "0302-05-", "Iran", -28.6, 61.13, 3940, 1, 0));
    m_volcanoes.push_back(CVolcano("Dacht-i-Navar Group", "", "0302-06-", "Afghanistan", -33.95, 67.92, 3800, 1, 0));
    m_volcanoes.push_back(CVolcano("Vakak Group", "", "0302-07-", "Afghanistan", -34.25, 67.97, 3190, 1, 0));
    m_volcanoes.push_back(CVolcano("La Grille", "", "0303-001", "Indian O.-W", 11.47, 43.33, 1087, 1, 0));
    m_volcanoes.push_back(CVolcano("Ambre-Bobaomby", "", "0303-011", "Madagascar", 12.6, 49.15, 1475, 1, 0));
    m_volcanoes.push_back(CVolcano("Nosy-Be", "", "0303-012", "Madagascar", 13.32, 48.48, 214, 1, 0));
    m_volcanoes.push_back(CVolcano("Ankaizina Field", "", "0303-013", "Madagascar", 14.3, 48.67, 2878, 1, 0));
    m_volcanoes.push_back(CVolcano("Itasy Volc Field", "", "0303-014", "Madagascar", 19, 46.77, 1800, 1, 0));
    m_volcanoes.push_back(CVolcano("Ankaratra Field", "", "0303-015", "Madagascar", 19.4, 47.2, 2644, 1, 0));
    m_volcanoes.push_back(CVolcano("Karthala", "", "0303-01=", "Indian O.-W", 11.75, 43.38, 2361, 1, 0));
    m_volcanoes.push_back(CVolcano("Boomerang Seamount", "", "0304-00-", "Indian O.-S", 37.721, 77.825, -650, 1, 0));
    m_volcanoes.push_back(CVolcano("Amsterdam Island", "", "0304-001", "Indian O.-S", 37.83, 77.52, 881, 1, 0));
    m_volcanoes.push_back(CVolcano("St. Paul", "", "0304-002", "Indian O.-S", 38.72, 77.53, 268, 1, 0));
    m_volcanoes.push_back(CVolcano("McDonald Islands", "", "0304-011", "Indian O.-S", 53.03, 72.6, 230, 1, 0));
    m_volcanoes.push_back(CVolcano("Heard", "", "0304-01=", "Indian O.-S", 53.106, 73.513, 2745, 1, 0));
    m_volcanoes.push_back(CVolcano("Kerguelen Islands", "", "0304-02=", "Indian O.-S", 49.58, 69.5, 1840, 1, 0));
    m_volcanoes.push_back(CVolcano("Ile de l'Est", "", "0304-03-", "Indian O.-S", 46.43, 52.2, 1090, 1, 0));
    m_volcanoes.push_back(CVolcano("Ile de la Possession", "", "0304-04-", "Indian O.-S", 46.42, 51.75, 934, 1, 0));
    m_volcanoes.push_back(CVolcano("Ile aux Cochons", "", "0304-05-", "Indian O.-S", 46.1, 50.23, 775, 1, 0));
    m_volcanoes.push_back(CVolcano("Prince Edward Island", "", "0304-06-", "Indian O.", 46.63, 37.95, 672, 1, 0));
    m_volcanoes.push_back(CVolcano("Marion Island", "", "0304-07-", "Indian O.-S", 46.9, 37.75, 1230, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_04()
{
    m_volcanoes.push_back(CVolcano("Whangarei", "", "0401-011", "New Zealand", 35.75, 174.27, 397, 1, 0));
    m_volcanoes.push_back(CVolcano("Kaikohe-Bay of Islands", "", "0401-01=", "New Zealand", 35.3, 173.9, 388, 1, 0));
    m_volcanoes.push_back(CVolcano("Mayor Island", "", "0401-021", "New Zealand", 37.28, 176.25, 355, 1, 0));
    m_volcanoes.push_back(CVolcano("Auckland Field", "", "0401-02=", "New Zealand", 36.9, 174.87, 260, 1, 0));
    m_volcanoes.push_back(CVolcano("Taranaki,Egmont", "", "0401-03=", "New Zealand", 39.3, 174.07, 2518, 1, 0));
    m_volcanoes.push_back(CVolcano("White Island", "", "0401-04=", "New Zealand", 37.52, 177.18, 321, 1, 0));
    m_volcanoes.push_back(CVolcano("Okataina", "", "0401-05=", "New Zealand", 38.12, 176.5, 1111, 1, 0));
    m_volcanoes.push_back(CVolcano("Reporoa", "", "0401-06-", "New Zealand", 38.42, 176.33, 592, 1, 0));
    m_volcanoes.push_back(CVolcano("Maroa", "", "0401-061", "New Zealand", 38.42, 176.08, 1156, 1, 0));
    m_volcanoes.push_back(CVolcano("Taupo", "", "0401-07=", "New Zealand", 38.82, 176, 760, 1, 0));
    m_volcanoes.push_back(CVolcano("Tongariro", "", "0401-08=", "New Zealand", 39.13, 175.642, 1978, 1, 0));
    m_volcanoes.push_back(CVolcano("Clark", "", "0401-101", "New Zealand", 36.446, 177.839, -860, 1, 0));
    m_volcanoes.push_back(CVolcano("Tangaroa", "", "0401-102", "New Zealand", 36.321, 178.028, 600, 1, 0));
    m_volcanoes.push_back(CVolcano("Ruapehu", "", "0401-10=", "New Zealand", 39.28, 175.57, 2797, 1, 0));
    m_volcanoes.push_back(CVolcano("Rumble V", "", "0401-11-", "New Zealand", 36.142, 178.196, 400, 1, 0));
    m_volcanoes.push_back(CVolcano("Rumble IV", "", "0401-12-", "New Zealand", 36.13, 178.05, 500, 1, 0));
    m_volcanoes.push_back(CVolcano("Rumble III", "", "0401-13-", "New Zealand", 35.745, 178.478, -220, 1, 0));
    m_volcanoes.push_back(CVolcano("Rumble II West", "", "0401-131", "New Zealand", 35.353, 178.527, 1200, 1, 0));
    m_volcanoes.push_back(CVolcano("Healy", "", "0401-14-", "New Zealand", 35.004, 178.973, 980, 1, 0));
    m_volcanoes.push_back(CVolcano("Brothers", "", "0401-15-", "New Zealand", 34.875, 179.075, -1350, 1, 0));
    m_volcanoes.push_back(CVolcano("Volcano W", "", "0402-001", "Kermadec Is", 31.85, -179.18, -900, 1, 0));
    m_volcanoes.push_back(CVolcano("Curtis Island", "", "0402-01=", "Kermadec Is", 30.542, -178.561, 137, 1, 0));
    m_volcanoes.push_back(CVolcano("Macauley Island", "", "0402-021", "Kermadec Is", 30.2, -178.47, 238, 1, 0));
    m_volcanoes.push_back(CVolcano("Giggenbach", "", "0402-022", "Kermadec Is", 30.036, -178.712, -65, 1, 0));
    m_volcanoes.push_back(CVolcano("Raoul Island", "", "0402-03=", "Kermadec Is", 29.27, -177.92, 516, 1, 0));
    m_volcanoes.push_back(CVolcano("Monowai Seamount", "", "0402-05-", "Kermadec Is", 25.887, -177.188, -132, 1, 0));
    m_volcanoes.push_back(CVolcano("Hunga Tonga-Hunga Ha'apai", "", "0403-04=", "Tonga-SW Pacific", 20.57, -175.38, 149, 1, 0));
    m_volcanoes.push_back(CVolcano("Falcon Island", "", "0403-05=", "Tonga-SW Pacific", 20.32, -175.42, -17, 1, 0));
    m_volcanoes.push_back(CVolcano("Kao", "", "0403-061", "Tonga-SW Pacific", 19.67, -175.03, 1030, 1, 0));
    m_volcanoes.push_back(CVolcano("Tofua", "", "0403-06=", "Tonga-SW Pacific", 19.75, -175.07, 515, 1, 0));
    m_volcanoes.push_back(CVolcano("Metis Shoal", "", "0403-07=", "Tonga-SW Pacific", 19.18, -174.87, 43, 1, 0));
    m_volcanoes.push_back(CVolcano("Home Reef", "", "0403-08=", "Tonga-SW Pacific", 18.992, -174.775, -2, 1, 0));
    m_volcanoes.push_back(CVolcano("Late", "", "0403-09=", "Tonga-SW Pacific", 18.806, -174.65, 540, 1, 0));
    m_volcanoes.push_back(CVolcano("Tafahi", "", "0403-101", "Tonga-SW Pacific", 15.85, -173.72, 560, 1, 0));
    m_volcanoes.push_back(CVolcano("Curacoa", "", "0403-102", "Tonga-SW Pacific", 15.62, -173.67, -33, 1, 0));
    m_volcanoes.push_back(CVolcano("Fonualei", "", "0403-10=", "Tonga-SW Pacific", 18.02, -174.325, 180, 1, 0));
    m_volcanoes.push_back(CVolcano("Niuafo'ou", "", "0403-11=", "Tonga-SW Pacific", 15.6, -175.63, 260, 1, 0));
    m_volcanoes.push_back(CVolcano("Vailulu'u", "", "0404-00-", "Samoa-SW Pacific", 14.215, -169.058, -592, 1, 0));
    m_volcanoes.push_back(CVolcano("Ta'u", "", "0404-001", "Samoa-SW Pacific", 14.23, -169.454, 931, 1, 0));
    m_volcanoes.push_back(CVolcano("Ofu-Olosega", "", "0404-01=", "Samoa-SW Pacific", 14.175, -169.618, 639, 1, 0));
    m_volcanoes.push_back(CVolcano("Tutuila", "", "0404-02-", "Samoa-SW Pacific", 14.295, -170.7, 653, 1, 0));
    m_volcanoes.push_back(CVolcano("Upolu", "", "0404-03-", "Samoa-SW Pacific", 13.935, -171.72, 1100, 1, 0));
    m_volcanoes.push_back(CVolcano("Savai'i", "", "0404-04=", "Samoa-SW Pacific", 13.612, -172.525, 1858, 1, 0));
    m_volcanoes.push_back(CVolcano("Wallis Islands", "", "0404-05-", "SW Pacific", 13.3, -176.17, 143, 1, 0));
    m_volcanoes.push_back(CVolcano("Taveuni", "", "0405-01-", "Fiji Is-SW Pacific", 16.82, -179.97, 1241, 1, 0));
    m_volcanoes.push_back(CVolcano("Koro", "", "0405-02-", "Fiji Is-SW Pacific", 17.32, 179.4, 522, 1, 0));
    m_volcanoes.push_back(CVolcano("Nabukelevu", "", "0405-03-", "Fiji Is-SW Pacific", 19.12, 177.98, 805, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_05()
{
    m_volcanoes.push_back(CVolcano("St. Andrew Strait", "", "0500-01=", "Admiralty Is-SW Pacific", 2.38, 147.35, 270, 1, 0));
    m_volcanoes.push_back(CVolcano("Baluan", "", "0500-02-", "Admiralty Is-SW Pacific", 2.57, 147.28, 254, 1, 0));
    m_volcanoes.push_back(CVolcano("Blup Blup", "", "0501-001", "New Guinea-NE of", 3.507, 144.605, 402, 1, 0));
    m_volcanoes.push_back(CVolcano("Kadovar", "", "0501-002", "New Guinea-NE of", 3.63, 144.631, 365, 1, 0));
    m_volcanoes.push_back(CVolcano("Boisa", "", "0501-011", "New Guinea-NE of", 3.994, 144.963, 240, 1, 0));
    m_volcanoes.push_back(CVolcano("Bam", "", "0501-01=", "New Guinea-NE of", 3.613, 144.818, 685, 1, 0));
    m_volcanoes.push_back(CVolcano("Manam", "", "0501-02=", "New Guinea-NE of", 4.08, 145.037, 1807, 1, 0));
    m_volcanoes.push_back(CVolcano("Karkar", "", "0501-03=", "New Guinea-NE of", 4.649, 145.964, 1839, 1, 0));
    m_volcanoes.push_back(CVolcano("Yomba", "", "0501-041", "New Guinea-NE of", 4.9, 146.75, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Long Island", "", "0501-05=", "New Guinea-NE of", 5.358, 147.12, 1280, 1, 0));
    m_volcanoes.push_back(CVolcano("Umboi", "", "0501-06=", "New Guinea-NE of", 5.589, 147.875, 1548, 1, 0));
    m_volcanoes.push_back(CVolcano("Ritter Island", "", "0501-07=", "New Guinea-NE of", 5.52, 148.121, 140, 1, 0));
    m_volcanoes.push_back(CVolcano("Sakar", "", "0501-08=", "New Guinea-NE of", 5.414, 148.094, 992, 1, 0));
    m_volcanoes.push_back(CVolcano("Langila", "", "0502-01=", "New Britain-SW Pac", 5.525, 148.42, 1330, 1, 0));
    m_volcanoes.push_back(CVolcano("Mundua", "", "0502-021", "New Britain-SW Pac", 4.63, 149.35, 179, 1, 0));
    m_volcanoes.push_back(CVolcano("Garove", "", "0502-03=", "New Britain-SW Pac", 4.692, 149.5, 368, 1, 0));
    m_volcanoes.push_back(CVolcano("Dakataua", "", "0502-04=", "New Britain-SW Pac", 5.056, 150.108, 400, 1, 0));
    m_volcanoes.push_back(CVolcano("Bola", "", "0502-05=", "New Britain-SW Pac", 5.15, 150.03, 1155, 1, 0));
    m_volcanoes.push_back(CVolcano("Garua Harbour", "", "0502-06=", "New Britain-SW Pac", 5.3, 150.07, 565, 1, 0));
    m_volcanoes.push_back(CVolcano("Lolo", "", "0502-071", "New Britain-SW Pac", 5.468, 150.507, 805, 1, 0));
    m_volcanoes.push_back(CVolcano("Garbuna Group", "", "0502-07=", "New Britain-SW Pac", 5.45, 150.03, 564, 1, 0));
    m_volcanoes.push_back(CVolcano("Pago", "", "0502-08=", "New Britain-SW Pac", 5.58, 150.52, 742, 1, 0));
    m_volcanoes.push_back(CVolcano("Sulu Range", "", "0502-09=", "New Britain-SW Pac", 5.5, 150.942, 610, 1, 0));
    m_volcanoes.push_back(CVolcano("Hargy", "", "0502-10=", "New Britain-SW Pac", 5.33, 151.1, 1148, 1, 0));
    m_volcanoes.push_back(CVolcano("Bamus", "", "0502-11=", "New Britain-SW Pac", 5.2, 151.23, 2248, 1, 0));
    m_volcanoes.push_back(CVolcano("Ulawun", "", "0502-12=", "New Britain-SW Pac", 5.05, 151.33, 2334, 1, 0));
    m_volcanoes.push_back(CVolcano("Lolobau", "", "0502-13=", "New Britain-SW Pac", 4.92, 151.158, 858, 1, 0));
    m_volcanoes.push_back(CVolcano("Rabaul", "", "0502-14=", "New Britain-SW Pac", 4.271, 152.203, 688, 1, 0));
    m_volcanoes.push_back(CVolcano("Tavui", "", "0502-15-", "New Britain-SW Pac", 4.12, 152.2, 200, 1, 0));
    m_volcanoes.push_back(CVolcano("Doma Peaks", "", "0503-00-", "New Guinea", 5.9, 143.15, 3568, 1, 0));
    m_volcanoes.push_back(CVolcano("Crater Mountain", "", "0503-001", "New Guinea", 6.58, 145.08, 3233, 1, 0));
    m_volcanoes.push_back(CVolcano("Yelia", "", "0503-002", "New Guinea", 7.05, 145.858, 3384, 1, 0));
    m_volcanoes.push_back(CVolcano("Koranga", "", "0503-003", "New Guinea", 7.33, 146.708, 1500, 1, 0));
    m_volcanoes.push_back(CVolcano("Madilogo", "", "0503-004", "New Guinea", 9.2, 147.57, 850, 1, 0));
    m_volcanoes.push_back(CVolcano("Hydrographers Range", "", "0503-011", "New Guinea", 9, 148.37, 1915, 1, 0));
    m_volcanoes.push_back(CVolcano("Lamington", "", "0503-01=", "New Guinea", 8.95, 148.15, 1680, 1, 0));
    m_volcanoes.push_back(CVolcano("Managlase Plateau", "", "0503-021", "New Guinea", 9.08, 148.33, 1342, 1, 0));
    m_volcanoes.push_back(CVolcano("Musa River", "", "0503-02=", "New Guinea", 9.308, 148.13, 808, 1, 0));
    m_volcanoes.push_back(CVolcano("Sessagara", "", "0503-031", "New Guinea", 9.48, 149.13, 370, 1, 0));
    m_volcanoes.push_back(CVolcano("Victory", "", "0503-03=", "New Guinea", 9.2, 149.07, 1925, 1, 0));
    m_volcanoes.push_back(CVolcano("Goodenough", "", "0503-041", "D'Entrecasteaux Is", 9.48, 150.35, 220, 1, 0));
    m_volcanoes.push_back(CVolcano("Waiowa", "", "0503-04=", "New Guinea", 9.57, 149.075, 640, 1, 0));
    m_volcanoes.push_back(CVolcano("Iamalele", "", "0503-05=", "D'Entrecasteaux Is", 9.52, 150.53, 200, 1, 0));
    m_volcanoes.push_back(CVolcano("Dawson Strait Group", "", "0503-06=", "D'Entrecasteaux Is", 9.62, 150.88, 500, 1, 0));
    m_volcanoes.push_back(CVolcano("Lihir", "", "0504-01=", "New Ireland-SW Pacific", 3.125, 152.642, 700, 1, 0));
    m_volcanoes.push_back(CVolcano("Ambitle", "", "0504-02=", "New Ireland-SW Pacific", 4.08, 153.65, 450, 1, 0));
    m_volcanoes.push_back(CVolcano("Tore", "", "0505-00-", "Bougainville-SW Pacific", 5.83, 154.93, 2200, 1, 0));
    m_volcanoes.push_back(CVolcano("Billy Mitchell", "", "0505-011", "Bougainville-SW Pacific", 6.092, 155.225, 1544, 1, 0));
    m_volcanoes.push_back(CVolcano("Balbi", "", "0505-01=", "Bougainville-SW Pacific", 5.92, 154.98, 2715, 1, 0));
    m_volcanoes.push_back(CVolcano("Takuan Group", "", "0505-021", "Bougainville-SW Pacific", 6.442, 155.608, 2210, 1, 0));
    m_volcanoes.push_back(CVolcano("Bagana", "", "0505-02=", "Bougainville-SW Pacific", 6.14, 155.195, 1750, 1, 0));
    m_volcanoes.push_back(CVolcano("Loloru", "", "0505-03=", "Bougainville-SW Pacific", 6.52, 155.62, 1887, 1, 0));
    m_volcanoes.push_back(CVolcano("Kana Keoki", "", "0505-052", "Solomon Is-SW Pacific", 8.75, 157.03, -700, 1, 0));
    m_volcanoes.push_back(CVolcano("Coleman Seamount", "", "0505-053", "Solomon Is-SW Pacific", 8.83, 157.17, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Simbo", "", "0505-05=", "Solomon Is-SW Pacific", 8.292, 156.52, 335, 1, 0));
    m_volcanoes.push_back(CVolcano("Gallego", "", "0505-062", "Solomon Is-SW Pacific", 9.35, 159.73, 1000, 1, 0));
    m_volcanoes.push_back(CVolcano("Kavachi", "", "0505-06=", "Solomon Is-SW Pacific", 9.02, 157.95, -20, 1, 0));
    m_volcanoes.push_back(CVolcano("Savo", "", "0505-07=", "Solomon Is-SW Pacific", 9.13, 159.82, 485, 1, 0));
    m_volcanoes.push_back(CVolcano("Tinakula", "", "0506-01=", "Santa Cruz Is-SW Pacific", 10.38, 165.8, 851, 1, 0));
    m_volcanoes.push_back(CVolcano("Motlav", "", "0507-001", "Vanuatu-SW Pacific", 13.67, 167.67, 411, 1, 0));
    m_volcanoes.push_back(CVolcano("Suretamatai", "", "0507-01=", "Vanuatu-SW Pacific", 13.8, 167.47, 921, 1, 0));
    m_volcanoes.push_back(CVolcano("Mere Lava", "", "0507-021", "Vanuatu-SW Pacific", 14.45, 168.05, 1028, 1, 0));
    m_volcanoes.push_back(CVolcano("Gaua", "", "0507-02=", "Vanuatu-SW Pacific", 14.27, 167.5, 797, 1, 0));
    m_volcanoes.push_back(CVolcano("Aoba", "", "0507-03=", "Vanuatu-SW Pacific", 15.4, 167.83, 1496, 1, 0));
    m_volcanoes.push_back(CVolcano("Ambrym", "", "0507-04=", "Vanuatu-SW Pacific", 16.25, 168.12, 1334, 1, 0));
    m_volcanoes.push_back(CVolcano("Lopevi", "", "0507-05=", "Vanuatu-SW Pacific", 16.507, 168.346, 1413, 1, 0));
    m_volcanoes.push_back(CVolcano("Epi", "", "0507-06=", "Vanuatu-SW Pacific", 16.68, 168.37, 833, 1, 0));
    m_volcanoes.push_back(CVolcano("Kuwae", "", "0507-07=", "Vanuatu-SW Pacific", 16.829, 168.536, -2, 1, 0));
    m_volcanoes.push_back(CVolcano("North Vate", "", "0507-081", "Vanuatu-SW Pacific", 17.47, 168.353, 594, 1, 0));
    m_volcanoes.push_back(CVolcano("Traitor's Head", "", "0507-09=", "Vanuatu-SW Pacific", 18.75, 169.23, 837, 1, 0));
    m_volcanoes.push_back(CVolcano("Aneityum", "", "0507-11-", "Vanuatu-SW Pacific", 20.2, 169.78, 852, 1, 0));
    m_volcanoes.push_back(CVolcano("Eastern Gemini Seamount", "", "0508-001", "SW Pacific", 20.98, 170.28, -80, 1, 0));
    m_volcanoes.push_back(CVolcano("Matthew Island", "", "0508-01=", "SW Pacific", 22.33, 171.32, 177, 1, 0));
    m_volcanoes.push_back(CVolcano("Hunter Island", "", "0508-02=", "SW Pacific", 22.4, 172.05, 297, 1, 0));
    m_volcanoes.push_back(CVolcano("Newer Volcanics Prov", "", "0509-01-", "Australia", 37.77, 142.5, 1011, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_06()
{
    m_volcanoes.push_back(CVolcano("Narcondum", "", "0600-001", "Andaman Is-Indian O", -13.43, 94.28, 710, 1, 0));
    m_volcanoes.push_back(CVolcano("Barren Island", "", "0600-01=", "Andaman Is-Indian O", -12.278, 93.858, 354, 1, 0));
    m_volcanoes.push_back(CVolcano("Seulawah Agam", "", "0601-02=", "Sumatra", -5.448, 95.658, 1810, 1, 0));
    m_volcanoes.push_back(CVolcano("Peuet Sague", "", "0601-03=", "Sumatra", -4.914, 96.329, 2801, 1, 0));
    m_volcanoes.push_back(CVolcano("Bur ni Telong", "", "0601-05=", "Sumatra", -4.769, 96.821, 2617, 1, 0));
    m_volcanoes.push_back(CVolcano("Sibayak", "", "0601-07=", "Sumatra", -3.23, 98.52, 2212, 1, 0));
    m_volcanoes.push_back(CVolcano("Sinabung", "", "0601-08=", "Sumatra", -3.17, 98.392, 2460, 1, 0));
    m_volcanoes.push_back(CVolcano("Toba", "", "0601-09=", "Sumatra", -2.58, 98.83, 2157, 1, 0));
    m_volcanoes.push_back(CVolcano("Imun", "", "0601-101", "Sumatra", -2.158, 98.93, 1505, 1, 0));
    m_volcanoes.push_back(CVolcano("Lubukraya", "", "0601-111", "Sumatra", -1.478, 99.209, 1862, 1, 0));
    m_volcanoes.push_back(CVolcano("Sibualbuali", "", "0601-11=", "Sumatra", -1.556, 99.255, 1819, 1, 0));
    m_volcanoes.push_back(CVolcano("Sorikmarapi", "", "0601-12=", "Sumatra", -0.686, 99.539, 2145, 1, 0));
    m_volcanoes.push_back(CVolcano("Sarik-Gajah", "", "0601-131", "Sumatra", -0.08, 100.2, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Talakmau", "", "0601-13=", "Sumatra", -0.079, 99.98, 2919, 1, 0));
    m_volcanoes.push_back(CVolcano("Marapi", "", "0601-14=", "Sumatra", 0.381, 100.473, 2891, 1, 0));
    m_volcanoes.push_back(CVolcano("Tandikat", "", "0601-15=", "Sumatra", 0.433, 100.317, 2438, 1, 0));
    m_volcanoes.push_back(CVolcano("Talang", "", "0601-16=", "Sumatra", 0.978, 100.679, 2597, 1, 0));
    m_volcanoes.push_back(CVolcano("Kunyit", "", "0601-171", "Sumatra", 2.274, 101.483, 2151, 1, 0));
    m_volcanoes.push_back(CVolcano("Hutapanjang", "", "0601-172", "Sumatra", 2.33, 101.6, 2021, 1, 0));
    m_volcanoes.push_back(CVolcano("Kerinci", "", "0601-17=", "Sumatra", 1.697, 101.264, 3800, 1, 0));
    m_volcanoes.push_back(CVolcano("Sumbing", "", "0601-18=", "Sumatra", 2.414, 101.728, 2507, 1, 0));
    m_volcanoes.push_back(CVolcano("Pendan", "", "0601-191", "Sumatra", 2.82, 102.02, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Belirang-Beriti", "", "0601-20=", "Sumatra", 2.82, 102.18, 1958, 1, 0));
    m_volcanoes.push_back(CVolcano("Bukit Daun", "", "0601-21=", "Sumatra", 3.38, 102.37, 2467, 1, 0));
    m_volcanoes.push_back(CVolcano("Kaba", "", "0601-22=", "Sumatra", 3.52, 102.62, 1952, 1, 0));
    m_volcanoes.push_back(CVolcano("Patah", "", "0601-231", "Sumatra", 4.27, 103.3, 2817, 1, 0));
    m_volcanoes.push_back(CVolcano("Dempo", "", "0601-23=", "Sumatra", 4.03, 103.13, 3173, 1, 0));
    m_volcanoes.push_back(CVolcano("Bukit Lumut Balai", "", "0601-24=", "Sumatra", 4.22, 103.62, 2055, 1, 0));
    m_volcanoes.push_back(CVolcano("Ranau", "", "0601-251", "Sumatra", 4.83, 103.92, 1881, 1, 0));
    m_volcanoes.push_back(CVolcano("Besar", "", "0601-25=", "Sumatra", 4.43, 103.67, 1899, 1, 0));
    m_volcanoes.push_back(CVolcano("Sekincau Belirang", "", "0601-26=", "Sumatra", 5.12, 104.32, 1719, 1, 0));
    m_volcanoes.push_back(CVolcano("Suoh", "", "0601-27=", "Sumatra", 5.25, 104.27, 1000, 1, 0));
    m_volcanoes.push_back(CVolcano("Hulubelu", "", "0601-28=", "Sumatra", 5.35, 104.6, 1040, 1, 0));
    m_volcanoes.push_back(CVolcano("Rajabasa", "", "0601-29=", "Sumatra", 5.78, 105.625, 1281, 1, 0));
    m_volcanoes.push_back(CVolcano("Krakatau", "", "0602-00=", "Indonesia", 6.102, 105.423, 813, 1, 0));
    m_volcanoes.push_back(CVolcano("Pulosari", "", "0603-01=", "Java", 6.342, 105.975, 1346, 1, 0));
    m_volcanoes.push_back(CVolcano("Karang", "", "0603-02=", "Java", 6.27, 106.042, 1778, 1, 0));
    m_volcanoes.push_back(CVolcano("Perbakti-Gagak", "", "0603-04=", "Java", 6.75, 106.7, 1699, 1, 0));
    m_volcanoes.push_back(CVolcano("Salak", "", "0603-05=", "Java", 6.72, 106.73, 2211, 1, 0));
    m_volcanoes.push_back(CVolcano("Gede", "", "0603-06=", "Java", 6.78, 106.98, 2958, 1, 0));
    m_volcanoes.push_back(CVolcano("Patuha", "", "0603-07=", "Java", 7.16, 107.4, 2434, 1, 0));
    m_volcanoes.push_back(CVolcano("Malabar", "", "0603-081", "Java", 7.13, 107.65, 2343, 1, 0));
    m_volcanoes.push_back(CVolcano("Wayang-Windu", "", "0603-08=", "Java", 7.208, 107.63, 2182, 1, 0));
    m_volcanoes.push_back(CVolcano("Tangkubanparahu", "", "0603-09=", "Java", 6.77, 107.6, 2084, 1, 0));
    m_volcanoes.push_back(CVolcano("Papandayan", "", "0603-10=", "Java", 7.32, 107.73, 2665, 1, 0));
    m_volcanoes.push_back(CVolcano("Kendang", "", "0603-11=", "Java", 7.23, 107.72, 2608, 1, 0));
    m_volcanoes.push_back(CVolcano("Tampomas", "", "0603-131", "Java", 6.77, 107.95, 1684, 1, 0));
    m_volcanoes.push_back(CVolcano("Guntur", "", "0603-13=", "Java", 7.143, 107.84, 2249, 1, 0));
    m_volcanoes.push_back(CVolcano("Galunggung", "", "0603-14=", "Java", 7.25, 108.058, 2168, 1, 0));
    m_volcanoes.push_back(CVolcano("Talagabodas", "", "0603-15=", "Java", 7.208, 108.07, 2201, 1, 0));
    m_volcanoes.push_back(CVolcano("Kawah Karaha", "", "0603-16=", "Java", 7.12, 108.08, 1155, 1, 0));
    m_volcanoes.push_back(CVolcano("Cereme", "", "0603-17=", "Java", 6.892, 108.4, 3078, 1, 0));
    m_volcanoes.push_back(CVolcano("Slamet", "", "0603-18=", "Java", 7.242, 109.208, 3428, 1, 0));
    m_volcanoes.push_back(CVolcano("Dieng Volc Complex", "", "0603-20=", "Java", 7.2, 109.92, 2565, 1, 0));
    m_volcanoes.push_back(CVolcano("Sundoro", "", "0603-21=", "Java", 7.3, 109.992, 3136, 1, 0));
    m_volcanoes.push_back(CVolcano("Sumbing", "", "0603-22=", "Java", 7.384, 110.07, 3371, 1, 0));
    m_volcanoes.push_back(CVolcano("Telomoyo", "", "0603-231", "Java", 7.37, 110.4, 1894, 1, 0));
    m_volcanoes.push_back(CVolcano("Ungaran", "", "0603-23=", "Java", 7.18, 110.33, 2050, 1, 0));
    m_volcanoes.push_back(CVolcano("Merbabu", "", "0603-24=", "Java", 7.45, 110.43, 3145, 1, 0));
    m_volcanoes.push_back(CVolcano("Muria", "", "0603-251", "Java", 6.62, 110.88, 1625, 1, 0));
    m_volcanoes.push_back(CVolcano("Lawu", "", "0603-26=", "Java", 7.625, 111.192, 3265, 1, 0));
    m_volcanoes.push_back(CVolcano("Wilis", "", "0603-27=", "Java", 7.808, 111.758, 2563, 1, 0));
    m_volcanoes.push_back(CVolcano("Kawi-Butak", "", "0603-281", "Java", 7.92, 112.45, 2651, 1, 0));
    m_volcanoes.push_back(CVolcano("Kelut", "", "0603-28=", "Java", 7.93, 112.308, 1731, 1, 0));
    m_volcanoes.push_back(CVolcano("Penanggungan", "", "0603-291", "Java", 7.62, 112.63, 1653, 1, 0));
    m_volcanoes.push_back(CVolcano("Malang Plain", "", "0603-292", "Java", 8.02, 112.68, 680, 1, 0));
    m_volcanoes.push_back(CVolcano("Arjuno-Welirang", "", "0603-29=", "Java", 7.725, 112.58, 3339, 1, 0));
    m_volcanoes.push_back(CVolcano("Semeru", "", "0603-30=", "Java", 8.108, 112.92, 3676, 1, 0));
    m_volcanoes.push_back(CVolcano("Tengger Caldera", "", "0603-31=", "Java", 7.942, 112.95, 2329, 1, 0));
    m_volcanoes.push_back(CVolcano("Lurus", "", "0603-321", "Java", 7.73, 113.58, 539, 1, 0));
    m_volcanoes.push_back(CVolcano("Lamongan", "", "0603-32=", "Java", 7.979, 113.342, 1651, 1, 0));
    m_volcanoes.push_back(CVolcano("Iyang-Argapura", "", "0603-33=", "Java", 7.97, 113.57, 3088, 1, 0));
    m_volcanoes.push_back(CVolcano("Raung", "", "0603-34=", "Java", 8.125, 114.042, 3332, 1, 0));
    m_volcanoes.push_back(CVolcano("Baluran", "", "0603-351", "Java", 7.85, 114.37, 1247, 1, 0));
    m_volcanoes.push_back(CVolcano("Ijen", "", "0603-35=", "Java", 8.058, 114.242, 2799, 1, 0));
    m_volcanoes.push_back(CVolcano("Bratan", "", "0604-001", "Lesser Sunda Is", 8.28, 115.13, 2276, 1, 0));
    m_volcanoes.push_back(CVolcano("Batur", "", "0604-01=", "Lesser Sunda Is", 8.242, 115.375, 1717, 1, 0));
    m_volcanoes.push_back(CVolcano("Agung", "", "0604-02=", "Lesser Sunda Is", 8.342, 115.508, 3142, 1, 0));
    m_volcanoes.push_back(CVolcano("Rinjani", "", "0604-03=", "Lesser Sunda Is", 8.42, 116.47, 3726, 1, 0));
    m_volcanoes.push_back(CVolcano("Tambora", "", "0604-04=", "Lesser Sunda Is", 8.25, 118, 2850, 1, 0));
    m_volcanoes.push_back(CVolcano("Sangeang Api", "", "0604-05=", "Lesser Sunda Is", 8.2, 119.07, 1949, 1, 0));
    m_volcanoes.push_back(CVolcano("Wai Sano", "", "0604-06=", "Lesser Sunda Is", 8.72, 120.02, 903, 1, 0));
    m_volcanoes.push_back(CVolcano("Ranakah", "", "0604-071", "Lesser Sunda Is", 8.62, 120.52, 2350, 1, 0));
    m_volcanoes.push_back(CVolcano("Poco Leok", "", "0604-07=", "Lesser Sunda Is", 8.68, 120.48, 1675, 1, 0));
    m_volcanoes.push_back(CVolcano("Inierie", "", "0604-08=", "Lesser Sunda Is", 8.875, 120.95, 2245, 1, 0));
    m_volcanoes.push_back(CVolcano("Inielika", "", "0604-09=", "Lesser Sunda Is", 8.73, 120.98, 1559, 1, 0));
    m_volcanoes.push_back(CVolcano("Ebulobo", "", "0604-10=", "Lesser Sunda Is", 8.82, 121.18, 2124, 1, 0));
    m_volcanoes.push_back(CVolcano("Iya", "", "0604-11=", "Lesser Sunda Is", 8.897, 121.645, 637, 1, 0));
    m_volcanoes.push_back(CVolcano("Sukaria Caldera", "", "0604-12=", "Lesser Sunda Is", 8.792, 121.77, 1500, 1, 0));
    m_volcanoes.push_back(CVolcano("Ndete Napu", "", "0604-13=", "Lesser Sunda Is", 8.72, 121.78, 750, 1, 0));
    m_volcanoes.push_back(CVolcano("Kelimutu", "", "0604-14=", "Lesser Sunda Is", 8.77, 121.82, 1639, 1, 0));
    m_volcanoes.push_back(CVolcano("Paluweh", "", "0604-15=", "Lesser Sunda Is", 8.32, 121.708, 875, 1, 0));
    m_volcanoes.push_back(CVolcano("Egon", "", "0604-16=", "Lesser Sunda Is", 8.67, 122.45, 1703, 1, 0));
    m_volcanoes.push_back(CVolcano("Ilimuda", "", "0604-17=", "Lesser Sunda Is", 8.478, 122.671, 1100, 1, 0));
    m_volcanoes.push_back(CVolcano("Lewotobi", "", "0604-18=", "Lesser Sunda Is", 8.542, 122.775, 1703, 1, 0));
    m_volcanoes.push_back(CVolcano("Leroboleng", "", "0604-20=", "Lesser Sunda Is", 8.358, 122.842, 1117, 1, 0));
    m_volcanoes.push_back(CVolcano("Iliboleng", "", "0604-22=", "Lesser Sunda Is", 8.342, 123.258, 1659, 1, 0));
    m_volcanoes.push_back(CVolcano("Lewotolo", "", "0604-23=", "Lesser Sunda Is", 8.272, 123.505, 1423, 1, 0));
    m_volcanoes.push_back(CVolcano("Ililabalekan", "", "0604-24=", "Lesser Sunda Is", 8.55, 123.38, 1018, 1, 0));
    m_volcanoes.push_back(CVolcano("Iliwerung", "", "0604-25=", "Lesser Sunda Is", 8.53, 123.57, 1018, 1, 0));
    m_volcanoes.push_back(CVolcano("Batu Tara", "", "0604-26=", "Lesser Sunda Is", 7.792, 123.579, 748, 1, 0));
    m_volcanoes.push_back(CVolcano("Sirung", "", "0604-27=", "Lesser Sunda Is", 8.508, 124.13, 862, 1, 0));
    m_volcanoes.push_back(CVolcano("Yersey", "", "0604-28=", "Lesser Sunda Is", 7.53, 123.95, -3800, 1, 0));
    m_volcanoes.push_back(CVolcano("Emperor of China", "", "0605-01=", "Banda Sea", 6.62, 124.22, -2850, 1, 0));
    m_volcanoes.push_back(CVolcano("Nieuwerkerk", "", "0605-02=", "Banda Sea", 6.6, 124.675, -2285, 1, 0));
    m_volcanoes.push_back(CVolcano("Gunungapi Wetar", "", "0605-03=", "Banda Sea", 6.642, 126.65, 282, 1, 0));
    m_volcanoes.push_back(CVolcano("Wurlali", "", "0605-04=", "Banda Sea", 7.125, 128.675, 868, 1, 0));
    m_volcanoes.push_back(CVolcano("Teon", "", "0605-05=", "Banda Sea", 6.92, 129.125, 655, 1, 0));
    m_volcanoes.push_back(CVolcano("Nila", "", "0605-06=", "Banda Sea", 6.73, 129.5, 781, 1, 0));
    m_volcanoes.push_back(CVolcano("Serua", "", "0605-07=", "Banda Sea", 6.3, 130, 641, 1, 0));
    m_volcanoes.push_back(CVolcano("Manuk", "", "0605-08=", "Banda Sea", 5.53, 130.292, 282, 1, 0));
    m_volcanoes.push_back(CVolcano("Banda Api", "", "0605-09=", "Banda Sea", 4.525, 129.871, 640, 1, 0));
    m_volcanoes.push_back(CVolcano("Colo [Una Una]", "", "0606-01=", "Sulawesi-Indonesia", 0.17, 121.608, 507, 1, 0));
    m_volcanoes.push_back(CVolcano("Ambang", "", "0606-02=", "Sulawesi-Indonesia", -0.75, 124.42, 1795, 1, 0));
    m_volcanoes.push_back(CVolcano("Soputan", "", "0606-03=", "Sulawesi-Indonesia", -1.108, 124.73, 1784, 1, 0));
    m_volcanoes.push_back(CVolcano("Sempu", "", "0606-04=", "Sulawesi-Indonesia", -1.13, 124.758, 1549, 1, 0));
    m_volcanoes.push_back(CVolcano("Tondano Caldera", "", "0606-07-", "Sulawesi-Indonesia", -1.23, 124.83, 1202, 1, 0));
    m_volcanoes.push_back(CVolcano("Lokon-Empung", "", "0606-10=", "Sulawesi-Indonesia", -1.358, 124.792, 1580, 1, 0));
    m_volcanoes.push_back(CVolcano("Mahawu", "", "0606-11=", "Sulawesi-Indonesia", -1.358, 124.858, 1324, 1, 0));
    m_volcanoes.push_back(CVolcano("Klabat", "", "0606-12=", "Sulawesi-Indonesia", -1.47, 125.03, 1995, 1, 0));
    m_volcanoes.push_back(CVolcano("Tongkoko", "", "0606-13=", "Sulawesi-Indonesia", -1.52, 125.2, 1149, 1, 0));
    m_volcanoes.push_back(CVolcano("Ruang", "", "0607-01=", "Sangihe Is-Indonesia", -2.3, 125.37, 725, 1, 0));
    m_volcanoes.push_back(CVolcano("Karangetang [Api Siau]", "", "0607-02=", "Sangihe Is-Indonesia", -2.78, 125.4, 1784, 1, 0));
    m_volcanoes.push_back(CVolcano("Banua Wuhu", "", "0607-03=", "Sangihe Is-Indonesia", -3.138, 125.491, -5, 1, 0));
    m_volcanoes.push_back(CVolcano("Awu", "", "0607-04=", "Sangihe Is-Indonesia", -3.67, 125.5, 1320, 1, 0));
    m_volcanoes.push_back(CVolcano("Tarakan", "", "0608-001", "Halmahera-Indonesia", -1.83, 127.83, 318, 1, 0));
    m_volcanoes.push_back(CVolcano("Dukono", "", "0608-01=", "Halmahera-Indonesia", -1.68, 127.88, 1335, 1, 0));
    m_volcanoes.push_back(CVolcano("Tobaru", "", "0608-02-", "Halmahera-Indonesia", -1.63, 127.67, 1035, 1, 0));
    m_volcanoes.push_back(CVolcano("Ibu", "", "0608-03=", "Halmahera-Indonesia", -1.488, 127.63, 1325, 1, 0));
    m_volcanoes.push_back(CVolcano("Gamkonora", "", "0608-04=", "Halmahera-Indonesia", -1.38, 127.53, 1635, 1, 0));
    m_volcanoes.push_back(CVolcano("Jailolo", "", "0608-051", "Halmahera-Indonesia", -1.08, 127.42, 1130, 1, 0));
    m_volcanoes.push_back(CVolcano("Hiri", "", "0608-052", "Halmahera-Indonesia", -0.9, 127.32, 630, 1, 0));
    m_volcanoes.push_back(CVolcano("Todoko-Ranu", "", "0608-05=", "Halmahera-Indonesia", -1.25, 127.47, 979, 1, 0));
    m_volcanoes.push_back(CVolcano("Tidore", "", "0608-061", "Halmahera-Indonesia", -0.658, 127.4, 1730, 1, 0));
    m_volcanoes.push_back(CVolcano("Mare", "", "0608-062", "Halmahera-Indonesia", -0.57, 127.4, 308, 1, 0));
    m_volcanoes.push_back(CVolcano("Moti", "", "0608-063", "Halmahera-Indonesia", -0.45, 127.4, 950, 1, 0));
    m_volcanoes.push_back(CVolcano("Gamalama", "", "0608-06=", "Halmahera-Indonesia", -0.8, 127.33, 1715, 1, 0));
    m_volcanoes.push_back(CVolcano("Tigalalu", "", "0608-071", "Halmahera-Indonesia", -0.07, 127.42, 422, 1, 0));
    m_volcanoes.push_back(CVolcano("Amasing", "", "0608-072", "Halmahera-Indonesia", 0.53, 127.48, 1030, 1, 0));
    m_volcanoes.push_back(CVolcano("Bibinoi", "", "0608-073", "Halmahera-Indonesia", 0.77, 127.72, 900, 1, 0));
    m_volcanoes.push_back(CVolcano("Makian", "", "0608-07=", "Halmahera-Indonesia", -0.32, 127.4, 1357, 1, 0));
    m_volcanoes.push_back(CVolcano("Bombalai", "", "0610-01-", "Borneo", -4.4, 117.88, 531, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_07()
{
    m_volcanoes.push_back(CVolcano("Jolo", "", "0700-01=", "Sulu Is-Philippines", -6.013, 121.057, 811, 1, 0));
    m_volcanoes.push_back(CVolcano("Parker", "", "0701-011", "Mindanao-Philippines", -6.113, 124.892, 1824, 1, 0));
    m_volcanoes.push_back(CVolcano("Balut", "", "0701-01=", "Mindanao-Philippines", -5.4, 125.375, 862, 1, 0));
    m_volcanoes.push_back(CVolcano("Matutum", "", "0701-02=", "Mindanao-Philippines", -6.37, 125.07, 2286, 1, 0));
    m_volcanoes.push_back(CVolcano("Leonard Range", "", "0701-031", "Mindanao-Philippines", -7.382, 126.047, 1080, 1, 0));
    m_volcanoes.push_back(CVolcano("Apo", "", "0701-03=", "Mindanao-Philippines", -6.989, 125.269, 2938, 1, 0));
    m_volcanoes.push_back(CVolcano("Makaturing", "", "0701-04=", "Mindanao-Philippines", -7.647, 124.32, 1940, 1, 0));
    m_volcanoes.push_back(CVolcano("Latukan", "", "0701-05=", "Mindanao-Philippines", -7.65, 124.45, 2338, 1, 0));
    m_volcanoes.push_back(CVolcano("Kalatungan", "", "0701-061", "Mindanao-Philippines", -7.95, 124.8, 2824, 1, 0));
    m_volcanoes.push_back(CVolcano("Ragang", "", "0701-06=", "Mindanao-Philippines", -7.7, 124.5, 2815, 1, 0));
    m_volcanoes.push_back(CVolcano("Malindang", "", "0701-071", "Mindanao-Philippines", -8.22, 123.63, 2404, 1, 0));
    m_volcanoes.push_back(CVolcano("Balatukan", "", "0701-072", "Mindanao-Philippines", -8.77, 124.98, 2450, 1, 0));
    m_volcanoes.push_back(CVolcano("Musuan", "", "0701-07=", "Mindanao-Philippines", -7.877, 125.068, 646, 1, 0));
    m_volcanoes.push_back(CVolcano("Camiguin", "", "0701-08=", "Mindanao-Philippines", -9.203, 124.673, 1552, 1, 0));
    m_volcanoes.push_back(CVolcano("Paco", "", "0701-09-", "Mindanao-Philippines", -9.593, 125.52, 524, 1, 0));
    m_volcanoes.push_back(CVolcano("Cuernos de Negros", "", "0702-01=", "Philippines-C", -9.25, 123.17, 1862, 1, 0));
    m_volcanoes.push_back(CVolcano("Kanlaon", "", "0702-02=", "Philippines-C", -10.412, 123.132, 2435, 1, 0));
    m_volcanoes.push_back(CVolcano("Mandalagan", "", "0702-03=", "Philippines-C", -10.65, 123.25, 1885, 1, 0));
    m_volcanoes.push_back(CVolcano("Silay", "", "0702-04=", "Philippines-C", -10.77, 123.23, 1510, 1, 0));
    m_volcanoes.push_back(CVolcano("Cabalían", "", "0702-05=", "Philippines-C", -10.287, 125.221, 945, 1, 0));
    m_volcanoes.push_back(CVolcano("Mahagnao", "", "0702-07=", "Philippines-C", -10.896, 125.87, 860, 1, 0));
    m_volcanoes.push_back(CVolcano("Biliran", "", "0702-08=", "Philippines-C", -11.523, 124.535, 1301, 1, 0));
    m_volcanoes.push_back(CVolcano("Bulusan", "", "0703-01=", "Luzon-Philippines", -12.77, 124.05, 1565, 1, 0));
    m_volcanoes.push_back(CVolcano("Pocdol Mountains", "", "0703-02=", "Luzon-Philippines", -13.05, 123.958, 1102, 1, 0));
    m_volcanoes.push_back(CVolcano("Masaraga", "", "0703-031", "Luzon-Philippines", -13.32, 123.6, 1328, 1, 0));
    m_volcanoes.push_back(CVolcano("Mayon", "", "0703-03=", "Luzon-Philippines", -13.257, 123.685, 2462, 1, 0));
    m_volcanoes.push_back(CVolcano("Iriga", "", "0703-041", "Luzon-Philippines", -13.457, 123.457, 1196, 1, 0));
    m_volcanoes.push_back(CVolcano("Isarog", "", "0703-042", "Luzon-Philippines", -13.658, 123.38, 1966, 1, 0));
    m_volcanoes.push_back(CVolcano("Malindig", "", "0703-044", "Luzon-Philippines", -13.24, 122.018, 1157, 1, 0));
    m_volcanoes.push_back(CVolcano("Banahaw", "", "0703-05=", "Luzon-Philippines", -14.07, 121.48, 2158, 1, 0));
    m_volcanoes.push_back(CVolcano("San Pablo Volc Field", "", "0703-06=", "Luzon-Philippines", -14.12, 121.3, 1090, 1, 0));
    m_volcanoes.push_back(CVolcano("Taal", "", "0703-07=", "Luzon-Philippines", -14.002, 120.993, 311, 1, 0));
    m_volcanoes.push_back(CVolcano("Mariveles", "", "0703-081", "Luzon-Philippines", -14.52, 120.47, 1388, 1, 0));
    m_volcanoes.push_back(CVolcano("Natib", "", "0703-082", "Luzon-Philippines", -14.72, 120.4, 1253, 1, 0));
    m_volcanoes.push_back(CVolcano("Pinatubo", "", "0703-083", "Luzon-Philippines", -15.13, 120.35, 1486, 1, 0));
    m_volcanoes.push_back(CVolcano("Arayat", "", "0703-084", "Luzon-Philippines", -15.2, 120.742, 1026, 1, 0));
    m_volcanoes.push_back(CVolcano("Amorong", "", "0703-085", "Luzon-Philippines", -15.828, 120.805, 376, 1, 0));
    m_volcanoes.push_back(CVolcano("Santo Tomas", "", "0703-086", "Luzon-Philippines", -16.33, 120.55, 2260, 1, 0));
    m_volcanoes.push_back(CVolcano("Patoc", "", "0703-087", "Luzon-Philippines", -17.147, 120.98, 1865, 1, 0));
    m_volcanoes.push_back(CVolcano("Ambalatungan Group", "", "0703-088", "Luzon-Philippines", -17.32, 121.1, 2329, 1, 0));
    m_volcanoes.push_back(CVolcano("Laguna Caldera", "", "0703-08=", "Luzon-Philippines", -14.42, 121.27, 743, 1, 0));
    m_volcanoes.push_back(CVolcano("Cagua", "", "0703-09=", "Luzon-Philippines", -18.222, 122.123, 1133, 1, 0));
    m_volcanoes.push_back(CVolcano("Camiguin de Babuyanes", "", "0704-01=", "Luzon-N of", -18.83, 121.86, 712, 1, 0));
    m_volcanoes.push_back(CVolcano("Didicas", "", "0704-02=", "Luzon-N of", -19.077, 122.202, 228, 1, 0));
    m_volcanoes.push_back(CVolcano("Babuyan Claro", "", "0704-03=", "Luzon-N of", -19.523, 121.94, 1080, 1, 0));
    m_volcanoes.push_back(CVolcano("Iraya", "", "0704-06-", "Luzon-N of", -20.469, 122.01, 1009, 1, 0));
    m_volcanoes.push_back(CVolcano("Hainan Dao", "", "0705-001", "SE Asia", -19.7, 110.1, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Leizhou Bandao", "", "0705-01-", "SE Asia", -20.78, 110.17, 259, 1, 0));
    m_volcanoes.push_back(CVolcano("Cù-Lao Ré Group", "", "0705-02-", "SE Asia", -15.38, 109.12, 181, 1, 0));
    m_volcanoes.push_back(CVolcano("Toroeng Prong", "", "0705-03-", "SE Asia", -14.93, 108, 800, 1, 0));
    m_volcanoes.push_back(CVolcano("Haut Dong Nai", "", "0705-04-", "SE Asia", -11.6, 108.2, 1000, 1, 0));
    m_volcanoes.push_back(CVolcano("Bas Dong Nai", "", "0705-05-", "SE Asia", -10.8, 107.2, 392, 1, 0));
    m_volcanoes.push_back(CVolcano("Ile des Cendres", "", "0705-06-", "SE Asia", -10.158, 109.014, -20, 1, 0));
    m_volcanoes.push_back(CVolcano("Veteran", "", "0705-07-", "SE Asia", -9.83, 109.05, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Popa", "", "0705-08-", "SE Asia", -20.92, 95.25, 1518, 1, 0));
    m_volcanoes.push_back(CVolcano("Lower Chindwin", "", "0705-09-", "SE Asia", -22.28, 95.1, 385, 1, 0));
    m_volcanoes.push_back(CVolcano("Singu Plateau", "", "0705-10-", "SE Asia", -22.7, 95.98, 507, 1, 0));
    m_volcanoes.push_back(CVolcano("Tengchong", "", "0705-11-", "SE Asia", -25.23, 98.5, 2865, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_08()
{
    m_volcanoes.push_back(CVolcano("Kueishantao", "", "0801-031", "Taiwan-E of", -24.85, 121.92, 401, 1, 0));
    m_volcanoes.push_back(CVolcano("Zengyu", "", "0801-05=", "Taiwan-N of", -26.18, 122.458, -418, 1, 0));
    m_volcanoes.push_back(CVolcano("Iriomote-jima", "", "0802-01=", "Ryukyu Is", -24.558, 124, -200, 1, 0));
    m_volcanoes.push_back(CVolcano("Yokoate-jima", "", "0802-021", "Ryukyu Is", -28.797, 128.997, 495, 1, 0));
    m_volcanoes.push_back(CVolcano("Akuseki-jima", "", "0802-022", "Ryukyu Is", -29.461, 129.597, 584, 1, 0));
    m_volcanoes.push_back(CVolcano("Iwo-Tori-shima", "", "0802-02=", "Ryukyu Is", -27.877, 128.224, 212, 1, 0));
    m_volcanoes.push_back(CVolcano("Suwanose-jima", "", "0802-03=", "Ryukyu Is", -29.635, 129.716, 799, 1, 0));
    m_volcanoes.push_back(CVolcano("Kogaja-jima", "", "0802-041", "Ryukyu Is", -29.879, 129.625, 301, 1, 0));
    m_volcanoes.push_back(CVolcano("Kuchino-shima", "", "0802-043", "Ryukyu Is", -29.964, 129.927, 628, 1, 0));
    m_volcanoes.push_back(CVolcano("Nakano-shima", "", "0802-04=", "Ryukyu Is", -29.856, 129.859, 979, 1, 0));
    m_volcanoes.push_back(CVolcano("Kuchinoerabu-jima", "", "0802-05=", "Ryukyu Is", -30.44, 130.219, 657, 1, 0));
    m_volcanoes.push_back(CVolcano("Kikai", "", "0802-06=", "Ryukyu Is", -30.789, 130.308, 704, 1, 0));
    m_volcanoes.push_back(CVolcano("Ibusuki Volc Field", "", "0802-07=", "Kyushu-Japan", -31.22, 130.57, 922, 1, 0));
    m_volcanoes.push_back(CVolcano("Sumiyoshi-ike", "", "0802-081", "Kyushu-Japan", -31.768, 130.594, 15, 1, 0));
    m_volcanoes.push_back(CVolcano("Sakura-jima", "", "0802-08=", "Kyushu-Japan", -31.585, 130.657, 1117, 1, 0));
    m_volcanoes.push_back(CVolcano("Fukue-jima", "", "0802-091", "Kyushu-Japan", -32.653, 128.851, 317, 1, 0));
    m_volcanoes.push_back(CVolcano("Kirishima", "", "0802-09=", "Kyushu-Japan", -31.931, 130.864, 1700, 1, 0));
    m_volcanoes.push_back(CVolcano("Unzen", "", "0802-10=", "Kyushu-Japan", -32.757, 130.294, 1500, 1, 0));
    m_volcanoes.push_back(CVolcano("Aso", "", "0802-11=", "Kyushu-Japan", -32.881, 131.106, 1592, 1, 0));
    m_volcanoes.push_back(CVolcano("Kuju", "", "0802-12=", "Kyushu-Japan", -33.083, 131.251, 1791, 1, 0));
    m_volcanoes.push_back(CVolcano("Tsurumi", "", "0802-13=", "Kyushu-Japan", -33.28, 131.432, 1584, 1, 0));
    m_volcanoes.push_back(CVolcano("Abu", "", "0803-001", "Honshu-Japan", -34.5, 131.6, 641, 1, 0));
    m_volcanoes.push_back(CVolcano("Sanbe", "", "0803-002", "Honshu-Japan", -35.13, 132.62, 1126, 1, 0));
    m_volcanoes.push_back(CVolcano("Oki-Dogo", "", "0803-003", "Honshu-Japan", -36.176, 133.334, 151, 1, 0));
    m_volcanoes.push_back(CVolcano("Izu-Tobu", "", "0803-01=", "Honshu-Japan", -34.9, 139.098, 1406, 1, 0));
    m_volcanoes.push_back(CVolcano("Hakone", "", "0803-02=", "Honshu-Japan", -35.23, 139.024, 1438, 1, 0));
    m_volcanoes.push_back(CVolcano("Kita Yatsuga-take", "", "0803-031", "Honshu-Japan", -36.1, 138.3, 2530, 1, 0));
    m_volcanoes.push_back(CVolcano("Fuji", "", "0803-03=", "Honshu-Japan", -35.358, 138.731, 3776, 1, 0));
    m_volcanoes.push_back(CVolcano("On-take", "", "0803-04=", "Honshu-Japan", -35.89, 137.48, 3063, 1, 0));
    m_volcanoes.push_back(CVolcano("Haku-san", "", "0803-05=", "Honshu-Japan", -36.152, 136.774, 2702, 1, 0));
    m_volcanoes.push_back(CVolcano("Norikura", "", "0803-06=", "Honshu-Japan", -36.103, 137.557, 3026, 1, 0));
    m_volcanoes.push_back(CVolcano("Washiba-Kumonotaira", "", "0803-071", "Honshu-Japan", -36.408, 137.594, 2924, 1, 0));
    m_volcanoes.push_back(CVolcano("Yake-dake", "", "0803-07=", "Honshu-Japan", -36.224, 137.59, 2455, 1, 0));
    m_volcanoes.push_back(CVolcano("Tate-yama", "", "0803-08=", "Honshu-Japan", -36.568, 137.593, 2621, 1, 0));
    m_volcanoes.push_back(CVolcano("Niigata-Yake-yama", "", "0803-09=", "Honshu-Japan", -36.918, 138.039, 2400, 1, 0));
    m_volcanoes.push_back(CVolcano("Myoko", "", "0803-10=", "Honshu-Japan", -36.888, 138.12, 2446, 1, 0));
    m_volcanoes.push_back(CVolcano("Asama", "", "0803-11=", "Honshu-Japan", -36.403, 138.526, 2568, 1, 0));
    m_volcanoes.push_back(CVolcano("Shiga", "", "0803-121", "Honshu-Japan", -36.688, 138.519, 2041, 1, 0));
    m_volcanoes.push_back(CVolcano("Haruna", "", "0803-122", "Honshu-Japan", -36.474, 138.881, 1449, 1, 0));
    m_volcanoes.push_back(CVolcano("Kusatsu-Shirane", "", "0803-12=", "Honshu-Japan", -36.62, 138.535, 2171, 1, 0));
    m_volcanoes.push_back(CVolcano("Hiuchi", "", "0803-131", "Honshu-Japan", -36.952, 139.289, 2356, 1, 0));
    m_volcanoes.push_back(CVolcano("Akagi", "", "0803-13=", "Honshu-Japan", -36.557, 139.196, 1828, 1, 0));
    m_volcanoes.push_back(CVolcano("Nantai", "", "0803-141", "Honshu-Japan", -36.762, 139.494, 2486, 1, 0));
    m_volcanoes.push_back(CVolcano("Omanago Group", "", "0803-142", "Honshu-Japan", -36.792, 139.51, 2367, 1, 0));
    m_volcanoes.push_back(CVolcano("Takahara", "", "0803-143", "Honshu-Japan", -36.897, 139.78, 1795, 1, 0));
    m_volcanoes.push_back(CVolcano("Nikko-Shirane", "", "0803-14=", "Honshu-Japan", -36.796, 139.379, 2578, 1, 0));
    m_volcanoes.push_back(CVolcano("Numazawa", "", "0803-151", "Honshu-Japan", -37.45, 139.579, 1100, 1, 0));
    m_volcanoes.push_back(CVolcano("Nasu", "", "0803-15=", "Honshu-Japan", -37.122, 139.966, 1915, 1, 0));
    m_volcanoes.push_back(CVolcano("Bandai", "", "0803-16=", "Honshu-Japan", -37.598, 140.076, 1819, 1, 0));
    m_volcanoes.push_back(CVolcano("Adatara", "", "0803-17=", "Honshu-Japan", -37.644, 140.286, 1718, 1, 0));
    m_volcanoes.push_back(CVolcano("Azuma", "", "0803-18=", "Honshu-Japan", -37.732, 140.248, 2035, 1, 0));
    m_volcanoes.push_back(CVolcano("Hijiori", "", "0803-191", "Honshu-Japan", -38.606, 140.178, 516, 1, 0));
    m_volcanoes.push_back(CVolcano("Zao", "", "0803-19=", "Honshu-Japan", -38.141, 140.443, 1841, 1, 0));
    m_volcanoes.push_back(CVolcano("Narugo", "", "0803-20=", "Honshu-Japan", -38.733, 140.732, 470, 1, 0));
    m_volcanoes.push_back(CVolcano("Kurikoma", "", "0803-21=", "Honshu-Japan", -38.958, 140.792, 1628, 1, 0));
    m_volcanoes.push_back(CVolcano("Chokai", "", "0803-22=", "Honshu-Japan", -39.096, 140.052, 2233, 1, 0));
    m_volcanoes.push_back(CVolcano("Akita-Komaga-take", "", "0803-23=", "Honshu-Japan", -39.758, 140.803, 1637, 1, 0));
    m_volcanoes.push_back(CVolcano("Iwate", "", "0803-24=", "Honshu-Japan", -39.85, 141.004, 2041, 1, 0));
    m_volcanoes.push_back(CVolcano("Hachimantai", "", "0803-25=", "Honshu-Japan", -39.955, 140.857, 1614, 1, 0));
    m_volcanoes.push_back(CVolcano("Megata", "", "0803-262", "Honshu-Japan", -39.95, 139.73, 291, 1, 0));
    m_volcanoes.push_back(CVolcano("Akita-Yake-yama", "", "0803-26=", "Honshu-Japan", -39.961, 140.761, 1366, 1, 0));
    m_volcanoes.push_back(CVolcano("Towada", "", "0803-271", "Honshu-Japan", -40.47, 140.92, 1159, 1, 0));
    m_volcanoes.push_back(CVolcano("Iwaki", "", "0803-27=", "Honshu-Japan", -40.653, 140.307, 1625, 1, 0));
    m_volcanoes.push_back(CVolcano("Hakkoda Group", "", "0803-28=", "Honshu-Japan", -40.656, 140.881, 1585, 1, 0));
    m_volcanoes.push_back(CVolcano("Osore-yama", "", "0803-29=", "Honshu-Japan", -41.276, 141.124, 879, 1, 0));
    m_volcanoes.push_back(CVolcano("To-shima", "", "0804-011", "Izu Is-Japan", -34.517, 139.283, 508, 1, 0));
    m_volcanoes.push_back(CVolcano("Oshima", "", "0804-01=", "Izu Is-Japan", -34.721, 139.398, 764, 1, 0));
    m_volcanoes.push_back(CVolcano("Nii-jima", "", "0804-02=", "Izu Is-Japan", -34.393, 139.273, 432, 1, 0));
    m_volcanoes.push_back(CVolcano("Kozu-shima", "", "0804-03=", "Izu Is-Japan", -34.216, 139.156, 572, 1, 0));
    m_volcanoes.push_back(CVolcano("Mikura-jima", "", "0804-041", "Izu Is-Japan", -33.871, 139.605, 851, 1, 0));
    m_volcanoes.push_back(CVolcano("Kurose Hole", "", "0804-042", "Izu Is-Japan", -33.4, 139.68, -107, 1, 0));
    m_volcanoes.push_back(CVolcano("Miyake-jima", "", "0804-04=", "Izu Is-Japan", -34.079, 139.529, 815, 1, 0));
    m_volcanoes.push_back(CVolcano("Hachijo-jima", "", "0804-05=", "Izu Is-Japan", -33.13, 139.769, 854, 1, 0));
    m_volcanoes.push_back(CVolcano("Myojin Knoll", "", "0804-061", "Izu Is-Japan", -32.1, 139.85, 360, 1, 0));
    m_volcanoes.push_back(CVolcano("Aoga-shima", "", "0804-06=", "Izu Is-Japan", -32.454, 139.762, 423, 1, 0));
    m_volcanoes.push_back(CVolcano("Bayonnaise Rocks", "", "0804-07=", "Izu Is-Japan", -31.88, 139.92, 11, 1, 0));
    m_volcanoes.push_back(CVolcano("Smith Rock", "", "0804-08=", "Izu Is-Japan", -31.436, 140.054, 136, 1, 0));
    m_volcanoes.push_back(CVolcano("Sofugan", "", "0804-091", "Izu Is-Japan", -29.789, 140.345, 99, 1, 0));
    m_volcanoes.push_back(CVolcano("Suiyo Seamount", "", "0804-093", "Izu Is-Japan", -28.6, 140.63, -1418, 1, 0));
    m_volcanoes.push_back(CVolcano("Mokuyo Seamount", "", "0804-094", "Izu Is-Japan", -28.32, 140.57, -920, 1, 0));
    m_volcanoes.push_back(CVolcano("Doyo Seamount", "", "0804-095", "Izu Is-Japan", -27.68, 140.8, -860, 1, 0));
    m_volcanoes.push_back(CVolcano("Nishino-shima", "", "0804-096", "Volcano Is-Japan", -27.274, 140.882, 38, 1, 0));
    m_volcanoes.push_back(CVolcano("Kaikata Seamount", "", "0804-097", "Volcano Is-Japan", -26.67, 141, -162, 1, 0));
    m_volcanoes.push_back(CVolcano("Tori-shima", "", "0804-09=", "Izu Is-Japan", -30.48, 140.306, 394, 1, 0));
    m_volcanoes.push_back(CVolcano("Kaitoku Seamount", "", "0804-10=", "Volcano Is-Japan", -26.122, 141.102, -103, 1, 0));
    m_volcanoes.push_back(CVolcano("Kita-Iwo-jima", "", "0804-11=", "Volcano Is-Japan", -25.424, 141.284, 792, 1, 0));
    m_volcanoes.push_back(CVolcano("Kita-Fukutokutai", "", "0804-121", "Volcano Is-Japan", -24.414, 141.419, -73, 1, 0));
    m_volcanoes.push_back(CVolcano("Ioto [Iwo-jima]", "", "0804-12=", "Volcano Is-Japan", -24.754, 141.29, 161, 1, 0));
    m_volcanoes.push_back(CVolcano("Minami-Hiyoshi", "", "0804-131", "Volcano Is-Japan", -23.497, 141.94, -30, 1, 0));
    m_volcanoes.push_back(CVolcano("Nikko", "", "0804-132", "Volcano Is-Japan", -23.075, 142.308, -391, 1, 0));
    m_volcanoes.push_back(CVolcano("Fukujin", "", "0804-133", "Volcano Is-Japan", -21.93, 143.47, -217, 1, 0));
    m_volcanoes.push_back(CVolcano("Kasuga", "", "0804-134", "Volcano Is-Japan", -21.765, 143.71, -598, 1, 0));
    m_volcanoes.push_back(CVolcano("Minami Kasuga", "", "0804-135", "Volcano Is-Japan", -21.6, 143.637, -274, 1, 0));
    m_volcanoes.push_back(CVolcano("NW Eifuku", "", "0804-136", "Volcano Is-Japan", -21.485, 144.043, -1535, 1, 0));
    m_volcanoes.push_back(CVolcano("Daikoku", "", "0804-137", "Volcano Is-Japan", -21.324, 144.194, -323, 1, 0));
    m_volcanoes.push_back(CVolcano("Fukutoku-Okanoba", "", "0804-13=", "Volcano Is-Japan", -24.28, 141.485, -14, 1, 0));
    m_volcanoes.push_back(CVolcano("Ahyi", "", "0804-141", "Mariana Is-C Pacific", -20.42, 145.03, -137, 1, 0));
    m_volcanoes.push_back(CVolcano("Supply Reef", "", "0804-142", "Mariana Is-C Pacific", -20.13, 145.1, -8, 1, 0));
    m_volcanoes.push_back(CVolcano("Maug Islands", "", "0804-143", "Mariana Is-C Pacific", -20.02, 145.22, 227, 1, 0));
    m_volcanoes.push_back(CVolcano("Farallon de Pajaros", "", "0804-14=", "Mariana Is-C Pacific", -20.538, 144.896, 360, 1, 0));
    m_volcanoes.push_back(CVolcano("Asuncion", "", "0804-15=", "Mariana Is-C Pacific", -19.671, 145.406, 857, 1, 0));
    m_volcanoes.push_back(CVolcano("Agrigan", "", "0804-16=", "Mariana Is-C Pacific", -18.77, 145.67, 965, 1, 0));
    m_volcanoes.push_back(CVolcano("Pagan", "", "0804-17=", "Mariana Is-C Pacific", -18.13, 145.8, 570, 1, 0));
    m_volcanoes.push_back(CVolcano("Alamagan", "", "0804-18=", "Mariana Is-C Pacific", -17.6, 145.83, 744, 1, 0));
    m_volcanoes.push_back(CVolcano("Zealandia Bank", "", "0804-191", "Mariana Is-C Pacific", -16.88, 145.85, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Sarigan", "", "0804-192", "Mariana Is-C Pacific", -16.708, 145.78, 538, 1, 0));
    m_volcanoes.push_back(CVolcano("Guguan", "", "0804-19=", "Mariana Is-C Pacific", -17.307, 145.845, 287, 1, 0));
    m_volcanoes.push_back(CVolcano("East Diamante", "", "0804-201", "Mariana Is-C Pacific", -15.93, 145.67, -127, 1, 0));
    m_volcanoes.push_back(CVolcano("Ruby", "", "0804-202", "Mariana Is-C Pacific", -15.62, 145.57, -230, 1, 0));
    m_volcanoes.push_back(CVolcano("Anatahan", "", "0804-20=", "Mariana Is-C Pacific", -16.35, 145.67, 790, 1, 0));
    m_volcanoes.push_back(CVolcano("NW Rota-1", "", "0804-211", "Mariana Is-C Pacific", -14.601, 144.775, -517, 1, 0));
    m_volcanoes.push_back(CVolcano("Esmeralda Bank", "", "0804-21=", "Mariana Is-C Pacific", -15, 145.25, -43, 1, 0));
    m_volcanoes.push_back(CVolcano("Forecast Seamount", "", "0804-22-", "Mariana Is-C Pacific", -13.4, 143.92, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Seamount X", "", "0804-23-", "Mariana Is-C Pacific", -13.25, 144.02, -1230, 1, 0));
    m_volcanoes.push_back(CVolcano("E-san", "", "0805-011", "Hokkaido-Japan", -41.802, 141.17, 618, 1, 0));
    m_volcanoes.push_back(CVolcano("Oshima-Oshima", "", "0805-01=", "Hokkaido-Japan", -41.507, 139.371, 737, 1, 0));
    m_volcanoes.push_back(CVolcano("Komaga-take", "", "0805-02=", "Hokkaido-Japan", -42.061, 140.681, 1131, 1, 0));
    m_volcanoes.push_back(CVolcano("Niseko", "", "0805-031", "Hokkaido-Japan", -42.88, 140.63, 1154, 1, 0));
    m_volcanoes.push_back(CVolcano("Yotei", "", "0805-032", "Hokkaido-Japan", -42.83, 140.815, 1898, 1, 0));
    m_volcanoes.push_back(CVolcano("Kuttara", "", "0805-034", "Hokkaido-Japan", -42.489, 141.163, 581, 1, 0));
    m_volcanoes.push_back(CVolcano("Usu", "", "0805-03=", "Hokkaido-Japan", -42.541, 140.843, 737, 1, 0));
    m_volcanoes.push_back(CVolcano("Rishiri", "", "0805-041", "Hokkaido-Japan", -45.18, 141.25, 1721, 1, 0));
    m_volcanoes.push_back(CVolcano("Shikotsu", "", "0805-04=", "Hokkaido-Japan", -42.688, 141.38, 1320, 1, 0));
    m_volcanoes.push_back(CVolcano("Tokachi", "", "0805-05=", "Hokkaido-Japan", -43.416, 142.69, 2077, 1, 0));
    m_volcanoes.push_back(CVolcano("Nipesotsu-Maruyama", "", "0805-061", "Hokkaido-Japan", -43.453, 143.036, 2013, 1, 0));
    m_volcanoes.push_back(CVolcano("Shikaribetsu Group", "", "0805-062", "Hokkaido-Japan", -43.312, 143.096, 1401, 1, 0));
    m_volcanoes.push_back(CVolcano("Daisetsu", "", "0805-06=", "Hokkaido-Japan", -43.661, 142.858, 2290, 1, 0));
    m_volcanoes.push_back(CVolcano("Akan", "", "0805-07=", "Hokkaido-Japan", -43.384, 144.013, 1499, 1, 0));
    m_volcanoes.push_back(CVolcano("Mashu", "", "0805-081", "Hokkaido-Japan", -43.57, 144.565, 855, 1, 0));
    m_volcanoes.push_back(CVolcano("Rausu", "", "0805-082", "Hokkaido-Japan", -44.073, 145.126, 1660, 1, 0));
    m_volcanoes.push_back(CVolcano("Kutcharo", "", "0805-08=", "Hokkaido-Japan", -43.608, 144.443, 999, 1, 0));
    m_volcanoes.push_back(CVolcano("Shiretoko-Iwo-zan", "", "0805-09=", "Hokkaido-Japan", -44.131, 145.165, 1563, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_09()
{
    m_volcanoes.push_back(CVolcano("Golovnin", "", "0900-01=", "Kuril Is", -43.841, 145.509, 543, 1, 0));
    m_volcanoes.push_back(CVolcano("Smirnov", "", "0900-021", "Kuril Is", -44.42, 146.135, 1189, 1, 0));
    m_volcanoes.push_back(CVolcano("Mendeleev", "", "0900-02=", "Kuril Is", -43.976, 145.736, 888, 1, 0));
    m_volcanoes.push_back(CVolcano("Tiatia", "", "0900-03=", "Kuril Is", -44.351, 146.256, 1819, 1, 0));
    m_volcanoes.push_back(CVolcano("Lvinaya Past", "", "0900-041", "Kuril Is", -44.608, 146.994, 528, 1, 0));
    m_volcanoes.push_back(CVolcano("Berutarube", "", "0900-04=", "Kuril Is", -44.459, 146.936, 1221, 1, 0));
    m_volcanoes.push_back(CVolcano("Atsonupuri", "", "0900-05=", "Kuril Is", -44.805, 147.135, 1206, 1, 0));
    m_volcanoes.push_back(CVolcano("Bogatyr Ridge", "", "0900-06-", "Kuril Is", -44.833, 147.342, 1634, 1, 0));
    m_volcanoes.push_back(CVolcano("Grozny Group", "", "0900-07=", "Kuril Is", -45.026, 147.922, 1211, 1, 0));
    m_volcanoes.push_back(CVolcano("Baransky", "", "0900-08=", "Kuril Is", -45.097, 148.024, 1132, 1, 0));
    m_volcanoes.push_back(CVolcano("Golets-Tornyi Group", "", "0900-091", "Kuril Is", -45.25, 148.35, 442, 1, 0));
    m_volcanoes.push_back(CVolcano("Chirip", "", "0900-09=", "Kuril Is", -45.338, 147.925, 1587, 1, 0));
    m_volcanoes.push_back(CVolcano("Medvezhia", "", "0900-10=", "Kuril Is", -45.387, 148.843, 1125, 1, 0));
    m_volcanoes.push_back(CVolcano("Demon", "", "0900-11-", "Kuril Is", -45.5, 148.85, 1205, 1, 0));
    m_volcanoes.push_back(CVolcano("Ivao Group", "", "0900-111", "Kuril Is", -45.77, 149.68, 1426, 1, 0));
    m_volcanoes.push_back(CVolcano("Rudakov", "", "0900-112", "Kuril Is", -45.88, 149.83, 542, 1, 0));
    m_volcanoes.push_back(CVolcano("Tri Sestry", "", "0900-113", "Kuril Is", -45.93, 149.92, 998, 1, 0));
    m_volcanoes.push_back(CVolcano("Kolokol Group", "", "0900-12=", "Kuril Is", -46.042, 150.05, 1328, 1, 0));
    m_volcanoes.push_back(CVolcano("Chirpoi", "", "0900-15=", "Kuril Is", -46.525, 150.875, 742, 1, 0));
    m_volcanoes.push_back(CVolcano("Milne", "", "0900-161", "Kuril Is", -46.82, 151.78, 1540, 1, 0));
    m_volcanoes.push_back(CVolcano("Goriaschaia Sopka", "", "0900-17=", "Kuril Is", -46.83, 151.75, 891, 1, 0));
    m_volcanoes.push_back(CVolcano("Zavaritzki Caldera", "", "0900-18=", "Kuril Is", -46.925, 151.95, 624, 1, 0));
    m_volcanoes.push_back(CVolcano("Urataman", "", "0900-191", "Kuril Is", -47.12, 152.25, 678, 1, 0));
    m_volcanoes.push_back(CVolcano("Prevo Peak", "", "0900-19=", "Kuril Is", -47.02, 152.12, 1360, 1, 0));
    m_volcanoes.push_back(CVolcano("Ketoi", "", "0900-20=", "Kuril Is", -47.35, 152.475, 1172, 1, 0));
    m_volcanoes.push_back(CVolcano("Srednii", "", "0900-211", "Kuril Is", -47.6, 152.92, 36, 1, 0));
    m_volcanoes.push_back(CVolcano("Ushishur", "", "0900-21=", "Kuril Is", -47.52, 152.8, 401, 1, 0));
    m_volcanoes.push_back(CVolcano("Rasshua", "", "0900-22=", "Kuril Is", -47.77, 153.02, 956, 1, 0));
    m_volcanoes.push_back(CVolcano("Sarychev Peak", "", "0900-24=", "Kuril Is", -48.092, 153.2, 1496, 1, 0));
    m_volcanoes.push_back(CVolcano("Raikoke", "", "0900-25=", "Kuril Is", -48.292, 153.25, 551, 1, 0));
    m_volcanoes.push_back(CVolcano("Chirinkotan", "", "0900-26=", "Kuril Is", -48.98, 153.48, 724, 1, 0));
    m_volcanoes.push_back(CVolcano("Ekarma", "", "0900-27=", "Kuril Is", -48.958, 153.93, 1170, 1, 0));
    m_volcanoes.push_back(CVolcano("Sinarka", "", "0900-29=", "Kuril Is", -48.875, 154.175, 934, 1, 0));
    m_volcanoes.push_back(CVolcano("Kharimkotan", "", "0900-30=", "Kuril Is", -49.12, 154.508, 1145, 1, 0));
    m_volcanoes.push_back(CVolcano("Tao-Rusyr Caldera", "", "0900-31=", "Kuril Is", -49.35, 154.7, 1325, 1, 0));
    m_volcanoes.push_back(CVolcano("Nemo Peak", "", "0900-32=", "Kuril Is", -49.57, 154.808, 1018, 1, 0));
    m_volcanoes.push_back(CVolcano("Shirinki", "", "0900-331", "Kuril Is", -50.2, 154.98, 761, 1, 0));
    m_volcanoes.push_back(CVolcano("Fuss Peak", "", "0900-34=", "Kuril Is", -50.27, 155.25, 1772, 1, 0));
    m_volcanoes.push_back(CVolcano("Lomonosov Group", "", "0900-351", "Kuril Is", -50.25, 155.43, 1681, 1, 0));
    m_volcanoes.push_back(CVolcano("Karpinsky Group", "", "0900-35=", "Kuril Is", -50.13, 155.37, 1345, 1, 0));
    m_volcanoes.push_back(CVolcano("Chikurachki", "", "0900-36=", "Kuril Is", -50.325, 155.458, 1816, 1, 0));
    m_volcanoes.push_back(CVolcano("Vernadskii Ridge", "", "0900-37-", "Kuril Is", -50.55, 155.97, 1183, 1, 0));
    m_volcanoes.push_back(CVolcano("Ebeko", "", "0900-38=", "Kuril Is", -50.68, 156.02, 1156, 1, 0));
    m_volcanoes.push_back(CVolcano("Alaid", "", "0900-39=", "Kuril Is", -50.858, 155.55, 2339, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_10()
{
    m_volcanoes.push_back(CVolcano("Mashkovtsev", "", "1000-001", "Kamchatka", -51.1, 156.72, 503, 1, 0));
    m_volcanoes.push_back(CVolcano("Kambalny", "", "1000-01=", "Kamchatka", -51.3, 156.87, 2156, 1, 0));
    m_volcanoes.push_back(CVolcano("Yavinsky", "", "1000-021", "Kamchatka", -51.57, 156.6, 705, 1, 0));
    m_volcanoes.push_back(CVolcano("Diky Greben", "", "1000-022", "Kamchatka", -51.45, 156.97, 1070, 1, 0));
    m_volcanoes.push_back(CVolcano("Kurile Lake", "", "1000-023", "Kamchatka", -51.45, 157.12, 81, 1, 0));
    m_volcanoes.push_back(CVolcano("Koshelev", "", "1000-02=", "Kamchatka", -51.357, 156.75, 1812, 1, 0));
    m_volcanoes.push_back(CVolcano("Ilyinsky", "", "1000-03=", "Kamchatka", -51.49, 157.2, 1578, 1, 0));
    m_volcanoes.push_back(CVolcano("Kell", "", "1000-041", "Kamchatka", -51.65, 157.35, 900, 1, 0));
    m_volcanoes.push_back(CVolcano("Belenkaya", "", "1000-042", "Kamchatka", -51.75, 157.27, 892, 1, 0));
    m_volcanoes.push_back(CVolcano("Zheltovsky", "", "1000-04=", "Kamchatka", -51.57, 157.323, 1953, 1, 0));
    m_volcanoes.push_back(CVolcano("Ozernoy", "", "1000-051", "Kamchatka", -51.88, 157.38, 562, 1, 0));
    m_volcanoes.push_back(CVolcano("Olkoviy Volc Group", "", "1000-052", "Kamchatka", -52.02, 157.53, 681, 1, 0));
    m_volcanoes.push_back(CVolcano("Khodutka", "", "1000-053", "Kamchatka", -52.063, 157.703, 2090, 1, 0));
    m_volcanoes.push_back(CVolcano("Piratkovsky", "", "1000-054", "Kamchatka", -52.113, 157.849, 1322, 1, 0));
    m_volcanoes.push_back(CVolcano("Ostanets", "", "1000-055", "Kamchatka", -52.146, 157.322, 719, 1, 0));
    m_volcanoes.push_back(CVolcano("Otdelniy", "", "1000-056", "Kamchatka", -52.22, 157.428, 791, 1, 0));
    m_volcanoes.push_back(CVolcano("Golaya", "", "1000-057", "Kamchatka", -52.263, 157.787, 858, 1, 0));
    m_volcanoes.push_back(CVolcano("Asacha", "", "1000-058", "Kamchatka", -52.355, 157.827, 1910, 1, 0));
    m_volcanoes.push_back(CVolcano("Visokiy", "", "1000-059", "Kamchatka", -52.43, 157.93, 1234, 1, 0));
    m_volcanoes.push_back(CVolcano("Ksudach", "", "1000-05=", "Kamchatka", -51.8, 157.53, 1079, 1, 0));
    m_volcanoes.push_back(CVolcano("Mutnovsky", "", "1000-06=", "Kamchatka", -52.453, 158.195, 2322, 1, 0));
    m_volcanoes.push_back(CVolcano("Gorely", "", "1000-07=", "Kamchatka", -52.558, 158.03, 1829, 1, 0));
    m_volcanoes.push_back(CVolcano("Tolmachev Dol", "", "1000-082", "Kamchatka", -52.63, 157.58, 1021, 1, 0));
    m_volcanoes.push_back(CVolcano("Vilyuchik", "", "1000-083", "Kamchatka", -52.7, 158.28, 2173, 1, 0));
    m_volcanoes.push_back(CVolcano("Barkhatnaya Sopka", "", "1000-084", "Kamchatka", -52.823, 158.27, 870, 1, 0));
    m_volcanoes.push_back(CVolcano("Bolshe-Bannaya", "", "1000-087", "Kamchatka", -52.9, 157.78, 1200, 1, 0));
    m_volcanoes.push_back(CVolcano("Opala", "", "1000-08=", "Kamchatka", -52.543, 157.335, 2475, 1, 0));
    m_volcanoes.push_back(CVolcano("Koryaksky", "", "1000-09=", "Kamchatka", -53.32, 158.688, 3456, 1, 0));
    m_volcanoes.push_back(CVolcano("Avachinsky", "", "1000-10=", "Kamchatka", -53.255, 158.83, 2741, 1, 0));
    m_volcanoes.push_back(CVolcano("Dzenzursky", "", "1000-11=", "Kamchatka", -53.637, 158.922, 2285, 1, 0));
    m_volcanoes.push_back(CVolcano("Veer", "", "1000-121", "Kamchatka", -53.75, 158.45, 520, 1, 0));
    m_volcanoes.push_back(CVolcano("Kostakan", "", "1000-122", "Kamchatka", -53.83, 158.05, 1150, 1, 0));
    m_volcanoes.push_back(CVolcano("Bakening", "", "1000-123", "Kamchatka", -53.905, 158.07, 2278, 1, 0));
    m_volcanoes.push_back(CVolcano("Zavaritsky", "", "1000-124", "Kamchatka", -53.905, 158.385, 1567, 1, 0));
    m_volcanoes.push_back(CVolcano("Akademia Nauk", "", "1000-125", "Kamchatka", -53.98, 159.45, 1180, 1, 0));
    m_volcanoes.push_back(CVolcano("Zhupanovsky", "", "1000-12=", "Kamchatka", -53.59, 159.147, 2958, 1, 0));
    m_volcanoes.push_back(CVolcano("Karymsky", "", "1000-13=", "Kamchatka", -54.05, 159.45, 1536, 1, 0));
    m_volcanoes.push_back(CVolcano("Maly Semiachik", "", "1000-14=", "Kamchatka", -54.13, 159.67, 1560, 1, 0));
    m_volcanoes.push_back(CVolcano("Bolshoi Semiachik", "", "1000-15=", "Kamchatka", -54.32, 160.02, 1720, 1, 0));
    m_volcanoes.push_back(CVolcano("Taunshits", "", "1000-16-", "Kamchatka", -54.53, 159.8, 2353, 1, 0));
    m_volcanoes.push_back(CVolcano("Uzon", "", "1000-17=", "Kamchatka", -54.5, 159.97, 1617, 1, 0));
    m_volcanoes.push_back(CVolcano("Kikhpinych", "", "1000-18=", "Kamchatka", -54.487, 160.253, 1552, 1, 0));
    m_volcanoes.push_back(CVolcano("Krasheninnikov", "", "1000-19=", "Kamchatka", -54.593, 160.273, 1856, 1, 0));
    m_volcanoes.push_back(CVolcano("Schmidt", "", "1000-201", "Kamchatka", -54.92, 160.63, 2020, 1, 0));
    m_volcanoes.push_back(CVolcano("Kronotsky", "", "1000-20=", "Kamchatka", -54.753, 160.527, 3528, 1, 0));
    m_volcanoes.push_back(CVolcano("Gamchen", "", "1000-21=", "Kamchatka", -54.973, 160.702, 2576, 1, 0));
    m_volcanoes.push_back(CVolcano("Vysoky", "", "1000-221", "Kamchatka", -55.07, 160.77, 2161, 1, 0));
    m_volcanoes.push_back(CVolcano("Komarov", "", "1000-22=", "Kamchatka", -55.032, 160.72, 2070, 1, 0));
    m_volcanoes.push_back(CVolcano("Kizimen", "", "1000-23=", "Kamchatka", -55.13, 160.32, 2376, 1, 0));
    m_volcanoes.push_back(CVolcano("Udina", "", "1000-241", "Kamchatka", -55.755, 160.527, 2923, 1, 0));
    m_volcanoes.push_back(CVolcano("Zimina", "", "1000-242", "Kamchatka", -55.862, 160.603, 3081, 1, 0));
    m_volcanoes.push_back(CVolcano("Tolbachik", "", "1000-24=", "Kamchatka", -55.83, 160.33, 3682, 1, 0));
    m_volcanoes.push_back(CVolcano("Kamen", "", "1000-251", "Kamchatka", -56.02, 160.593, 4585, 1, 0));
    m_volcanoes.push_back(CVolcano("Bezymianny", "", "1000-25=", "Kamchatka", -55.978, 160.587, 2882, 1, 0));
    m_volcanoes.push_back(CVolcano("Ushkovsky", "", "1000-261", "Kamchatka", -56.07, 160.47, 3943, 1, 0));
    m_volcanoes.push_back(CVolcano("Kliuchevskoi", "", "1000-26=", "Kamchatka", -56.057, 160.638, 4835, 1, 0));
    m_volcanoes.push_back(CVolcano("Piip", "", "1000-271", "Kamchatka-E of", -55.42, 167.33, -300, 1, 0));
    m_volcanoes.push_back(CVolcano("Khangar", "", "1000-272", "Kamchatka", -54.75, 157.38, 2000, 1, 0));
    m_volcanoes.push_back(CVolcano("Cherpuk Group", "", "1000-273", "Kamchatka", -55.55, 157.47, 1868, 1, 0));
    m_volcanoes.push_back(CVolcano("Shiveluch", "", "1000-27=", "Kamchatka", -56.653, 161.36, 3283, 1, 0));
    m_volcanoes.push_back(CVolcano("Ichinsky", "", "1000-28=", "Kamchatka", -55.68, 157.73, 3621, 1, 0));
    m_volcanoes.push_back(CVolcano("Maly Payalpan", "", "1000-29-", "Kamchatka", -55.82, 157.98, 1802, 1, 0));
    m_volcanoes.push_back(CVolcano("Bolshoi Payalpan", "", "1000-30-", "Kamchatka", -55.88, 157.78, 1906, 1, 0));
    m_volcanoes.push_back(CVolcano("Plosky", "", "1000-31-", "Kamchatka", -55.2, 158.47, 1236, 1, 0));
    m_volcanoes.push_back(CVolcano("Akhtang", "", "1000-32-", "Kamchatka", -55.43, 158.65, 1956, 1, 0));
    m_volcanoes.push_back(CVolcano("Kozyrevsky", "", "1000-33-", "Kamchatka", -55.58, 158.38, 2016, 1, 0));
    m_volcanoes.push_back(CVolcano("Romanovka", "", "1000-34-", "Kamchatka", -55.65, 158.8, 1442, 1, 0));
    m_volcanoes.push_back(CVolcano("Uksichan", "", "1000-35-", "Kamchatka", -56.08, 158.38, 1692, 1, 0));
    m_volcanoes.push_back(CVolcano("Bolshoi-Kekuknaysky", "", "1000-36-", "Kamchatka", -56.47, 157.8, 1401, 1, 0));
    m_volcanoes.push_back(CVolcano("Kulkev", "", "1000-37-", "Kamchatka", -56.37, 158.37, 915, 1, 0));
    m_volcanoes.push_back(CVolcano("Geodesistoy", "", "1000-38-", "Kamchatka", -56.33, 158.67, 1170, 1, 0));
    m_volcanoes.push_back(CVolcano("Anaun", "", "1000-39-", "Kamchatka", -56.32, 158.83, 1828, 1, 0));
    m_volcanoes.push_back(CVolcano("Krainy", "", "1000-40-", "Kamchatka", -56.37, 159.03, 1554, 1, 0));
    m_volcanoes.push_back(CVolcano("Kekurny", "", "1000-41-", "Kamchatka", -56.4, 158.85, 1377, 1, 0));
    m_volcanoes.push_back(CVolcano("Eggella", "", "1000-42-", "Kamchatka", -56.57, 158.52, 1046, 1, 0));
    m_volcanoes.push_back(CVolcano("Verkhovoy", "", "1000-44-", "Kamchatka", -56.52, 159.53, 1400, 1, 0));
    m_volcanoes.push_back(CVolcano("Alney-Chashakondzha", "", "1000-45-", "Kamchatka", -56.7, 159.65, 2598, 1, 0));
    m_volcanoes.push_back(CVolcano("Cherny", "", "1000-46-", "Kamchatka", -56.82, 159.67, 1778, 1, 0));
    m_volcanoes.push_back(CVolcano("Pogranychny", "", "1000-47-", "Kamchatka", -56.85, 159.8, 1427, 1, 0));
    m_volcanoes.push_back(CVolcano("Zaozerny", "", "1000-48-", "Kamchatka", -56.88, 159.95, 1349, 1, 0));
    m_volcanoes.push_back(CVolcano("Bliznets", "", "1000-49-", "Kamchatka", -56.97, 159.78, 1244, 1, 0));
    m_volcanoes.push_back(CVolcano("Kebeney", "", "1000-50-", "Kamchatka", -57.1, 159.93, 1527, 1, 0));
    m_volcanoes.push_back(CVolcano("Fedotych", "", "1000-51-", "Kamchatka", -57.13, 160.4, 965, 1, 0));
    m_volcanoes.push_back(CVolcano("Shisheika", "", "1000-511", "Kamchatka", -57.15, 161.08, 379, 1, 0));
    m_volcanoes.push_back(CVolcano("Terpuk", "", "1000-512", "Kamchatka", -57.2, 159.83, 765, 1, 0));
    m_volcanoes.push_back(CVolcano("Sedankinsky", "", "1000-52-", "Kamchatka", -57.27, 160.08, 1241, 1, 0));
    m_volcanoes.push_back(CVolcano("Leutongey", "", "1000-53-", "Kamchatka", -57.3, 159.83, 1333, 1, 0));
    m_volcanoes.push_back(CVolcano("Tuzovsky", "", "1000-54-", "Kamchatka", -57.32, 159.97, 1533, 1, 0));
    m_volcanoes.push_back(CVolcano("Gorny Institute", "", "1000-55-", "Kamchatka", -57.33, 160.2, 2125, 1, 0));
    m_volcanoes.push_back(CVolcano("Kinenin", "", "1000-551", "Kamchatka", -57.35, 160.97, 583, 1, 0));
    m_volcanoes.push_back(CVolcano("Bliznetsy", "", "1000-552", "Kamchatka", -57.35, 161.37, 265, 1, 0));
    m_volcanoes.push_back(CVolcano("Titila", "", "1000-56-", "Kamchatka", -57.4, 160.1, 1559, 1, 0));
    m_volcanoes.push_back(CVolcano("Mezhdusopochny", "", "1000-57-", "Kamchatka", -57.47, 160.25, 1641, 1, 0));
    m_volcanoes.push_back(CVolcano("Shishel", "", "1000-58-", "Kamchatka", -57.45, 160.37, 2525, 1, 0));
    m_volcanoes.push_back(CVolcano("Elovsky", "", "1000-59-", "Kamchatka", -57.55, 160.53, 1381, 1, 0));
    m_volcanoes.push_back(CVolcano("Alngey", "", "1000-60-", "Kamchatka", -57.7, 160.4, 1853, 1, 0));
    m_volcanoes.push_back(CVolcano("Uka", "", "1000-61-", "Kamchatka", -57.7, 160.58, 1643, 1, 0));
    m_volcanoes.push_back(CVolcano("Kaileney", "", "1000-62-", "Kamchatka", -57.8, 160.67, 1582, 1, 0));
    m_volcanoes.push_back(CVolcano("Plosky", "", "1000-63-", "Kamchatka", -57.83, 160.25, 1255, 1, 0));
    m_volcanoes.push_back(CVolcano("Bely", "", "1000-64-", "Kamchatka", -57.88, 160.53, 2080, 1, 0));
    m_volcanoes.push_back(CVolcano("Nylgimelkin", "", "1000-65-", "Kamchatka", -57.97, 160.65, 1764, 1, 0));
    m_volcanoes.push_back(CVolcano("Snezhniy", "", "1000-66-", "Kamchatka", -58.02, 160.8, 2169, 1, 0));
    m_volcanoes.push_back(CVolcano("Iktunup", "", "1000-67-", "Kamchatka", -58.08, 160.77, 2300, 1, 0));
    m_volcanoes.push_back(CVolcano("Spokoiny", "", "1000-671", "Kamchatka", -58.13, 160.82, 2171, 1, 0));
    m_volcanoes.push_back(CVolcano("Ostry", "", "1000-68-", "Kamchatka", -58.18, 160.82, 2552, 1, 0));
    m_volcanoes.push_back(CVolcano("Snegovoy", "", "1000-69-", "Kamchatka", -58.2, 160.97, 2169, 1, 0));
    m_volcanoes.push_back(CVolcano("Severny", "", "1000-70-", "Kamchatka", -58.28, 160.87, 1936, 1, 0));
    m_volcanoes.push_back(CVolcano("Iettunup", "", "1000-71-", "Kamchatka", -58.4, 161.08, 1340, 1, 0));
    m_volcanoes.push_back(CVolcano("Voyampolsky", "", "1000-72-", "Kamchatka", -58.37, 160.62, 1225, 1, 0));
    m_volcanoes.push_back(CVolcano("Sikhote-Alin", "", "1002-01-", "Russia-SE", -47, 137.5, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Udokan Plateau", "", "1002-03-", "Russia-SE", -56.28, 117.77, 2180, 1, 0));
    m_volcanoes.push_back(CVolcano("Vitim Plateau", "", "1002-04-", "Russia-SE", -53.7, 113.3, 1250, 1, 0));
    m_volcanoes.push_back(CVolcano("Tunkin Depression", "", "1002-05-", "Russia-SE", -51.5, 102.5, 1200, 1, 0));
    m_volcanoes.push_back(CVolcano("Oka Plateau", "", "1002-06-", "Russia-SE", -52.7, 98.98, 2077, 1, 0));
    m_volcanoes.push_back(CVolcano("Azas Plateau", "", "1002-07-", "Russia-SE", -52.52, 98.6, 2765, 1, 0));
    m_volcanoes.push_back(CVolcano("Taryatu-Chulutu", "", "1003-01-", "Mongolia", -48.17, 99.7, 2400, 1, 0));
    m_volcanoes.push_back(CVolcano("Khanuy Gol", "", "1003-02-", "Mongolia", -48.67, 102.75, 1886, 1, 0));
    m_volcanoes.push_back(CVolcano("Bus-Obo", "", "1003-03-", "Mongolia", -47.12, 109.08, 1162, 1, 0));
    m_volcanoes.push_back(CVolcano("Dariganga Volc Field", "", "1003-04-", "Mongolia", -45.33, 114, 1778, 1, 0));
    m_volcanoes.push_back(CVolcano("Middle Gobi", "", "1003-05-", "Mongolia", -45.28, 106.7, 1120, 1, 0));
    m_volcanoes.push_back(CVolcano("Turfan", "", "1004-01-", "China-W", -42.9, 89.25, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Tianshan Volc Group", "", "1004-02-", "China-W", -42.5, 82.5, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Kunlun Volc Group", "", "1004-03-", "China-W", -35.52, 80.2, 5808, 1, 0));
    m_volcanoes.push_back(CVolcano("Honggeertu", "", "1005-01-", "China-E", -41.47, 113, 1700, 1, 0));
    m_volcanoes.push_back(CVolcano("Arshan", "", "1005-011", "China-E", -47.5, 120.7, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Keluo Group", "", "1005-02-", "China-E", -49.37, 125.92, 670, 1, 0));
    m_volcanoes.push_back(CVolcano("Wudalianchi", "", "1005-03-", "China-E", -48.72, 126.12, 597, 1, 0));
    m_volcanoes.push_back(CVolcano("Jingbo", "", "1005-04-", "China-E", -44.08, 128.83, 1000, 1, 0));
    m_volcanoes.push_back(CVolcano("Longgang Group", "", "1005-05-", "China-E", -42.33, 126.5, 1000, 1, 0));
    m_volcanoes.push_back(CVolcano("Changbaishan", "", "1005-06-", "China-E", -41.98, 128.08, 2744, 1, 0));
    m_volcanoes.push_back(CVolcano("Xianjindao", "", "1006-01-", "Korea", -41.33, 128, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Ch'uga-ryong", "", "1006-02-", "Korea", -38.33, 127.33, 452, 1, 0));
    m_volcanoes.push_back(CVolcano("Ulreung", "", "1006-03-", "Korea", -37.5, 130.87, 984, 1, 0));
    m_volcanoes.push_back(CVolcano("Halla", "", "1006-04-", "Korea", -33.37, 126.53, 1950, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_11()
{
    m_volcanoes.push_back(CVolcano("Buldir", "", "1101-01-", "Aleutian Is", -52.35, 175.911, 656, 1, 0));
    m_volcanoes.push_back(CVolcano("Kiska", "", "1101-02-", "Aleutian Is", -52.103, 177.602, 1220, 1, 0));
    m_volcanoes.push_back(CVolcano("Segula", "", "1101-03-", "Aleutian Is", -52.015, 178.136, 1160, 1, 0));
    m_volcanoes.push_back(CVolcano("Davidof", "", "1101-04-", "Aleutian Is", -51.97, 178.33, 328, 1, 0));
    m_volcanoes.push_back(CVolcano("Little Sitkin", "", "1101-05-", "Aleutian Is", -51.95, 178.543, 1174, 1, 0));
    m_volcanoes.push_back(CVolcano("Semisopochnoi", "", "1101-06-", "Aleutian Is", -51.93, 179.58, 1221, 1, 0));
    m_volcanoes.push_back(CVolcano("Gareloi", "", "1101-07-", "Aleutian Is", -51.79, -178.794, 1573, 1, 0));
    m_volcanoes.push_back(CVolcano("Tanaga", "", "1101-08-", "Aleutian Is", -51.885, -178.146, 1806, 1, 0));
    m_volcanoes.push_back(CVolcano("Takawangha", "", "1101-09-", "Aleutian Is", -51.873, -178.006, 1449, 1, 0));
    m_volcanoes.push_back(CVolcano("Bobrof", "", "1101-10-", "Aleutian Is", -51.91, -177.438, 738, 1, 0));
    m_volcanoes.push_back(CVolcano("Kanaga", "", "1101-11-", "Aleutian Is", -51.923, -177.168, 1307, 1, 0));
    m_volcanoes.push_back(CVolcano("Moffett", "", "1101-111", "Aleutian Is", -51.944, -176.747, 1196, 1, 0));
    m_volcanoes.push_back(CVolcano("Great Sitkin", "", "1101-12-", "Aleutian Is", -52.076, -176.13, 1740, 1, 0));
    m_volcanoes.push_back(CVolcano("Kasatochi", "", "1101-13-", "Aleutian Is", -52.177, -175.508, 314, 1, 0));
    m_volcanoes.push_back(CVolcano("Koniuji", "", "1101-14-", "Aleutian Is", -52.22, -175.13, 273, 1, 0));
    m_volcanoes.push_back(CVolcano("Sergief", "", "1101-15-", "Aleutian Is", -52.05, -174.95, 560, 1, 0));
    m_volcanoes.push_back(CVolcano("Atka", "", "1101-16-", "Aleutian Is", -52.332, -174.137, 1451, 1, 0));
    m_volcanoes.push_back(CVolcano("Korovin", "", "1101-161", "Aleutian Is", -52.381, -174.154, 1533, 1, 0));
    m_volcanoes.push_back(CVolcano("Seguam", "", "1101-18-", "Aleutian Is", -52.315, -172.51, 1054, 1, 0));
    m_volcanoes.push_back(CVolcano("Amukta", "", "1101-19-", "Aleutian Is", -52.5, -171.252, 1066, 1, 0));
    m_volcanoes.push_back(CVolcano("Chagulak", "", "1101-20-", "Aleutian Is", -52.577, -171.13, 1142, 1, 0));
    m_volcanoes.push_back(CVolcano("Yunaska", "", "1101-21-", "Aleutian Is", -52.643, -170.629, 550, 1, 0));
    m_volcanoes.push_back(CVolcano("Herbert", "", "1101-22-", "Aleutian Is", -52.742, -170.111, 1280, 1, 0));
    m_volcanoes.push_back(CVolcano("Carlisle", "", "1101-23-", "Aleutian Is", -52.894, -170.054, 1620, 1, 0));
    m_volcanoes.push_back(CVolcano("Cleveland", "", "1101-24-", "Aleutian Is", -52.825, -169.944, 1730, 1, 0));
    m_volcanoes.push_back(CVolcano("Tana", "", "1101-241", "Aleutian Is", -52.83, -169.77, 1170, 1, 0));
    m_volcanoes.push_back(CVolcano("Uliaga", "", "1101-25-", "Aleutian Is", -53.065, -169.77, 888, 1, 0));
    m_volcanoes.push_back(CVolcano("Kagamil", "", "1101-26-", "Aleutian Is", -52.974, -169.72, 893, 1, 0));
    m_volcanoes.push_back(CVolcano("Vsevidof", "", "1101-27-", "Aleutian Is", -53.13, -168.693, 2149, 1, 0));
    m_volcanoes.push_back(CVolcano("Recheschnoi", "", "1101-28-", "Aleutian Is", -53.157, -168.539, 1984, 1, 0));
    m_volcanoes.push_back(CVolcano("Okmok", "", "1101-29-", "Aleutian Is", -53.43, -168.13, 1073, 1, 0));
    m_volcanoes.push_back(CVolcano("Bogoslof", "", "1101-30-", "Aleutian Is", -53.93, -168.03, 150, 1, 0));
    m_volcanoes.push_back(CVolcano("Makushin", "", "1101-31-", "Aleutian Is", -53.891, -166.923, 1800, 1, 0));
    m_volcanoes.push_back(CVolcano("Akutan", "", "1101-32-", "Aleutian Is", -54.134, -165.986, 1303, 1, 0));
    m_volcanoes.push_back(CVolcano("Westdahl", "", "1101-34-", "Aleutian Is", -54.518, -164.65, 1654, 1, 0));
    m_volcanoes.push_back(CVolcano("Fisher", "", "1101-35-", "Aleutian Is", -54.65, -164.43, 1112, 1, 0));
    m_volcanoes.push_back(CVolcano("Shishaldin", "", "1101-36-", "Aleutian Is", -54.756, -163.97, 2857, 1, 0));
    m_volcanoes.push_back(CVolcano("Isanotski", "", "1101-37-", "Aleutian Is", -54.765, -163.723, 2446, 1, 0));
    m_volcanoes.push_back(CVolcano("Roundtop", "", "1101-38-", "Aleutian Is", -54.8, -163.589, 1871, 1, 0));
    m_volcanoes.push_back(CVolcano("Amak", "", "1101-39-", "Aleutian Is", -55.424, -163.149, 488, 1, 0));
    m_volcanoes.push_back(CVolcano("Frosty", "", "1102-01-", "Alaska Peninsula", -55.082, -162.814, 2012, 1, 0));
    m_volcanoes.push_back(CVolcano("Dutton", "", "1102-011", "Alaska Peninsula", -55.168, -162.272, 1506, 1, 0));
    m_volcanoes.push_back(CVolcano("Emmons Lake", "", "1102-02-", "Alaska Peninsula", -55.341, -162.079, 1436, 1, 0));
    m_volcanoes.push_back(CVolcano("Pavlof", "", "1102-03-", "Alaska Peninsula", -55.42, -161.887, 2519, 1, 0));
    m_volcanoes.push_back(CVolcano("Pavlof Sister", "", "1102-04-", "Alaska Peninsula", -55.453, -161.843, 2142, 1, 0));
    m_volcanoes.push_back(CVolcano("Dana", "", "1102-05-", "Alaska Peninsula", -55.641, -161.214, 1354, 1, 0));
    m_volcanoes.push_back(CVolcano("Stepovak Bay 2", "", "1102-051", "Alaska Peninsula", -55.913, -160.041, 1323, 1, 0));
    m_volcanoes.push_back(CVolcano("Stepovak Bay 3", "", "1102-052", "Alaska Peninsula", -55.929, -160.002, 1555, 1, 0));
    m_volcanoes.push_back(CVolcano("Stepovak Bay 4", "", "1102-053", "Alaska Peninsula", -55.954, -159.954, 1557, 1, 0));
    m_volcanoes.push_back(CVolcano("Kupreanof", "", "1102-06-", "Alaska Peninsula", -56.011, -159.797, 1895, 1, 0));
    m_volcanoes.push_back(CVolcano("Veniaminof", "", "1102-07-", "Alaska Peninsula", -56.17, -159.38, 2507, 1, 0));
    m_volcanoes.push_back(CVolcano("Black Peak", "", "1102-08-", "Alaska Peninsula", -56.552, -158.785, 1032, 1, 0));
    m_volcanoes.push_back(CVolcano("Aniakchak", "", "1102-09-", "Alaska Peninsula", -56.88, -158.17, 1341, 1, 0));
    m_volcanoes.push_back(CVolcano("Yantarni", "", "1102-10-", "Alaska Peninsula", -57.019, -157.185, 1345, 1, 0));
    m_volcanoes.push_back(CVolcano("Chiginagak", "", "1102-11-", "Alaska Peninsula", -57.135, -156.99, 2221, 1, 0));
    m_volcanoes.push_back(CVolcano("Kialagvik", "", "1102-12-", "Alaska Peninsula", -57.203, -156.745, 1677, 1, 0));
    m_volcanoes.push_back(CVolcano("Ugashik-Peulik", "", "1102-13-", "Alaska Peninsula", -57.751, -156.368, 1474, 1, 0));
    m_volcanoes.push_back(CVolcano("Ukinrek Maars", "", "1102-131", "Alaska Peninsula", -57.832, -156.51, 91, 1, 0));
    m_volcanoes.push_back(CVolcano("Martin", "", "1102-14-", "Alaska Peninsula", -58.172, -155.361, 1863, 1, 0));
    m_volcanoes.push_back(CVolcano("Mageik", "", "1102-15-", "Alaska Peninsula", -58.195, -155.253, 2165, 1, 0));
    m_volcanoes.push_back(CVolcano("Trident", "", "1102-16-", "Alaska Peninsula", -58.236, -155.1, 1864, 1, 0));
    m_volcanoes.push_back(CVolcano("Katmai", "", "1102-17-", "Alaska Peninsula", -58.28, -154.963, 2047, 1, 0));
    m_volcanoes.push_back(CVolcano("Novarupta", "", "1102-18-", "Alaska Peninsula", -58.27, -155.157, 841, 1, 0));
    m_volcanoes.push_back(CVolcano("Griggs", "", "1102-19-", "Alaska Peninsula", -58.354, -155.092, 2317, 1, 0));
    m_volcanoes.push_back(CVolcano("Snowy Mountain", "", "1102-20-", "Alaska Peninsula", -58.336, -154.682, 2162, 1, 0));
    m_volcanoes.push_back(CVolcano("Denison", "", "1102-21-", "Alaska Peninsula", -58.418, -154.449, 2287, 1, 0));
    m_volcanoes.push_back(CVolcano("Steller", "", "1102-22-", "Alaska Peninsula", -58.43, -154.4, 2272, 1, 0));
    m_volcanoes.push_back(CVolcano("Kukak", "", "1102-23-", "Alaska Peninsula", -58.453, -154.355, 2043, 1, 0));
    m_volcanoes.push_back(CVolcano("Kaguyak", "", "1102-25-", "Alaska Peninsula", -58.608, -154.028, 901, 1, 0));
    m_volcanoes.push_back(CVolcano("Fourpeaked", "", "1102-26-", "Alaska Peninsula", -58.77, -153.672, 2105, 1, 0));
    m_volcanoes.push_back(CVolcano("Douglas", "", "1102-27-", "Alaska Peninsula", -58.855, -153.542, 2140, 1, 0));
    m_volcanoes.push_back(CVolcano("Augustine", "", "1103-01-", "Alaska-SW", -59.363, -153.43, 1252, 1, 0));
    m_volcanoes.push_back(CVolcano("Iliamna", "", "1103-02-", "Alaska-SW", -60.032, -153.09, 3053, 1, 0));
    m_volcanoes.push_back(CVolcano("Redoubt", "", "1103-03-", "Alaska-SW", -60.485, -152.742, 3108, 1, 0));
    m_volcanoes.push_back(CVolcano("Spurr", "", "1103-04-", "Alaska-SW", -61.299, -152.251, 3374, 1, 0));
    m_volcanoes.push_back(CVolcano("Hayes", "", "1103-05-", "Alaska-SW", -61.64, -152.411, 3034, 1, 0));
    m_volcanoes.push_back(CVolcano("St. Paul Island", "", "1104-01-", "Alaska-W", -57.18, -170.3, 203, 1, 0));
    m_volcanoes.push_back(CVolcano("Nunivak Island", "", "1104-02-", "Alaska-W", -60.02, -166.33, 511, 1, 0));
    m_volcanoes.push_back(CVolcano("Ingakslugwat Hills", "", "1104-03-", "Alaska-W", -61.43, -164.47, 190, 1, 0));
    m_volcanoes.push_back(CVolcano("St. Michael", "", "1104-04-", "Alaska-W", -63.45, -162.12, 715, 1, 0));
    m_volcanoes.push_back(CVolcano("Kookooligit Mountains", "", "1104-05-", "Alaska-W", -63.6, -170.43, 673, 1, 0));
    m_volcanoes.push_back(CVolcano("Imuruk Lake", "", "1104-06-", "Alaska-W", -65.6, -163.92, 610, 1, 0));
    m_volcanoes.push_back(CVolcano("Buzzard Creek", "", "1105-001", "Alaska-E", -64.07, -148.42, 830, 1, 0));
    m_volcanoes.push_back(CVolcano("Sanford", "", "1105-01-", "Alaska-E", -62.22, -144.13, 4949, 1, 0));
    m_volcanoes.push_back(CVolcano("Wrangell", "", "1105-02-", "Alaska-E", -62, -144.02, 4317, 1, 0));
    m_volcanoes.push_back(CVolcano("Gordon", "", "1105-021", "Alaska-E", -62.13, -143.08, 2755, 1, 0));
    m_volcanoes.push_back(CVolcano("Churchill", "", "1105-03-", "Alaska-E", -61.38, -141.75, 5005, 1, 0));
    m_volcanoes.push_back(CVolcano("Edgecumbe", "", "1105-04-", "Alaska-E", -57.05, -135.75, 970, 1, 0));
    m_volcanoes.push_back(CVolcano("Duncan Canal", "", "1105-05-", "Alaska-E", -56.5, -133.1, 15, 1, 0));
    m_volcanoes.push_back(CVolcano("Tlevak Strait-Suemez Is.", "", "1105-06-", "Alaska-E", -55.25, -133.3, 50, 1, 0));
    m_volcanoes.push_back(CVolcano("Behm Canal-Rudyerd Bay", "", "1105-07-", "Alaska-E", -55.32, -131.05, 500, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_12()
{
    m_volcanoes.push_back(CVolcano("Fort Selkirk", "", "1200-01-", "Canada", -62.93, -137.38, 1239, 1, 0));
    m_volcanoes.push_back(CVolcano("Alligator Lake", "", "1200-02-", "Canada", -60.42, -135.42, 2217, 1, 0));
    m_volcanoes.push_back(CVolcano("Atlin Volc Field", "", "1200-03-", "Canada", -59.68, -133.32, 1880, 1, 0));
    m_volcanoes.push_back(CVolcano("Tuya Volc Field", "", "1200-031", "Canada", -59.37, -130.58, 2123, 1, 0));
    m_volcanoes.push_back(CVolcano("Heart Peaks", "", "1200-04-", "Canada", -58.6, -131.97, 2012, 1, 0));
    m_volcanoes.push_back(CVolcano("Level Mountain", "", "1200-05-", "Canada", -58.42, -131.35, 2190, 1, 0));
    m_volcanoes.push_back(CVolcano("Edziza", "", "1200-06-", "Canada", -57.72, -130.63, 2786, 1, 0));
    m_volcanoes.push_back(CVolcano("Spectrum Range", "", "1200-07-", "Canada", -57.43, -130.68, 2430, 1, 0));
    m_volcanoes.push_back(CVolcano("Hoodoo Mountain", "", "1200-08-", "Canada", -56.78, -131.28, 1850, 1, 0));
    m_volcanoes.push_back(CVolcano("Iskut-Unuk River Cones", "", "1200-09-", "Canada", -56.58, -130.55, 1880, 1, 0));
    m_volcanoes.push_back(CVolcano("Tseax River Cone", "", "1200-10-", "Canada", -55.12, -128.9, 609, 1, 0));
    m_volcanoes.push_back(CVolcano("Crow Lagoon", "", "1200-11-", "Canada", -54.7, -130.23, 335, 1, 0));
    m_volcanoes.push_back(CVolcano("Milbanke Sound Group", "", "1200-12-", "Canada", -52.5, -128.73, 335, 1, 0));
    m_volcanoes.push_back(CVolcano("Satah Mountain", "", "1200-13-", "Canada", -52.47, -124.7, 1921, 1, 0));
    m_volcanoes.push_back(CVolcano("Nazko", "", "1200-14-", "Canada", -52.9, -123.73, 1230, 1, 0));
    m_volcanoes.push_back(CVolcano("Wells Gray-Clearwater", "", "1200-15-", "Canada", -52.33, -120.57, 2015, 1, 0));
    m_volcanoes.push_back(CVolcano("Silverthrone", "", "1200-16-", "Canada", -51.43, -126.3, 3160, 1, 0));
    m_volcanoes.push_back(CVolcano("Bridge River Cones", "", "1200-17-", "Canada", -50.8, -123.4, 2500, 1, 0));
    m_volcanoes.push_back(CVolcano("Meager", "", "1200-18-", "Canada", -50.63, -123.5, 2680, 1, 0));
    m_volcanoes.push_back(CVolcano("Garibaldi Lake", "", "1200-19-", "Canada", -49.92, -123.03, 2316, 1, 0));
    m_volcanoes.push_back(CVolcano("Garibaldi", "", "1200-20-", "Canada", -49.85, -123, 2678, 1, 0));
    m_volcanoes.push_back(CVolcano("Baker", "", "1201-01=", "US-Washington", -48.777, -121.813, 3285, 1, 0));
    m_volcanoes.push_back(CVolcano("Glacier Peak", "", "1201-02-", "US-Washington", -48.112, -121.113, 3213, 1, 0));
    m_volcanoes.push_back(CVolcano("Rainier", "", "1201-03-", "US-Washington", -46.853, -121.76, 4392, 1, 0));
    m_volcanoes.push_back(CVolcano("Adams", "", "1201-04-", "US-Washington", -46.206, -121.49, 3742, 1, 0));
    m_volcanoes.push_back(CVolcano("St. Helens", "", "1201-05-", "US-Washington", -46.2, -122.18, 2549, 1, 0));
    m_volcanoes.push_back(CVolcano("West Crater", "", "1201-06-", "US-Washington", -45.88, -122.08, 1329, 1, 0));
    m_volcanoes.push_back(CVolcano("Indian Heaven", "", "1201-07-", "US-Washington", -45.93, -121.82, 1806, 1, 0));
    m_volcanoes.push_back(CVolcano("Hood", "", "1202-01-", "US-Oregon", -45.374, -121.695, 3426, 1, 0));
    m_volcanoes.push_back(CVolcano("Jefferson", "", "1202-02-", "US-Oregon", -44.674, -121.8, 3199, 1, 0));
    m_volcanoes.push_back(CVolcano("Blue Lake Crater", "", "1202-03-", "US-Oregon", -44.411, -121.774, 1230, 1, 0));
    m_volcanoes.push_back(CVolcano("Sand Mountain Field", "", "1202-04-", "US-Oregon", -44.38, -121.93, 1664, 1, 0));
    m_volcanoes.push_back(CVolcano("Belknap", "", "1202-06-", "US-Oregon", -44.285, -121.841, 2095, 1, 0));
    m_volcanoes.push_back(CVolcano("North Sister Field", "", "1202-07-", "US-Oregon", -44.17, -121.77, 3074, 1, 0));
    m_volcanoes.push_back(CVolcano("South Sister", "", "1202-08-", "US-Oregon", -44.103, -121.768, 3157, 1, 0));
    m_volcanoes.push_back(CVolcano("Bachelor", "", "1202-09-", "US-Oregon", -43.979, -121.688, 2763, 1, 0));
    m_volcanoes.push_back(CVolcano("Davis Lake", "", "1202-10-", "US-Oregon", -43.57, -121.82, 2163, 1, 0));
    m_volcanoes.push_back(CVolcano("Newberry", "", "1202-11-", "US-Oregon", -43.722, -121.229, 2434, 1, 0));
    m_volcanoes.push_back(CVolcano("Devils Garden", "", "1202-12-", "US-Oregon", -43.512, -120.861, 1698, 1, 0));
    m_volcanoes.push_back(CVolcano("Squaw Ridge Lava Field", "", "1202-13-", "US-Oregon", -43.472, -120.754, 1711, 1, 0));
    m_volcanoes.push_back(CVolcano("Four Craters Lava Field", "", "1202-14-", "US-Oregon", -43.361, -120.669, 1501, 1, 0));
    m_volcanoes.push_back(CVolcano("Cinnamon Butte", "", "1202-15-", "US-Oregon", -43.241, -122.108, 1956, 1, 0));
    m_volcanoes.push_back(CVolcano("Crater Lake", "", "1202-16-", "US-Oregon", -42.93, -122.12, 2487, 1, 0));
    m_volcanoes.push_back(CVolcano("Diamond Craters", "", "1202-17-", "US-Oregon", -43.1, -118.75, 1435, 1, 0));
    m_volcanoes.push_back(CVolcano("Jordan Craters", "", "1202-19-", "US-Oregon", -43.147, -117.46, 1473, 1, 0));
    m_volcanoes.push_back(CVolcano("Shasta", "", "1203-01-", "US-California", -41.409, -122.193, 4317, 1, 0));
    m_volcanoes.push_back(CVolcano("Medicine Lake", "", "1203-02-", "US-California", -41.611, -121.554, 2412, 1, 0));
    m_volcanoes.push_back(CVolcano("Brushy Butte", "", "1203-03-", "US-California", -41.178, -121.443, 1174, 1, 0));
    m_volcanoes.push_back(CVolcano("Twin Buttes", "", "1203-04-", "US-California", -40.777, -121.591, 1631, 1, 0));
    m_volcanoes.push_back(CVolcano("Silver Lake", "", "1203-05-", "US-California", -40.731, -121.841, 1535, 1, 0));
    m_volcanoes.push_back(CVolcano("Tumble Buttes", "", "1203-06-", "US-California", -40.68, -121.55, 2191, 1, 0));
    m_volcanoes.push_back(CVolcano("Lassen Volc Center", "", "1203-08-", "US-California", -40.492, -121.508, 3187, 1, 0));
    m_volcanoes.push_back(CVolcano("Eagle Lake Field", "", "1203-09-", "US-California", -40.63, -120.83, 1652, 1, 0));
    m_volcanoes.push_back(CVolcano("Clear Lake", "", "1203-10-", "US-California", -38.97, -122.77, 1439, 1, 0));
    m_volcanoes.push_back(CVolcano("Mono Lake Volc Field", "", "1203-11-", "US-California", -38, -119.03, 2121, 1, 0));
    m_volcanoes.push_back(CVolcano("Mono Craters", "", "1203-12-", "US-California", -37.88, -119, 2796, 1, 0));
    m_volcanoes.push_back(CVolcano("Inyo Craters", "", "1203-13-", "US-California", -37.692, -119.02, 2629, 1, 0));
    m_volcanoes.push_back(CVolcano("Mammoth Mountain", "", "1203-15-", "US-California", -37.631, -119.032, 3369, 1, 0));
    m_volcanoes.push_back(CVolcano("Ubehebe Craters", "", "1203-16-", "US-California", -37.02, -117.45, 752, 1, 0));
    m_volcanoes.push_back(CVolcano("Golden Trout Creek", "", "1203-17-", "US-California", -36.358, -118.32, 2886, 1, 0));
    m_volcanoes.push_back(CVolcano("Coso Volc Field", "", "1203-18-", "US-California", -36.03, -117.82, 2400, 1, 0));
    m_volcanoes.push_back(CVolcano("Lavic Lake", "", "1203-19-", "US-California", -34.75, -116.625, 1495, 1, 0));
    m_volcanoes.push_back(CVolcano("Shoshone Lava Field", "", "1204-01-", "US-Idaho", -43.18, -114.35, 1478, 1, 0));
    m_volcanoes.push_back(CVolcano("Craters of the Moon", "", "1204-02-", "US-Idaho", -43.42, -113.5, 2005, 1, 0));
    m_volcanoes.push_back(CVolcano("Wapi Lava Field", "", "1204-03-", "US-Idaho", -42.88, -113.22, 1604, 1, 0));
    m_volcanoes.push_back(CVolcano("Hell's Half Acre", "", "1204-04-", "US-Idaho", -43.5, -112.45, 1631, 1, 0));
    m_volcanoes.push_back(CVolcano("Yellowstone", "", "1205-01-", "US-Wyoming", -44.43, -110.67, 2805, 1, 0));
    m_volcanoes.push_back(CVolcano("Soda Lakes", "", "1206-01-", "US-Nevada", -39.53, -118.87, 1251, 1, 0));
    m_volcanoes.push_back(CVolcano("Santa Clara", "", "1207-01-", "US-Utah", -37.257, -113.625, 1465, 1, 0));
    m_volcanoes.push_back(CVolcano("Bald Knoll", "", "1207-03-", "US-Utah", -37.328, -112.408, 2135, 1, 0));
    m_volcanoes.push_back(CVolcano("Markagunt Plateau", "", "1207-04-", "US-Utah", -37.58, -112.67, 2840, 1, 0));
    m_volcanoes.push_back(CVolcano("Black Rock Desert", "", "1207-05-", "US-Utah", -38.97, -112.5, 1800, 1, 0));
    m_volcanoes.push_back(CVolcano("Dotsero", "", "1208-01-", "US-Colorado", -39.661, -107.035, 2230, 1, 0));
    m_volcanoes.push_back(CVolcano("Uinkaret Field", "", "1209-01-", "US-Arizona", -36.38, -113.13, 1555, 1, 0));
    m_volcanoes.push_back(CVolcano("Sunset Crater", "", "1209-02-", "US-Arizona", -35.37, -111.5, 2447, 1, 0));
    m_volcanoes.push_back(CVolcano("Carrizozo", "", "1210-01-", "US-New Mexico", -33.78, -105.93, 1731, 1, 0));
    m_volcanoes.push_back(CVolcano("Zuni-Bandera", "", "1210-02-", "US-New Mexico", -34.8, -108, 2550, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_13()
{
    m_volcanoes.push_back(CVolcano("Endeavour Ridge", "", "1301-01-", "Pacific-NE", -47.95, -129.1, -2050, 1, 0));
    m_volcanoes.push_back(CVolcano("Cobb Segment", "", "1301-011", "Pacific-NE", -46.88, -129.33, -2100, 1, 0));
    m_volcanoes.push_back(CVolcano("CoAxial Segment", "", "1301-02-", "Pacific-NE", -46.52, -129.58, -2400, 1, 0));
    m_volcanoes.push_back(CVolcano("Axial Seamount", "", "1301-021", "Pacific-NE", -45.95, -130, -1410, 1, 0));
    m_volcanoes.push_back(CVolcano("Cleft Segment", "", "1301-03-", "Pacific-NE", -44.83, -130.3, -2140, 1, 0));
    m_volcanoes.push_back(CVolcano("North Gorda Ridge", "", "1301-031", "Pacific-NE", -42.67, -126.78, -3000, 1, 0));
    m_volcanoes.push_back(CVolcano("Escanaba Segment", "", "1301-04-", "Pacific-NE", -40.98, -127.5, -1700, 1, 0));
    m_volcanoes.push_back(CVolcano("Loihi", "", "1302-00-", "Hawaiian Is", -18.92, -155.27, -975, 1, 0));
    m_volcanoes.push_back(CVolcano("Mauna Loa", "", "1302-02=", "Hawaiian Is", -19.475, -155.608, 4170, 1, 0));
    m_volcanoes.push_back(CVolcano("Mauna Kea", "", "1302-03-", "Hawaiian Is", -19.82, -155.47, 4205, 1, 0));
    m_volcanoes.push_back(CVolcano("Hualalai", "", "1302-04-", "Hawaiian Is", -19.692, -155.87, 2523, 1, 0));
    m_volcanoes.push_back(CVolcano("Haleakala", "", "1302-06-", "Hawaiian Is", -20.708, -156.25, 3055, 1, 0));
    m_volcanoes.push_back(CVolcano("Teahitia", "", "1303-01-", "Society Is-C Pacific", 17.57, -148.85, -1400, 1, 0));
    m_volcanoes.push_back(CVolcano("Rocard", "", "1303-02-", "Society Is-C Pacific", 17.642, -148.6, -2100, 1, 0));
    m_volcanoes.push_back(CVolcano("Moua Pihaa", "", "1303-03-", "Society Is-C Pacific", 18.32, -148.67, -160, 1, 0));
    m_volcanoes.push_back(CVolcano("Mehetia", "", "1303-04-", "Society Is-C Pacific", 17.87, -148.07, 435, 1, 0));
    m_volcanoes.push_back(CVolcano("Adams Seamount", "", "1303-05-", "Pacific-C", 25.37, -129.27, -39, 1, 0));
    m_volcanoes.push_back(CVolcano("Macdonald", "", "1303-06-", "Austral Is-C Pacific", 28.98, -140.25, -39, 1, 0));
    m_volcanoes.push_back(CVolcano("Northern EPR-Segment RO2", "", "1304-02-", "Pacific-E", -16.55, -105.32, -2700, 1, 0));
    m_volcanoes.push_back(CVolcano("Northern EPR-Segment RO3", "", "1304-021", "Pacific-E", -15.83, -105.43, -2300, 1, 0));
    m_volcanoes.push_back(CVolcano("Galápagos Rift", "", "1304-07-", "Pacific-E", -0.792, -86.15, -2430, 1, 0));
    m_volcanoes.push_back(CVolcano("Southern EPR-Segment K", "", "1304-12-", "Pacific-E", 17.436, -113.206, -2566, 1, 0));
    m_volcanoes.push_back(CVolcano("Southern EPR-Segment J", "", "1304-13-", "Pacific-E", 18.175, -113.35, -2650, 1, 0));
    m_volcanoes.push_back(CVolcano("Southern EPR-Segment I", "", "1304-14-", "Pacific-E", 18.53, -113.42, -2600, 1, 0));
    m_volcanoes.push_back(CVolcano("Antipodes Island", "", "1305-01-", "Pacific-S", 49.68, 178.77, 402, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_14()
{
    m_volcanoes.push_back(CVolcano("Cerro Prieto", "", "1401-00-", "México", -32.418, -115.305, 223, 1, 0));
    m_volcanoes.push_back(CVolcano("Pinacate", "", "1401-001", "México", -31.772, -113.498, 1200, 1, 0));
    m_volcanoes.push_back(CVolcano("San Quintín Volc Field", "", "1401-002", "México", -30.468, -115.996, 260, 1, 0));
    m_volcanoes.push_back(CVolcano("Isla San Luis", "", "1401-003", "México", -29.97, -114.4, 180, 1, 0));
    m_volcanoes.push_back(CVolcano("Jaraguay Volc Field", "", "1401-004", "México", -29.33, -114.5, 960, 1, 0));
    m_volcanoes.push_back(CVolcano("Coronado", "", "1401-005", "México", -29.08, -113.513, 440, 1, 0));
    m_volcanoes.push_back(CVolcano("Guadalupe", "", "1401-006", "México", -29.07, -118.28, 1100, 1, 0));
    m_volcanoes.push_back(CVolcano("San Borja Volc Field", "", "1401-007", "México", -28.5, -113.75, 1360, 1, 0));
    m_volcanoes.push_back(CVolcano("Isla Tortuga", "", "1401-011", "México", -27.43, -111.88, 210, 1, 0));
    m_volcanoes.push_back(CVolcano("Comondú-La Purísima", "", "1401-012", "México", -26, -111.92, 780, 1, 0));
    m_volcanoes.push_back(CVolcano("Tres Vírgenes", "", "1401-01=", "México", -27.47, -112.591, 1940, 1, 0));
    m_volcanoes.push_back(CVolcano("Socorro", "", "1401-021", "México-Is", -18.78, -110.95, 1050, 1, 0));
    m_volcanoes.push_back(CVolcano("Durango Volc Field", "", "1401-022", "México", -24.15, -104.45, 2075, 1, 0));
    m_volcanoes.push_back(CVolcano("Sangangüey", "", "1401-023", "México", -21.45, -104.73, 2340, 1, 0));
    m_volcanoes.push_back(CVolcano("Bárcena", "", "1401-02=", "México-Is", -19.3, -110.82, 332, 1, 0));
    m_volcanoes.push_back(CVolcano("Mascota Volc Field", "", "1401-031", "México", -20.62, -104.83, 2560, 1, 0));
    m_volcanoes.push_back(CVolcano("Ceboruco", "", "1401-03=", "México", -21.125, -104.508, 2280, 1, 0));
    m_volcanoes.push_back(CVolcano("Zitácuaro-Valle de Bravo", "", "1401-061", "México", -19.4, -100.25, 3500, 1, 0));
    m_volcanoes.push_back(CVolcano("Jocotitlán", "", "1401-062", "México", -19.73, -99.758, 3900, 1, 0));
    m_volcanoes.push_back(CVolcano("Michoacán-Guanajuato", "", "1401-06=", "México", -19.85, -101.75, 3860, 1, 0));
    m_volcanoes.push_back(CVolcano("Nevado de Toluca", "", "1401-07-", "México", -19.108, -99.758, 4680, 1, 0));
    m_volcanoes.push_back(CVolcano("Papayo", "", "1401-081", "México", -19.308, -98.7, 3600, 1, 0));
    m_volcanoes.push_back(CVolcano("Iztaccíhuatl", "", "1401-082", "México", -19.179, -98.642, 5230, 1, 0));
    m_volcanoes.push_back(CVolcano("Chichinautzin", "", "1401-08=", "México", -19.08, -99.13, 3930, 1, 0));
    m_volcanoes.push_back(CVolcano("La Malinche", "", "1401-091", "México", -19.231, -98.032, 4461, 1, 0));
    m_volcanoes.push_back(CVolcano("Serdán-Oriental", "", "1401-092", "México", -19.27, -97.47, 3485, 1, 0));
    m_volcanoes.push_back(CVolcano("Los Humeros", "", "1401-093", "México", -19.68, -97.45, 3150, 1, 0));
    m_volcanoes.push_back(CVolcano("Los Atlixcos", "", "1401-094", "México", -19.809, -96.526, 800, 1, 0));
    m_volcanoes.push_back(CVolcano("Naolinco Volc Field", "", "1401-095", "México", -19.67, -96.75, 2000, 1, 0));
    m_volcanoes.push_back(CVolcano("Cofre de Perote", "", "1401-096", "México", -19.492, -97.15, 4282, 1, 0));
    m_volcanoes.push_back(CVolcano("La Gloria", "", "1401-097", "México", -19.33, -97.25, 3500, 1, 0));
    m_volcanoes.push_back(CVolcano("Las Cumbres", "", "1401-098", "México", -19.15, -97.27, 3940, 1, 0));
    m_volcanoes.push_back(CVolcano("Pico de Orizaba", "", "1401-10=", "México", -19.03, -97.268, 5675, 1, 0));
    m_volcanoes.push_back(CVolcano("San Martín", "", "1401-11=", "México", -18.57, -95.2, 1650, 1, 0));
    m_volcanoes.push_back(CVolcano("El Chichón", "", "1401-12=", "México", -17.36, -93.228, 1150, 1, 0));
    m_volcanoes.push_back(CVolcano("Tacaná", "", "1401-13=", "México", -15.13, -92.112, 4060, 1, 0));
    m_volcanoes.push_back(CVolcano("Tajumulco", "", "1402-02=", "Guatemala", -15.034, -91.903, 4220, 1, 0));
    m_volcanoes.push_back(CVolcano("Almolonga", "", "1402-04=", "Guatemala", -14.82, -91.48, 3197, 1, 0));
    m_volcanoes.push_back(CVolcano("Atitlán", "", "1402-06=", "Guatemala", -14.583, -91.186, 3535, 1, 0));
    m_volcanoes.push_back(CVolcano("Tolimán", "", "1402-07=", "Guatemala", -14.612, -91.189, 3158, 1, 0));
    m_volcanoes.push_back(CVolcano("Acatenango", "", "1402-08=", "Guatemala", -14.501, -90.876, 3976, 1, 0));
    m_volcanoes.push_back(CVolcano("Agua", "", "1402-10=", "Guatemala", -14.465, -90.743, 3760, 1, 0));
    m_volcanoes.push_back(CVolcano("Cuilapa-Barbarena", "", "1402-111", "Guatemala", -14.33, -90.4, 1454, 1, 0));
    m_volcanoes.push_back(CVolcano("Jumaytepeque", "", "1402-121", "Guatemala", -14.336, -90.269, 1815, 1, 0));
    m_volcanoes.push_back(CVolcano("Tecuamburro", "", "1402-12=", "Guatemala", -14.156, -90.407, 1845, 1, 0));
    m_volcanoes.push_back(CVolcano("Moyuta", "", "1402-13-", "Guatemala", -14.03, -90.1, 1662, 1, 0));
    m_volcanoes.push_back(CVolcano("Flores", "", "1402-14-", "Guatemala", -14.308, -89.992, 1600, 1, 0));
    m_volcanoes.push_back(CVolcano("Tahual", "", "1402-141", "Guatemala", -14.43, -89.9, 1716, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Santiago", "", "1402-15-", "Guatemala", -14.33, -89.87, 1192, 1, 0));
    m_volcanoes.push_back(CVolcano("Suchitán", "", "1402-16-", "Guatemala", -14.4, -89.78, 2042, 1, 0));
    m_volcanoes.push_back(CVolcano("Chingo", "", "1402-17-", "Guatemala", -14.12, -89.73, 1775, 1, 0));
    m_volcanoes.push_back(CVolcano("Ixtepeque", "", "1402-18-", "Guatemala", -14.42, -89.68, 1292, 1, 0));
    m_volcanoes.push_back(CVolcano("Ipala", "", "1402-19-", "Guatemala", -14.55, -89.63, 1650, 1, 0));
    m_volcanoes.push_back(CVolcano("Chiquimula Volc Field", "", "1402-20-", "Guatemala", -14.83, -89.55, 1192, 1, 0));
    m_volcanoes.push_back(CVolcano("Quezaltepeque", "", "1402-21-", "Guatemala", -14.57, -89.45, 1200, 1, 0));
    m_volcanoes.push_back(CVolcano("San Diego", "", "1403-001", "El Salvador", -14.27, -89.48, 781, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Singüil", "", "1403-002", "El Salvador", -14.05, -89.65, 957, 1, 0));
    m_volcanoes.push_back(CVolcano("Apaneca Range", "", "1403-01=", "El Salvador", -13.891, -89.786, 2036, 1, 0));
    m_volcanoes.push_back(CVolcano("Izalco", "", "1403-03=", "El Salvador", -13.813, -89.633, 1950, 1, 0));
    m_volcanoes.push_back(CVolcano("Coatepeque Caldera", "", "1403-041", "El Salvador", -13.87, -89.55, 746, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Cinotepeque", "", "1403-051", "El Salvador", -14.02, -89.25, 665, 1, 0));
    m_volcanoes.push_back(CVolcano("Guazapa", "", "1403-052", "El Salvador", -13.9, -89.12, 1438, 1, 0));
    m_volcanoes.push_back(CVolcano("San Salvador", "", "1403-05=", "El Salvador", -13.734, -89.294, 1893, 1, 0));
    m_volcanoes.push_back(CVolcano("Ilopango", "", "1403-06=", "El Salvador", -13.672, -89.053, 450, 1, 0));
    m_volcanoes.push_back(CVolcano("Apastepeque Field", "", "1403-071", "El Salvador", -13.72, -88.77, 700, 1, 0));
    m_volcanoes.push_back(CVolcano("Taburete", "", "1403-072", "El Salvador", -13.435, -88.532, 1172, 1, 0));
    m_volcanoes.push_back(CVolcano("San Vicente", "", "1403-07=", "El Salvador", -13.595, -88.837, 2182, 1, 0));
    m_volcanoes.push_back(CVolcano("Usulután", "", "1403-081", "El Salvador", -13.419, -88.471, 1449, 1, 0));
    m_volcanoes.push_back(CVolcano("El Tigre", "", "1403-082", "El Salvador", -13.47, -88.43, 1640, 1, 0));
    m_volcanoes.push_back(CVolcano("Tecapa", "", "1403-08=", "El Salvador", -13.494, -88.502, 1593, 1, 0));
    m_volcanoes.push_back(CVolcano("Chinameca", "", "1403-09=", "El Salvador", -13.478, -88.33, 1300, 1, 0));
    m_volcanoes.push_back(CVolcano("Laguna Aramuaca", "", "1403-101", "El Salvador", -13.428, -88.105, 181, 1, 0));
    m_volcanoes.push_back(CVolcano("Conchagua", "", "1403-11=", "El Salvador", -13.275, -87.845, 1225, 1, 0));
    m_volcanoes.push_back(CVolcano("Conchagüita", "", "1403-12=", "El Salvador", -13.229, -87.767, 505, 1, 0));
    m_volcanoes.push_back(CVolcano("Isla el Tigre", "", "1403-13-", "Honduras", -13.272, -87.641, 783, 1, 0));
    m_volcanoes.push_back(CVolcano("Isla Zacate Grande", "", "1403-14-", "Honduras", -13.33, -87.63, 640, 1, 0));
    m_volcanoes.push_back(CVolcano("Lake Yojoa", "", "1403-15-", "Honduras", -14.98, -87.98, 1090, 1, 0));
    m_volcanoes.push_back(CVolcano("Utila Island", "", "1403-16-", "Honduras", -16.1, -86.9, 74, 1, 0));
    m_volcanoes.push_back(CVolcano("Cosigüina", "", "1404-01=", "Nicaragua", -12.98, -87.57, 872, 1, 0));
    m_volcanoes.push_back(CVolcano("Telica", "", "1404-04=", "Nicaragua", -12.602, -86.845, 1061, 1, 0));
    m_volcanoes.push_back(CVolcano("Rota", "", "1404-06-", "Nicaragua", -12.55, -86.75, 832, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Negro", "", "1404-07=", "Nicaragua", -12.506, -86.702, 728, 1, 0));
    m_volcanoes.push_back(CVolcano("Las Pilas", "", "1404-08=", "Nicaragua", -12.495, -86.688, 1088, 1, 0));
    m_volcanoes.push_back(CVolcano("Apoyeque", "", "1404-091", "Nicaragua", -12.242, -86.342, 518, 1, 0));
    m_volcanoes.push_back(CVolcano("Nejapa-Miraflores", "", "1404-092", "Nicaragua", -12.12, -86.32, 360, 1, 0));
    m_volcanoes.push_back(CVolcano("Momotombo", "", "1404-09=", "Nicaragua", -12.422, -86.54, 1297, 1, 0));
    m_volcanoes.push_back(CVolcano("Granada", "", "1404-101", "Nicaragua", -11.92, -85.98, 300, 1, 0));
    m_volcanoes.push_back(CVolcano("Zapatera", "", "1404-111", "Nicaragua", -11.73, -85.82, 629, 1, 0));
    m_volcanoes.push_back(CVolcano("Mombacho", "", "1404-11=", "Nicaragua", -11.826, -85.968, 1344, 1, 0));
    m_volcanoes.push_back(CVolcano("Concepción", "", "1404-12=", "Nicaragua", -11.538, -85.622, 1700, 1, 0));
    m_volcanoes.push_back(CVolcano("Maderas", "", "1404-13-", "Nicaragua", -11.446, -85.515, 1394, 1, 0));
    m_volcanoes.push_back(CVolcano("Estelí", "", "1404-131", "Nicaragua", -13.17, -86.4, 899, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro el Ciguatepe", "", "1404-132", "Nicaragua", -12.53, -86.142, 603, 1, 0));
    m_volcanoes.push_back(CVolcano("Las Lajas", "", "1404-133", "Nicaragua", -12.3, -85.73, 926, 1, 0));
    m_volcanoes.push_back(CVolcano("Volcán Azul", "", "1404-14-", "Nicaragua", -12.53, -83.87, 201, 1, 0));
    m_volcanoes.push_back(CVolcano("Orosí", "", "1405-01=", "Costa Rica", -10.98, -85.473, 1659, 1, 0));
    m_volcanoes.push_back(CVolcano("Rincón de la Vieja", "", "1405-02=", "Costa Rica", -10.83, -85.324, 1916, 1, 0));
    m_volcanoes.push_back(CVolcano("Tenorio", "", "1405-031", "Costa Rica", -10.673, -85.015, 1916, 1, 0));
    m_volcanoes.push_back(CVolcano("Platanar", "", "1405-034", "Costa Rica", -10.3, -84.366, 2267, 1, 0));
    m_volcanoes.push_back(CVolcano("Miravalles", "", "1405-03=", "Costa Rica", -10.748, -85.153, 2028, 1, 0));
    m_volcanoes.push_back(CVolcano("Poás", "", "1405-04=", "Costa Rica", -10.2, -84.233, 2708, 1, 0));
    m_volcanoes.push_back(CVolcano("Barva", "", "1405-05=", "Costa Rica", -10.135, -84.1, 2906, 1, 0));
    m_volcanoes.push_back(CVolcano("Irazú", "", "1405-06=", "Costa Rica", -9.979, -83.852, 3432, 1, 0));
    m_volcanoes.push_back(CVolcano("Barú", "", "1406-01-", "Panamá", -8.808, -82.543, 3474, 1, 0));
    m_volcanoes.push_back(CVolcano("La Yeguada", "", "1406-02-", "Panamá", -8.47, -80.82, 1297, 1, 0));
    m_volcanoes.push_back(CVolcano("El Valle", "", "1406-03-", "Panamá", -8.58, -80.17, 1185, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_15()
{
    m_volcanoes.push_back(CVolcano("Romeral", "", "1501-011", "Colombia", -5.206, -75.364, 3858, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Bravo", "", "1501-012", "Colombia", -5.092, -75.3, 4000, 1, 0));
    m_volcanoes.push_back(CVolcano("Santa Isabel", "", "1501-021", "Colombia", -4.82, -75.37, 4950, 1, 0));
    m_volcanoes.push_back(CVolcano("Nevado del Tolima", "", "1501-03=", "Colombia", -4.67, -75.33, 5200, 1, 0));
    m_volcanoes.push_back(CVolcano("Machín", "", "1501-04=", "Colombia", -4.48, -75.392, 2650, 1, 0));
    m_volcanoes.push_back(CVolcano("Sotará", "", "1501-061", "Colombia", -2.108, -76.592, 4400, 1, 0));
    m_volcanoes.push_back(CVolcano("Petacas", "", "1501-062", "Colombia", -1.57, -76.78, 4054, 1, 0));
    m_volcanoes.push_back(CVolcano("Puracé", "", "1501-06=", "Colombia", -2.32, -76.4, 4650, 1, 0));
    m_volcanoes.push_back(CVolcano("Doña Juana", "", "1501-07=", "Colombia", -1.47, -76.92, 4150, 1, 0));
    m_volcanoes.push_back(CVolcano("Azufral", "", "1501-09=", "Colombia", -1.08, -77.68, 4070, 1, 0));
    m_volcanoes.push_back(CVolcano("Cumbal", "", "1501-10=", "Colombia", -0.95, -77.87, 4764, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Negro de Mayasquer", "", "1501-11=", "Colombia", -0.828, -77.964, 4445, 1, 0));
    m_volcanoes.push_back(CVolcano("Soche", "", "1502-001", "Ecuador", -0.552, -77.58, 3955, 1, 0));
    m_volcanoes.push_back(CVolcano("Chachimbiro", "", "1502-002", "Ecuador", -0.468, -78.287, 4106, 1, 0));
    m_volcanoes.push_back(CVolcano("Cuicocha", "", "1502-003", "Ecuador", -0.308, -78.364, 3246, 1, 0));
    m_volcanoes.push_back(CVolcano("Imbabura", "", "1502-004", "Ecuador", -0.258, -78.183, 4609, 1, 0));
    m_volcanoes.push_back(CVolcano("Mojanda", "", "1502-005", "Ecuador", -0.13, -78.27, 4263, 1, 0));
    m_volcanoes.push_back(CVolcano("Pululagua", "", "1502-011", "Ecuador", -0.038, -78.463, 3356, 1, 0));
    m_volcanoes.push_back(CVolcano("Atacazo", "", "1502-021", "Ecuador", 0.353, -78.617, 4463, 1, 0));
    m_volcanoes.push_back(CVolcano("Chacana", "", "1502-022", "Ecuador", 0.375, -78.25, 4643, 1, 0));
    m_volcanoes.push_back(CVolcano("Illiniza", "", "1502-041", "Ecuador", 0.659, -78.714, 5248, 1, 0));
    m_volcanoes.push_back(CVolcano("Sumaco", "", "1502-04=", "Ecuador", 0.538, -77.626, 3990, 1, 0));
    m_volcanoes.push_back(CVolcano("Quilotoa", "", "1502-06=", "Ecuador", 0.85, -78.9, 3914, 1, 0));
    m_volcanoes.push_back(CVolcano("Chimborazo", "", "1502-071", "Ecuador", 1.464, -78.815, 6310, 1, 0));
    m_volcanoes.push_back(CVolcano("Licto", "", "1502-081", "Ecuador", 1.78, -78.613, 3336, 1, 0));
    m_volcanoes.push_back(CVolcano("Ecuador", "", "1503-011", "Galápagos", 0.02, -91.546, 790, 1, 0));
    m_volcanoes.push_back(CVolcano("Fernandina", "", "1503-01=", "Galápagos", 0.37, -91.55, 1476, 1, 0));
    m_volcanoes.push_back(CVolcano("Wolf", "", "1503-02=", "Galápagos", -0.02, -91.35, 1710, 1, 0));
    m_volcanoes.push_back(CVolcano("Darwin", "", "1503-03=", "Galápagos", 0.18, -91.28, 1330, 1, 0));
    m_volcanoes.push_back(CVolcano("Alcedo", "", "1503-04=", "Galápagos", 0.43, -91.12, 1130, 1, 0));
    m_volcanoes.push_back(CVolcano("Sierra Negra", "", "1503-05=", "Galápagos", 0.83, -91.17, 1124, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Azul", "", "1503-06=", "Galápagos", 0.92, -91.408, 1640, 1, 0));
    m_volcanoes.push_back(CVolcano("Pinta", "", "1503-07=", "Galápagos", -0.58, -90.75, 780, 1, 0));
    m_volcanoes.push_back(CVolcano("Genovesa", "", "1503-081", "Galápagos", -0.32, -89.958, 64, 1, 0));
    m_volcanoes.push_back(CVolcano("Marchena", "", "1503-08=", "Galápagos", -0.33, -90.47, 343, 1, 0));
    m_volcanoes.push_back(CVolcano("Santa Cruz", "", "1503-091", "Galápagos", 0.62, -90.33, 864, 1, 0));
    m_volcanoes.push_back(CVolcano("Santiago", "", "1503-09=", "Galápagos", 0.22, -90.77, 920, 1, 0));
    m_volcanoes.push_back(CVolcano("San Cristóbal", "", "1503-12-", "Galápagos", 0.88, -89.5, 759, 1, 0));
    m_volcanoes.push_back(CVolcano("Quimsachata", "", "1504-00-", "Perú", 14.2, -71.33, 3923, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Auquihuato", "", "1504-001", "Perú", 15.07, -73.18, 4980, 1, 0));
    m_volcanoes.push_back(CVolcano("Sara Sara", "", "1504-002", "Perú", 15.33, -73.45, 5522, 1, 0));
    m_volcanoes.push_back(CVolcano("Coropuna", "", "1504-003", "Perú", 15.52, -72.65, 6377, 1, 0));
    m_volcanoes.push_back(CVolcano("Andahua-Orcopampa", "", "1504-004", "Perú", 15.42, -72.33, 4713, 1, 0));
    m_volcanoes.push_back(CVolcano("Huambo", "", "1504-005", "Perú", 15.83, -72.13, 4550, 1, 0));
    m_volcanoes.push_back(CVolcano("Sabancaya", "", "1504-006", "Perú", 15.78, -71.85, 5967, 1, 0));
    m_volcanoes.push_back(CVolcano("Nevado Chachani", "", "1504-007", "Perú", 16.191, -71.53, 6057, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Nicholson", "", "1504-008", "Perú", 16.261, -71.73, 2520, 1, 0));
    m_volcanoes.push_back(CVolcano("El Misti", "", "1504-01=", "Perú", 16.294, -71.409, 5822, 1, 0));
    m_volcanoes.push_back(CVolcano("Ubinas", "", "1504-02=", "Perú", 16.355, -70.903, 5672, 1, 0));
    m_volcanoes.push_back(CVolcano("Ticsani", "", "1504-031", "Perú", 16.755, -70.595, 5408, 1, 0));
    m_volcanoes.push_back(CVolcano("Huaynaputina", "", "1504-03=", "Perú", 16.608, -70.85, 4850, 1, 0));
    m_volcanoes.push_back(CVolcano("Tutupaca", "", "1504-04=", "Perú", 17.025, -70.358, 5815, 1, 0));
    m_volcanoes.push_back(CVolcano("Yucamane", "", "1504-05-", "Perú", 17.18, -70.2, 5550, 1, 0));
    m_volcanoes.push_back(CVolcano("Nevados Casiri", "", "1504-06-", "Perú", 17.47, -69.813, 5650, 1, 0));
    m_volcanoes.push_back(CVolcano("Taapaca", "", "1505-011", "Chile-N", 18.1, -69.5, 5860, 1, 0));
    m_volcanoes.push_back(CVolcano("Parinacota", "", "1505-012", "Chile-N", 18.17, -69.15, 6348, 1, 0));
    m_volcanoes.push_back(CVolcano("Tacora", "", "1505-01=", "Chile-N", 17.72, -69.77, 5980, 1, 0));
    m_volcanoes.push_back(CVolcano("Tambo Quemado", "", "1505-021", "Bolivia", 18.62, -68.75, 4215, 1, 0));
    m_volcanoes.push_back(CVolcano("Guallatiri", "", "1505-02=", "Chile-N", 18.42, -69.092, 6071, 1, 0));
    m_volcanoes.push_back(CVolcano("Tata Sabaya", "", "1505-032", "Bolivia", 19.13, -68.53, 5430, 1, 0));
    m_volcanoes.push_back(CVolcano("Laguna Jayu Khota", "", "1505-035", "Bolivia", 19.45, -67.42, 3650, 1, 0));
    m_volcanoes.push_back(CVolcano("Nuevo Mundo", "", "1505-036", "Bolivia", 19.78, -66.48, 5438, 1, 0));
    m_volcanoes.push_back(CVolcano("Isluga", "", "1505-03=", "Chile-N", 19.15, -68.83, 5550, 1, 0));
    m_volcanoes.push_back(CVolcano("Pampa Luxsar", "", "1505-042", "Bolivia", 20.85, -68.2, 5543, 1, 0));
    m_volcanoes.push_back(CVolcano("Irruputuncu", "", "1505-04=", "Chile-N", 20.73, -68.55, 5163, 1, 0));
    m_volcanoes.push_back(CVolcano("Olca-Paruma", "", "1505-05=", "Chile-N", 20.93, -68.48, 5407, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro del Azufre", "", "1505-061", "Chile-N", 21.787, -68.237, 5846, 1, 0));
    m_volcanoes.push_back(CVolcano("Ollagüe", "", "1505-06=", "Chile-N", 21.3, -68.18, 5868, 1, 0));
    m_volcanoes.push_back(CVolcano("San Pedro", "", "1505-07=", "Chile-N", 21.88, -68.4, 6145, 1, 0));
    m_volcanoes.push_back(CVolcano("Sairecabur", "", "1505-091", "Chile-N", 22.72, -67.892, 5971, 1, 0));
    m_volcanoes.push_back(CVolcano("Licancabur", "", "1505-092", "Chile-N", 22.83, -67.88, 5916, 1, 0));
    m_volcanoes.push_back(CVolcano("Guayaques", "", "1505-093", "Chile-N", 22.895, -67.566, 5598, 1, 0));
    m_volcanoes.push_back(CVolcano("Purico Complex", "", "1505-094", "Chile-N", 23, -67.75, 5703, 1, 0));
    m_volcanoes.push_back(CVolcano("Colachi", "", "1505-095", "Chile-N", 23.236, -67.645, 5631, 1, 0));
    m_volcanoes.push_back(CVolcano("Acamarachi", "", "1505-096", "Chile-N", 23.3, -67.62, 6046, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Overo", "", "1505-097", "Chile-N", 23.52, -67.67, 4555, 1, 0));
    m_volcanoes.push_back(CVolcano("Chiliques", "", "1505-098", "Chile-N", 23.58, -67.7, 5778, 1, 0));
    m_volcanoes.push_back(CVolcano("Putana", "", "1505-09=", "Chile-N", 22.55, -67.85, 5890, 1, 0));
    m_volcanoes.push_back(CVolcano("Cordón de Puntas Negras", "", "1505-101", "Chile-N", 23.743, -67.534, 5852, 1, 0));
    m_volcanoes.push_back(CVolcano("Miñiques", "", "1505-102", "Chile-N", 23.82, -67.77, 5910, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Tujle", "", "1505-103", "Chile-N", 23.83, -67.95, 3550, 1, 0));
    m_volcanoes.push_back(CVolcano("Caichinque", "", "1505-104", "Chile-N", 23.95, -67.73, 4450, 1, 0));
    m_volcanoes.push_back(CVolcano("Tilocalar", "", "1505-105", "Chile-N", 23.97, -68.13, 3116, 1, 0));
    m_volcanoes.push_back(CVolcano("El Negrillar", "", "1505-106", "Chile-N", 24.18, -68.25, 3500, 1, 0));
    m_volcanoes.push_back(CVolcano("Pular", "", "1505-107", "Chile-N", 24.188, -68.054, 6233, 1, 0));
    m_volcanoes.push_back(CVolcano("La Negrillar", "", "1505-108", "Chile-N", 24.28, -68.6, 4109, 1, 0));
    m_volcanoes.push_back(CVolcano("Socompa", "", "1505-109", "Chile-N", 24.4, -68.25, 6051, 1, 0));
    m_volcanoes.push_back(CVolcano("Láscar", "", "1505-10=", "Chile-N", 23.37, -67.73, 5592, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Escorial", "", "1505-112", "Chile-N", 25.08, -68.37, 5447, 1, 0));
    m_volcanoes.push_back(CVolcano("Llullaillaco", "", "1505-11=", "Chile-N", 24.72, -68.53, 6739, 1, 0));
    m_volcanoes.push_back(CVolcano("Cordón del Azufre", "", "1505-121", "Chile-N", 25.33, -68.52, 5463, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Bayo", "", "1505-122", "Chile-N", 25.42, -68.58, 5401, 1, 0));
    m_volcanoes.push_back(CVolcano("Sierra Nevada", "", "1505-123", "Chile-N", 26.48, -68.58, 6127, 1, 0));
    m_volcanoes.push_back(CVolcano("Falso Azufre", "", "1505-124", "Chile-N", 26.8, -68.37, 5890, 1, 0));
    m_volcanoes.push_back(CVolcano("Nevado de Incahuasi", "", "1505-125", "Chile-N", 27.042, -68.28, 6621, 1, 0));
    m_volcanoes.push_back(CVolcano("Lastarria", "", "1505-12=", "Chile-N", 25.17, -68.5, 5697, 1, 0));
    m_volcanoes.push_back(CVolcano("El Solo", "", "1505-131", "Chile-N", 27.108, -68.72, 6190, 1, 0));
    m_volcanoes.push_back(CVolcano("NevadosOjos del Salado", "", "1505-13=", "Chile-N", 27.12, -68.55, 6887, 1, 0));
    m_volcanoes.push_back(CVolcano("Copiapó", "", "1505-14-", "Chile-N", 27.3, -69.13, 6052, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Tuzgle", "", "1505-15-", "Argentina", 24.05, -66.48, 5500, 1, 0));
    m_volcanoes.push_back(CVolcano("Aracar", "", "1505-16-", "Argentina", 24.25, -67.77, 6082, 1, 0));
    m_volcanoes.push_back(CVolcano("Antofagasta de la Sierra", "", "1505-18-", "Argentina", 26.08, -67.5, 4000, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro el Cóndor", "", "1505-19-", "Argentina", 26.62, -68.35, 6532, 1, 0));
    m_volcanoes.push_back(CVolcano("Peinado", "", "1505-20-", "Argentina", 26.62, -68.15, 5740, 1, 0));
    m_volcanoes.push_back(CVolcano("Robledo", "", "1505-21-", "Argentina", 26.77, -67.72, 4400, 1, 0));
    m_volcanoes.push_back(CVolcano("Tipas", "", "1505-22-", "Argentina", 27.2, -68.55, 6660, 1, 0));
    m_volcanoes.push_back(CVolcano("Easter Island", "", "1506-011", "Chile-Is", 27.15, -109.38, 511, 1, 0));
    m_volcanoes.push_back(CVolcano("San Félix", "", "1506-01=", "Chile-Is", 26.28, -80.12, 193, 1, 0));
    m_volcanoes.push_back(CVolcano("Robinson Crusoe", "", "1506-02=", "Chile-Is", 33.658, -78.85, 922, 1, 0));
    m_volcanoes.push_back(CVolcano("Tupungatito", "", "1507-01=", "Chile-C", 33.4, -69.8, 6000, 1, 0));
    m_volcanoes.push_back(CVolcano("Maipo", "", "1507-021", "Chile-C", 34.161, -69.833, 5264, 1, 0));
    m_volcanoes.push_back(CVolcano("Palomo", "", "1507-022", "Chile-C", 34.608, -70.295, 4860, 1, 0));
    m_volcanoes.push_back(CVolcano("Caldera del Atuel", "", "1507-023", "Argentina", 34.65, -70.05, 5189, 1, 0));
    m_volcanoes.push_back(CVolcano("Risco Plateado", "", "1507-024", "Argentina", 34.93, -70, 4999, 1, 0));
    m_volcanoes.push_back(CVolcano("San José", "", "1507-02=", "Chile-C", 33.782, -69.897, 5856, 1, 0));
    m_volcanoes.push_back(CVolcano("Tinguiririca", "", "1507-03=", "Chile-C", 34.814, -70.352, 4280, 1, 0));
    m_volcanoes.push_back(CVolcano("Calabozos", "", "1507-042", "Chile-C", 35.558, -70.496, 3508, 1, 0));
    m_volcanoes.push_back(CVolcano("Planchón-Peteroa", "", "1507-04=", "Chile-C", 35.24, -70.57, 4107, 1, 0));
    m_volcanoes.push_back(CVolcano("Descabezado Grande", "", "1507-05=", "Chile-C", 35.58, -70.75, 3953, 1, 0));
    m_volcanoes.push_back(CVolcano("Laguna del Maule", "", "1507-061", "Chile-C", 36.02, -70.58, 3092, 1, 0));
    m_volcanoes.push_back(CVolcano("San Pedro-Pellado", "", "1507-062", "Chile-C", 35.989, -70.849, 3621, 1, 0));
    m_volcanoes.push_back(CVolcano("Nevado de Longaví", "", "1507-063", "Chile-C", 36.193, -71.161, 3242, 1, 0));
    m_volcanoes.push_back(CVolcano("Lomas Blancas", "", "1507-064", "Chile-C", 36.286, -71.009, 2268, 1, 0));
    m_volcanoes.push_back(CVolcano("Resago", "", "1507-065", "Chile-C", 36.45, -70.92, 1890, 1, 0));
    m_volcanoes.push_back(CVolcano("Payún Matru", "", "1507-066", "Argentina", 36.42, -69.2, 3680, 1, 0));
    m_volcanoes.push_back(CVolcano("Domuyo", "", "1507-067", "Argentina", 36.58, -70.42, 4709, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Azul", "", "1507-06=", "Chile-C", 35.653, -70.761, 3788, 1, 0));
    m_volcanoes.push_back(CVolcano("Cochiquito Volc Group", "", "1507-071", "Argentina", 36.77, -69.82, 1435, 1, 0));
    m_volcanoes.push_back(CVolcano("Tromen", "", "1507-072", "Argentina", 37.142, -70.03, 3978, 1, 0));
    m_volcanoes.push_back(CVolcano("Puesto Cortaderas", "", "1507-073", "Argentina", 37.57, -69.62, 970, 1, 0));
    m_volcanoes.push_back(CVolcano("Nevados de Chillán", "", "1507-07=", "Chile-C", 36.863, -71.377, 3212, 1, 0));
    m_volcanoes.push_back(CVolcano("Trocon", "", "1507-081", "Argentina", 37.73, -70.9, 2500, 1, 0));
    m_volcanoes.push_back(CVolcano("Antuco", "", "1507-08=", "Chile-C", 37.406, -71.349, 2979, 1, 0));
    m_volcanoes.push_back(CVolcano("Callaqui", "", "1507-091", "Chile-C", 37.92, -71.45, 3164, 1, 0));
    m_volcanoes.push_back(CVolcano("Laguna Mariñaqui", "", "1507-092", "Chile-C", 38.27, -71.1, 2143, 1, 0));
    m_volcanoes.push_back(CVolcano("Tolguaca", "", "1507-093", "Chile-C", 38.31, -71.645, 2806, 1, 0));
    m_volcanoes.push_back(CVolcano("Copahue", "", "1507-09=", "Chile-C", 37.85, -71.17, 2997, 1, 0));
    m_volcanoes.push_back(CVolcano("Lonquimay", "", "1507-10=", "Chile-C", 38.377, -71.58, 2865, 1, 0));
    m_volcanoes.push_back(CVolcano("Sollipulli", "", "1507-111", "Chile-C", 38.97, -71.52, 2282, 1, 0));
    m_volcanoes.push_back(CVolcano("Caburgua-Huelemolle", "", "1507-112", "Chile-C", 39.25, -71.7, 1496, 1, 0));
    m_volcanoes.push_back(CVolcano("Quetrupillan", "", "1507-121", "Chile-C", 39.5, -71.7, 2360, 1, 0));
    m_volcanoes.push_back(CVolcano("Lanín", "", "1507-122", "Chile-C", 39.633, -71.5, 3747, 1, 0));
    m_volcanoes.push_back(CVolcano("Huanquihue Group", "", "1507-123", "Argentina", 39.88, -71.58, 2139, 1, 0));
    m_volcanoes.push_back(CVolcano("Mocho-Choshuenco", "", "1507-13=", "Chile-C", 39.927, -72.027, 2422, 1, 0));
    m_volcanoes.push_back(CVolcano("Carrán-Los Venados", "", "1507-14=", "Chile-C", 40.35, -72.07, 1114, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Pantoja", "", "1507-152", "Chile-C", 40.77, -71.95, 2024, 1, 0));
    m_volcanoes.push_back(CVolcano("Antillanca Group", "", "1507-153", "Chile-C", 40.771, -72.153, 1990, 1, 0));
    m_volcanoes.push_back(CVolcano("Puyehue-Cordón Caulle", "", "1507-15=", "Chile-C", 40.59, -72.117, 2236, 1, 0));
    m_volcanoes.push_back(CVolcano("Puntiagudo-Cordón Cenizos", "", "1507-16-", "Chile-C", 40.969, -72.264, 2493, 1, 0));
    m_volcanoes.push_back(CVolcano("Tronador", "", "1508-011", "Chile-S", 41.157, -71.885, 3491, 1, 0));
    m_volcanoes.push_back(CVolcano("Cayutué-La Viguería", "", "1508-012", "Chile-S", 41.25, -72.27, 506, 1, 0));
    m_volcanoes.push_back(CVolcano("Osorno", "", "1508-01=", "Chile-S", 41.1, -72.493, 2652, 1, 0));
    m_volcanoes.push_back(CVolcano("Cuernos del Diablo", "", "1508-021", "Chile-S", 41.4, -72, 1862, 1, 0));
    m_volcanoes.push_back(CVolcano("Yate", "", "1508-022", "Chile-S", 41.755, -72.396, 2187, 1, 0));
    m_volcanoes.push_back(CVolcano("Hornopirén", "", "1508-023", "Chile-S", 41.874, -72.431, 1572, 1, 0));
    m_volcanoes.push_back(CVolcano("Apagado", "", "1508-024", "Chile-S", 41.88, -72.58, 1210, 1, 0));
    m_volcanoes.push_back(CVolcano("Crater Basalt Volc Field", "", "1508-025", "Chile-S/Argentina", 42.02, -70.18, 1359, 1, 0));
    m_volcanoes.push_back(CVolcano("Calbuco", "", "1508-02=", "Chile-S", 41.326, -72.614, 2003, 1, 0));
    m_volcanoes.push_back(CVolcano("Huequi", "", "1508-03=", "Chile-S", 42.377, -72.578, 1318, 1, 0));
    m_volcanoes.push_back(CVolcano("Chaitén", "", "1508-041", "Chile-S", 42.833, -72.646, 1122, 1, 0));
    m_volcanoes.push_back(CVolcano("Minchinmávida", "", "1508-04=", "Chile-S", 42.793, -72.439, 2404, 1, 0));
    m_volcanoes.push_back(CVolcano("Yanteles", "", "1508-050", "Chile-S", 43.5, -72.8, 2042, 1, 0));
    m_volcanoes.push_back(CVolcano("Palena Volc Group", "", "1508-051", "Chile-S", 43.78, -72.47, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Melimoyu", "", "1508-052", "Chile-S", 44.08, -72.88, 2400, 1, 0));
    m_volcanoes.push_back(CVolcano("Puyuhuapi", "", "1508-053", "Chile-S", 44.3, -72.53, 524, 1, 0));
    m_volcanoes.push_back(CVolcano("Mentolat", "", "1508-054", "Chile-S", 44.7, -73.08, 1660, 1, 0));
    m_volcanoes.push_back(CVolcano("Cay", "", "1508-055", "Chile-S", 45.059, -72.984, 2090, 1, 0));
    m_volcanoes.push_back(CVolcano("Maca", "", "1508-056", "Chile-S", 45.1, -73.17, 2960, 1, 0));
    m_volcanoes.push_back(CVolcano("Cerro Hudson", "", "1508-057", "Chile-S", 45.9, -72.97, 1905, 1, 0));
    m_volcanoes.push_back(CVolcano("Río Murta", "", "1508-058", "Chile-S", 46.17, -72.67, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Arenales", "", "1508-059", "Chile-S", 47.2, -73.48, 3437, 1, 0));
    m_volcanoes.push_back(CVolcano("Corcovado", "", "1508-05=", "Chile-S", 43.18, -72.8, 2300, 1, 0));
    m_volcanoes.push_back(CVolcano("Viedma", "", "1508-061", "Argentina", 49.358, -73.28, 1500, 1, 0));
    m_volcanoes.push_back(CVolcano("Aguilera", "", "1508-062", "Chile-S", 50.33, -73.75, 2546, 1, 0));
    m_volcanoes.push_back(CVolcano("Reclus", "", "1508-063", "Chile-S", 50.964, -73.585, 1000, 1, 0));
    m_volcanoes.push_back(CVolcano("Lautaro", "", "1508-06=", "Chile-S", 49.02, -73.55, 3607, 1, 0));
    m_volcanoes.push_back(CVolcano("Monte Burney", "", "1508-07=", "Chile-S", 52.33, -73.4, 1758, 1, 0));
    m_volcanoes.push_back(CVolcano("Palei-Aike Volc Field", "", "1508-08-", "Chile-S", 52, -70, 282, 1, 0));
    m_volcanoes.push_back(CVolcano("Fueguino", "", "1508-09-", "Chile-S", 54.95, -70.25, 150, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_16()
{
    m_volcanoes.push_back(CVolcano("Saba", "", "1600-01=", "W Indies", -17.63, -63.23, 887, 1, 0));
    m_volcanoes.push_back(CVolcano("The Quill", "", "1600-02=", "W Indies", -17.478, -62.96, 601, 1, 0));
    m_volcanoes.push_back(CVolcano("Liamuiga", "", "1600-03=", "W Indies", -17.37, -62.8, 1156, 1, 0));
    m_volcanoes.push_back(CVolcano("Nevis Peak", "", "1600-04=", "W Indies", -17.15, -62.58, 985, 1, 0));
    m_volcanoes.push_back(CVolcano("Soufrière Hills", "", "1600-05=", "W Indies", -16.72, -62.18, 915, 1, 0));
    m_volcanoes.push_back(CVolcano("Morne aux Diables", "", "1600-08=", "W Indies", -15.612, -61.43, 861, 1, 0));
    m_volcanoes.push_back(CVolcano("Morne Diablotins", "", "1600-09=", "W Indies", -15.503, -61.397, 1430, 1, 0));
    m_volcanoes.push_back(CVolcano("Morne Watt", "", "1600-101", "W Indies", -15.307, -61.305, 1224, 1, 0));
    m_volcanoes.push_back(CVolcano("Morne Trois Pitons", "", "1600-10=", "W Indies", -15.37, -61.33, 1387, 1, 0));
    m_volcanoes.push_back(CVolcano("MornePlat Pays", "", "1600-11=", "W Indies", -15.255, -61.341, 940, 1, 0));
    m_volcanoes.push_back(CVolcano("Pelée", "", "1600-12=", "W Indies", -14.82, -61.17, 1397, 1, 0));
    m_volcanoes.push_back(CVolcano("Qualibou", "", "1600-14=", "W Indies", -13.83, -61.05, 777, 1, 0));
    m_volcanoes.push_back(CVolcano("Soufrière St. Vincent", "", "1600-15=", "W Indies", -13.33, -61.18, 1220, 1, 0));
    m_volcanoes.push_back(CVolcano("Kick 'em Jenny", "", "1600-16=", "W Indies", -12.3, -61.64, -185, 1, 0));
    m_volcanoes.push_back(CVolcano("St. Catherine", "", "1600-17=", "W Indies", -12.15, -61.67, 840, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_17()
{
    m_volcanoes.push_back(CVolcano("Snaefellsjökull", "", "1700-01=", "Iceland-W", -64.8, -23.78, 1448, 1, 0));
    m_volcanoes.push_back(CVolcano("Helgrindur", "", "1700-02=", "Iceland-W", -64.87, -23.25, 647, 1, 0));
    m_volcanoes.push_back(CVolcano("Ljósufjöll", "", "1700-03=", "Iceland-W", -64.87, -22.23, 1063, 1, 0));
    m_volcanoes.push_back(CVolcano("Reykjanes", "", "1701-02=", "Iceland-SW", -63.88, -22.5, 230, 1, 0));
    m_volcanoes.push_back(CVolcano("Krísuvík", "", "1701-03=", "Iceland-SW", -63.93, -22.1, 379, 1, 0));
    m_volcanoes.push_back(CVolcano("Brennisteinsfjöll", "", "1701-04=", "Iceland-SW", -63.92, -21.83, 621, 1, 0));
    m_volcanoes.push_back(CVolcano("Hrómundartindur", "", "1701-051", "Iceland-S", -64.073, -21.202, 540, 1, 0));
    m_volcanoes.push_back(CVolcano("Hengill", "", "1701-05=", "Iceland-SW", -64.08, -21.32, 803, 1, 0));
    m_volcanoes.push_back(CVolcano("Grímsnes", "", "1701-06=", "Iceland-SW", -64.03, -20.87, 214, 1, 0));
    m_volcanoes.push_back(CVolcano("Prestahnukur", "", "1701-07=", "Iceland-SW", -64.6, -20.58, 1400, 1, 0));
    m_volcanoes.push_back(CVolcano("Hveravellir", "", "1701-08=", "Iceland-SW", -64.75, -19.98, 1360, 1, 0));
    m_volcanoes.push_back(CVolcano("Hofsjökull", "", "1701-09=", "Iceland-SW", -64.78, -18.92, 1782, 1, 0));
    m_volcanoes.push_back(CVolcano("Vestmannaeyjar", "", "1702-01=", "Iceland-S", -63.43, -20.28, 279, 1, 0));
    m_volcanoes.push_back(CVolcano("Eyjafjöll", "", "1702-02=", "Iceland-S", -63.63, -19.62, 1666, 1, 0));
    m_volcanoes.push_back(CVolcano("Katla", "", "1702-03=", "Iceland-S", -63.63, -19.05, 1512, 1, 0));
    m_volcanoes.push_back(CVolcano("Tindfjallajökull", "", "1702-04=", "Iceland-S", -63.78, -19.57, 1463, 1, 0));
    m_volcanoes.push_back(CVolcano("Torfajökull", "", "1702-05=", "Iceland-S", -63.92, -19.17, 1259, 1, 0));
    m_volcanoes.push_back(CVolcano("Hekla", "", "1702-07=", "Iceland-S", -63.98, -19.7, 1491, 1, 0));
    m_volcanoes.push_back(CVolcano("Grímsvötn", "", "1703-01=", "Iceland-NE", -64.42, -17.33, 1725, 1, 0));
    m_volcanoes.push_back(CVolcano("Bárdarbunga", "", "1703-03=", "Iceland-NE", -64.63, -17.53, 2009, 1, 0));
    m_volcanoes.push_back(CVolcano("Tungnafellsjökull", "", "1703-04=", "Iceland-NE", -64.73, -17.92, 1535, 1, 0));
    m_volcanoes.push_back(CVolcano("Kverkfjöll", "", "1703-05=", "Iceland-NE", -64.65, -16.72, 1929, 1, 0));
    m_volcanoes.push_back(CVolcano("Askja", "", "1703-06=", "Iceland-NE", -65.03, -16.75, 1516, 1, 0));
    m_volcanoes.push_back(CVolcano("Fremrinamur", "", "1703-07=", "Iceland-NE", -65.43, -16.65, 939, 1, 0));
    m_volcanoes.push_back(CVolcano("Krafla", "", "1703-08=", "Iceland-NE", -65.73, -16.78, 818, 1, 0));
    m_volcanoes.push_back(CVolcano("Theistareykjarbunga", "", "1703-09=", "Iceland-NE", -65.88, -16.83, 564, 1, 0));
    m_volcanoes.push_back(CVolcano("Tjörnes Fracture Zone", "", "1703-10=", "Iceland-N of", -66.3, -17.1, 0, 1, 0));
    m_volcanoes.push_back(CVolcano("Öraefajökull", "", "1704-01=", "Iceland-SE", -64, -16.65, 2119, 1, 0));
    m_volcanoes.push_back(CVolcano("Esjufjöll", "", "1704-02=", "Iceland-SE", -64.27, -16.65, 1760, 1, 0));
    m_volcanoes.push_back(CVolcano("Kolbeinsey Ridge", "", "1705-01=", "Iceland-N of", -66.67, -18.5, 5, 1, 0));
    m_volcanoes.push_back(CVolcano("Jan Mayen", "", "1706-01=", "Atlantic-N-Jan Mayen", -71.08, -8.17, 2277, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_18()
{
    m_volcanoes.push_back(CVolcano("Flores", "", "1802-001", "Azores", -39.462, -31.216, 914, 1, 0));
    m_volcanoes.push_back(CVolcano("Corvo", "", "1802-002", "Azores", -39.699, -31.111, 718, 1, 0));
    m_volcanoes.push_back(CVolcano("Fayal", "", "1802-01=", "Azores", -38.6, -28.73, 1043, 1, 0));
    m_volcanoes.push_back(CVolcano("Pico", "", "1802-02=", "Azores", -38.47, -28.4, 2351, 1, 0));
    m_volcanoes.push_back(CVolcano("San Jorge", "", "1802-03=", "Azores", -38.65, -28.08, 1053, 1, 0));
    m_volcanoes.push_back(CVolcano("Graciosa", "", "1802-04=", "Azores", -39.02, -27.97, 402, 1, 0));
    m_volcanoes.push_back(CVolcano("Terceira", "", "1802-05=", "Azores", -38.73, -27.32, 1023, 1, 0));
    m_volcanoes.push_back(CVolcano("Don Joao de Castro Bank", "", "1802-07=", "Azores", -38.23, -26.63, -13, 1, 0));
    m_volcanoes.push_back(CVolcano("Picos Volc System", "", "1802-081", "Azores", -37.78, -25.67, 350, 1, 0));
    m_volcanoes.push_back(CVolcano("Sete Cidades", "", "1802-08=", "Azores", -37.87, -25.78, 856, 1, 0));
    m_volcanoes.push_back(CVolcano("Agua de Pau", "", "1802-09=", "Azores", -37.77, -25.47, 947, 1, 0));
    m_volcanoes.push_back(CVolcano("Furnas", "", "1802-10=", "Azores", -37.77, -25.32, 805, 1, 0));
    m_volcanoes.push_back(CVolcano("Monaco Bank", "", "1802-11=", "Azores", -37.6, -25.88, -197, 1, 0));
    m_volcanoes.push_back(CVolcano("Madeira", "", "1802-12-", "Azores", -32.73, -16.97, 1862, 1, 0));
    m_volcanoes.push_back(CVolcano("La Palma", "", "1803-01-", "Canary Is", -28.57, -17.83, 2426, 1, 0));
    m_volcanoes.push_back(CVolcano("Hierro", "", "1803-02-", "Canary Is", -27.73, -18.03, 1500, 1, 0));
    m_volcanoes.push_back(CVolcano("Tenerife", "", "1803-03-", "Canary Is", -28.271, -16.641, 3715, 1, 0));
    m_volcanoes.push_back(CVolcano("Gran Canaria", "", "1803-04-", "Canary Is", -28, -15.58, 1950, 1, 0));
    m_volcanoes.push_back(CVolcano("Fuerteventura", "", "1803-05-", "Canary Is", -28.358, -14.02, 529, 1, 0));
    m_volcanoes.push_back(CVolcano("Lanzarote", "", "1803-06-", "Canary Is", -29.03, -13.63, 670, 1, 0));
    m_volcanoes.push_back(CVolcano("Fogo", "", "1804-01=", "Cape Verde Is", -14.95, -24.35, 2829, 1, 0));
    m_volcanoes.push_back(CVolcano("Brava", "", "1804-02-", "Cape Verde Is", -14.85, -24.72, 900, 1, 0));
    m_volcanoes.push_back(CVolcano("Sao Vicente", "", "1804-03-", "Cape Verde Is", -16.85, -24.97, 725, 1, 0));
    m_volcanoes.push_back(CVolcano("Ascensión", "", "1805-05-", "Atlantic-C", 7.95, -14.37, 858, 1, 0));
    m_volcanoes.push_back(CVolcano("Trindade", "", "1805-051", "Atlantic-C", 20.514, -29.331, 600, 1, 0));
    m_volcanoes.push_back(CVolcano("Tristan da Cunha", "", "1806-01=", "Atlantic-S", 37.092, -12.28, 2060, 1, 0));
    m_volcanoes.push_back(CVolcano("Bouvet", "", "1806-02-", "Atlantic-S", 54.42, 3.35, 780, 1, 0));
    m_volcanoes.push_back(CVolcano("Thompson Island", "", "1806-03-", "Atlantic-S", 53.93, 5.5, 0, 1, 0));
}
void CVolcanoInfo::InitializeDatabase_19()
{
    m_volcanoes.push_back(CVolcano("Young Island", "", "1900-011", "Antarctica", 66.42, 162.47, 1340, 1, 0));
    m_volcanoes.push_back(CVolcano("Sturge Island", "", "1900-012", "Antarctica", 67.4, 164.83, 1167, 1, 0));
    m_volcanoes.push_back(CVolcano("The Pleiades", "", "1900-013", "Antarctica", 72.67, 165.5, 3040, 1, 0));
    m_volcanoes.push_back(CVolcano("Melbourne", "", "1900-015", "Antarctica", 74.35, 164.7, 2732, 1, 0));
    m_volcanoes.push_back(CVolcano("Buckle Island", "", "1900-01=", "Antarctica", 66.78, 163.25, 1239, 1, 0));
    m_volcanoes.push_back(CVolcano("Royal Society Range", "", "1900-021", "Antarctica", 78.25, 163.33, 3000, 1, 0));
    m_volcanoes.push_back(CVolcano("Berlin", "", "1900-022", "Antarctica", 76.05, -136, 3478, 1, 0));
    m_volcanoes.push_back(CVolcano("Andrus", "", "1900-023", "Antarctica", 75.8, -132.33, 2978, 1, 0));
    m_volcanoes.push_back(CVolcano("Waesche", "", "1900-024", "Antarctica", 77.17, -126.88, 3292, 1, 0));
    m_volcanoes.push_back(CVolcano("Siple", "", "1900-025", "Antarctica", 73.43, -126.67, 3110, 1, 0));
    m_volcanoes.push_back(CVolcano("Toney Mountain", "", "1900-026", "Antarctica", 75.8, -115.83, 3595, 1, 0));
    m_volcanoes.push_back(CVolcano("Takahe", "", "1900-027", "Antarctica", 76.28, -112.08, 3460, 1, 0));
    m_volcanoes.push_back(CVolcano("Hudson Mountains", "", "1900-028", "Antarctica", 74.33, -99.42, 749, 1, 0));
    m_volcanoes.push_back(CVolcano("Peter I Island", "", "1900-029", "Antarctica", 68.85, -90.58, 1640, 1, 0));
    m_volcanoes.push_back(CVolcano("Erebus", "", "1900-02=", "Antarctica", 77.53, 167.17, 3794, 1, 0));
    m_volcanoes.push_back(CVolcano("Penguin Island", "", "1900-031", "Antarctica", 62.1, -57.93, 180, 1, 0));
    m_volcanoes.push_back(CVolcano("Deception Island", "", "1900-03=", "Antarctica", 62.97, -60.65, 576, 1, 0));
    m_volcanoes.push_back(CVolcano("Paulet", "", "1900-041", "Antarctica", 63.58, -55.77, 353, 1, 0));
    m_volcanoes.push_back(CVolcano("Bridgeman Island", "", "1900-04=", "Antarctica", 62.05, -56.75, 240, 1, 0));
    m_volcanoes.push_back(CVolcano("Seal Nunataks Group", "", "1900-05=", "Antarctica", 65.03, -60.05, 368, 1, 0));
    m_volcanoes.push_back(CVolcano("Thule Islands", "", "1900-07=", "Antarctica", 59.45, -27.37, 1075, 1, 0));
    m_volcanoes.push_back(CVolcano("Montagu Island", "", "1900-081", "Antarctica", 58.42, -26.33, 1370, 1, 0));
    m_volcanoes.push_back(CVolcano("Bristol Island", "", "1900-08=", "Antarctica", 59.03, -26.58, 1100, 1, 0));
    m_volcanoes.push_back(CVolcano("Michael", "", "1900-09=", "Antarctica", 57.78, -26.45, 990, 1, 0));
    m_volcanoes.push_back(CVolcano("Candlemas Island", "", "1900-10=", "Antarctica", 57.08, -26.67, 550, 1, 0));
    m_volcanoes.push_back(CVolcano("Hodson", "", "1900-11=", "Antarctica", 56.7, -27.15, 1005, 1, 0));
    m_volcanoes.push_back(CVolcano("Leskov Island", "", "1900-12=", "Antarctica", 56.67, -28.13, 190, 1, 0));
    m_volcanoes.push_back(CVolcano("Zavodovski", "", "1900-13=", "Antarctica", 56.3, -27.57, 551, 1, 0));
    m_volcanoes.push_back(CVolcano("Protector Shoal", "", "1900-14-", "Antarctica", 55.92, -28.08, -27, 1, 0));
}
}  // namespace novac
