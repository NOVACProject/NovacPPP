#ifndef NOVAC_PPP_MOLECULE_H
#define NOVAC_PPP_MOLECULE_H

#include <PPPLib/MFC/CString.h>

/** The class <b>CMolecule</b> is used to store information on a single
    trace gas species, such as it's molecular weight.
    This is also able to use these numbers to convert between columns in number
    of molecules and in mass.
*/
enum STANDARD_MOLECULE
{
    MOLEC_SO2,
    MOLEC_O3,
    MOLEC_BRO,
    MOLEC_NO2,
    MOLEC_HCHO
};

const double AVOGADROS_NUMBER = 6.02214179e23;

/** The class <b>CMolecule</b> is used to pass the properties of an individual
    molecule to different parts of the program. The main usage (for now) is
    passing on the molecular weight to convert from molec/cm2 to kg/m2.
*/

class CMolecule
{
public:

    CMolecule();

    CMolecule(STANDARD_MOLECULE molec);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The trivial name of this gas molecule */
    novac::CString m_name;

    /** The molecular weight of the molecule
        In grams per mol (g/mol) or (u/molecule) */
    double m_molecularWeight;

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Takes a column in molec/cm2 and converts it to kg/m2 */
    double Convert_MolecCm2_to_kgM2(double molec_cm2) const;

};

#endif  // NOVAC_PPP_MOLECULE_H