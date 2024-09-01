#include <PPPLib/Molecule.h>

CMolecule::CMolecule()
: m_name("SO2"), m_molecularWeight(64.0638)
{
}

CMolecule::CMolecule(STANDARD_MOLECULE molec)
{
    switch (molec)
    {
    case MOLEC_SO2: m_name.Format("SO2");	m_molecularWeight = 64.0638; break;
    case MOLEC_O3: m_name.Format("O3");	m_molecularWeight = 47.9982; break;
    case MOLEC_BRO: m_name.Format("BrO");	m_molecularWeight = 95.8980; break;
    case MOLEC_NO2: m_name.Format("NO2");	m_molecularWeight = 46.0055; break;
    case MOLEC_HCHO: m_name.Format("HCHO"); m_molecularWeight = 30.0206; break;
    default: m_name.Format("SO2");	m_molecularWeight = 64.0638; break;
    }
}

double CMolecule::Convert_MolecCm2_to_kgM2(double molec_cm2) const
{
    return molec_cm2 * (10.0) * this->m_molecularWeight / AVOGADROS_NUMBER;
}
