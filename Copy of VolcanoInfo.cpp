#include "StdAfx.h"
#include "volcanoinfo.h"
#include "Common/Common.h"

// The global database of volcanoes
CVolcanoInfo g_volcanoes;

CVolcanoInfo::CVolcanoInfo(void)
{
	// Setting up the information about the volcanoes
	int index = 0;

	m_name[index].Format("Arenal");
	m_simpleName[index].Format("arenal");
	m_number[index].Format("1405-033");
	m_country[index].Format("Costa Rica");
	m_peakLatitude[index]		= 10.46;
	m_peakLongitude[index]		= -84.7;
	m_peakHeight[index]			= 1670;
	m_hoursToGMT[index]			= -6; // Costa Rica
	m_observatory[index]		= 11;

	++index;
	m_name[index].Format("Po�s");
	m_simpleName[index].Format("poas");
	m_number[index].Format("1405-02=");
	m_country[index].Format("Costa Rica");
	m_peakLatitude[index]		= 10.20;
	m_peakLongitude[index]		= -84.23;
	m_peakHeight[index]			= 2708;
	m_hoursToGMT[index]			= -6;	// Costa Rica
	m_observatory[index]		= 12;

	++index;
	m_name[index].Format("Turrialba");
	m_simpleName[index].Format("turrialba");
	m_number[index].Format("1405-07=");
	m_country[index].Format("Costa Rica");
	m_peakLatitude[index]		= 10.025;
	m_peakLongitude[index]		= -83.767;
	m_peakHeight[index]			= 3340;
	m_hoursToGMT[index]			= -6;	// Costa Rica
	m_observatory[index]		= 12;

	++index;
	m_name[index].Format("Santa Ana");
	m_simpleName[index].Format("santa_ana");
	m_number[index].Format("1403-02=");
	m_country[index].Format("El Salvador");
	m_peakLatitude[index]		= 13.85;
	m_peakLongitude[index]		= -89.63;
	m_peakHeight[index]			= 2381;
	m_hoursToGMT[index]			= -6; // El Salvador
	m_observatory[index]		= 13;

	++index;
	m_name[index].Format("San Miguel");
	m_simpleName[index].Format("san_miguel");
	m_number[index].Format("1403-10=");
	m_country[index].Format("El Salvador");
	m_peakLatitude[index]		= 13.43;
	m_peakLongitude[index]		= -88.27;
	m_peakHeight[index]			= 2130;
	m_hoursToGMT[index]			= -6; // El Salvador
	m_observatory[index]		= 13;

	++index;
	m_name[index].Format("Popocat�petl");
	m_simpleName[index].Format("popocatepetl");
	m_number[index].Format("1401-09=");
	m_country[index].Format("M�xico");
	m_peakLatitude[index]		= 19.02;
	m_peakLongitude[index]		= -98.62;
	m_peakHeight[index]			= 5426;
	m_hoursToGMT[index]			= -6; // M�xico
	m_observatory[index]		= 17;

	++index;
	m_name[index].Format("Fuego de Colima");
	m_simpleName[index].Format("fuego_de_colima");
	m_number[index].Format("1401-04=");
	m_country[index].Format("M�xico");
	m_peakLatitude[index]		= 19.51;
	m_peakLongitude[index]		= -103.62;
	m_peakHeight[index]			= 3850;
	m_hoursToGMT[index]			= -6; // M�xico
	m_observatory[index]		= 17;

	++index;
	m_name[index].Format("San Crist�bal");
	m_simpleName[index].Format("san_cristobal");
	m_number[index].Format("1404-02=");
	m_country[index].Format("Nicaragua");
	m_peakLatitude[index]		= 12.70;
	m_peakLongitude[index]		= -87.0;
	m_peakHeight[index]			= 1745;
	m_hoursToGMT[index]			= -6;	// Nicaragua
	m_observatory[index]		= 4;

	++index;
	m_name[index].Format("Masaya");
	m_simpleName[index].Format("masaya");
	m_number[index].Format("1404-10=");
	m_country[index].Format("Nicaragua");
	m_peakLatitude[index]		= 11.98;
	m_peakLongitude[index]		= -86.16;
	m_peakHeight[index]			= 635;
	m_hoursToGMT[index]			= -6;	// Nicaragua
	m_observatory[index]		= 4;

	++index;
	m_name[index].Format("Galeras");
	m_simpleName[index].Format("galeras");
	m_number[index].Format("1501-08=");
	m_country[index].Format("Colombia");
	m_peakLatitude[index]		= 1.22;
	m_peakLongitude[index]		= -77.37;
	m_peakHeight[index]			= 4276;
	m_hoursToGMT[index]			= -5; // Colombia
	m_observatory[index]		= 5;

	++index;
	m_name[index].Format("Nevado del Ruiz");
	m_simpleName[index].Format("nevado_del_ruiz");
	m_number[index].Format("1501-02=");
	m_country[index].Format("Colombia");
	m_peakLatitude[index]		= 4.89;
	m_peakLongitude[index]		= -75.32;
	m_peakHeight[index]			= 5321;
	m_hoursToGMT[index]			= -5; // Colombia
	m_observatory[index]		= 5;

	++index;
	m_name[index].Format("Nevado del Huila");
	m_simpleName[index].Format("nevado_del_huila");
	m_number[index].Format("1501-05=");
	m_country[index].Format("Colombia");
	m_peakLatitude[index]		= 2.93;
	m_peakLongitude[index]		= -76.03;
	m_peakHeight[index]			= 5364 ;
	m_hoursToGMT[index]			= -5; // Colombia
	m_observatory[index]		= 5;

	++index;
	m_name[index].Format("Nyiragongo");
	m_simpleName[index].Format("nyiragongo");
	m_number[index].Format("0203-03=");
	m_country[index].Format("Democratic Republic of Congo");
	m_peakLatitude[index]		= -1.516732; //-1.40762;
	m_peakLongitude[index]		= 29.24668;  //29.2091;
	m_peakHeight[index]			= 3470;
	m_hoursToGMT[index]			= +2; // Dem. Rep. Congo
	m_observatory[index]		= 11;

	++index;
	m_name[index].Format("Nyamuragira");
	m_simpleName[index].Format("nyamuragira");
	m_number[index].Format("0203-02=");
	m_country[index].Format("Democratic Republic of Congo");
	m_peakLatitude[index]		= -1.41;
	m_peakLongitude[index]		= 29.20;
	m_peakHeight[index]			= 3058;
	m_hoursToGMT[index]			= +2; // Dem. Rep. Congo
	m_observatory[index]		= 11;

	++index;
	m_name[index].Format("Etna");
	m_simpleName[index].Format("etna");
	m_number[index].Format("0101-06=");
	m_country[index].Format("Italy");
	m_peakLatitude[index]		= 37.752;
	m_peakLongitude[index]		= 14.995;
	m_peakHeight[index]			= 3300;
	m_hoursToGMT[index]			= +1; // Italy
	m_observatory[index]		= 6;

	++index;
	m_name[index].Format("La Soufri�re");
	m_simpleName[index].Format("la_soufriere");
	m_number[index].Format("1600-06=");
	m_country[index].Format("France");
	m_peakLatitude[index]		= 16.05;
	m_peakLongitude[index]		= -61.67;
	m_peakHeight[index]			= 1467;
	m_hoursToGMT[index]			= 0; // France
	m_observatory[index]		= 9;

	++index;
	m_name[index].Format("Tungurahua");
	m_simpleName[index].Format("tungurahua");
	m_number[index].Format("1502-08=");
	m_country[index].Format("Ecuador");
	m_peakLatitude[index]		= -1.467;
	m_peakLongitude[index]		= -78.442;
	m_peakHeight[index]			= 5023;
	m_hoursToGMT[index]			= -5; // Ecuador
	m_observatory[index]		= 3;

	++index;
	m_name[index].Format("Cotopaxi");
	m_simpleName[index].Format("cotopaxi");
	m_number[index].Format("1502-05=");
	m_country[index].Format("Ecuador");
	m_peakLatitude[index]		= -0.677;
	m_peakLongitude[index]		= -78.436;
	m_peakHeight[index]			= 5911;
	m_hoursToGMT[index]			= -5;	// Ecuador
	m_observatory[index]		= 3;

	++index;
	m_name[index].Format("Pacaya");
	m_simpleName[index].Format("pacaya");
	m_number[index].Format("1402-11=");
	m_country[index].Format("Guatemala");
	m_peakLatitude[index]		= 14.381;
	m_peakLongitude[index]		= 90.601;
	m_peakHeight[index]			= 2552;
	m_hoursToGMT[index]			= -6; // Guatemala
	m_observatory[index]		= 8;

	++index;
	m_name[index].Format("Piton de la Fournaise");
	m_simpleName[index].Format("piton_de_la_fournaise");
	m_number[index].Format("0303-02=");
	m_country[index].Format("France");
	m_peakLatitude[index]		= -21.23;
	m_peakLongitude[index]		= 55.71;
	m_peakHeight[index]			= 2632;
	m_hoursToGMT[index]			= +4; // Reunion
	m_observatory[index]		= 9;

	++index;
	m_name[index].Format("Santiaguito");
	m_simpleName[index].Format("santiaguito");
	m_number[index].Format("1402-03=");
	m_country[index].Format("Guatemala");
	m_peakLatitude[index]		= 14.756;
	m_peakLongitude[index]		= 91.552;
	m_peakHeight[index]			= 3772;
	m_hoursToGMT[index]			= -6; // Guatemala
	m_observatory[index]		= 8;

	++index;
	m_name[index].Format("Fuego (Guatemala)");
	m_simpleName[index].Format("fuego_guatemala");
	m_number[index].Format("1402-09=");
	m_country[index].Format("Guatemala");
	m_peakLatitude[index]		= 14.473;
	m_peakLongitude[index]		= -90.880;
	m_peakHeight[index]			= 3763;
	m_hoursToGMT[index]			= -6; // Guatemala
	m_observatory[index]		= 8;

	++index;
	m_name[index].Format("Vulcano");
	m_simpleName[index].Format("vulcano");
	m_number[index].Format("0101-05=");
	m_country[index].Format("Italy");
	m_peakLatitude[index]		= 38.404;
	m_peakLongitude[index]		= 14.986;
	m_peakHeight[index]			= 500;
	m_hoursToGMT[index]			= +1; // Italy
	m_observatory[index]		= 7;

	++index;
	m_name[index].Format("Stromboli");
	m_simpleName[index].Format("stromboli");
	m_number[index].Format("0101-04=");
	m_country[index].Format("Italy");
	m_peakLatitude[index]		= 38.789;
	m_peakLongitude[index]		=  15.213;
	m_peakHeight[index]			= 924;
	m_hoursToGMT[index]			= +1; // Italy
	m_observatory[index]		= 7;
	

	++index;
	m_name[index].Format("Yasur");
	m_simpleName[index].Format("yasur");
	m_number[index].Format("1405-07=");
	m_country[index].Format("Indonesia");
	m_peakLatitude[index]		= -19.53;
	m_peakLongitude[index]		= 169.442;
	m_peakHeight[index]			= 361;
	m_hoursToGMT[index]			= +11; // Indonesia
	m_observatory[index]		= 2;


	++index;
	m_name[index].Format("Villarrica");
	m_simpleName[index].Format("villarrica");
	m_number[index].Format("1507-12=");
	m_country[index].Format("Chile");
	m_peakLatitude[index]		= -39.42;
	m_peakLongitude[index]		= -71.93;
	m_peakHeight[index]			= 2847;
	m_hoursToGMT[index]			= -7; // Chile
	m_observatory[index]		= 2;

	++index;
	m_name[index].Format("Llaima");
	m_simpleName[index].Format("llaima");
	m_number[index].Format("1507-11=");
	m_country[index].Format("Chile");
	m_peakLatitude[index]		= -38.692;
	m_peakLongitude[index]		= - 71.729;
	m_peakHeight[index]			= 3125;
	m_hoursToGMT[index]			= -7; // Chile
	m_observatory[index]		= 2;

	++index;
	m_name[index].Format("Antisana");
	m_simpleName[index].Format("antisana");
	m_number[index].Format("1502-03=");
	m_country[index].Format("Ecuador");
	m_peakLatitude[index]		= -0.481;
	m_peakLongitude[index]		= -78.141;
	m_peakHeight[index]			= 5753;
	m_hoursToGMT[index]			= -5; // Ecuador
	m_observatory[index]		= 3;

	++index;
	m_name[index].Format("El Reventador");
	m_simpleName[index].Format("el_reventador");
	m_number[index].Format("1502-01=");
	m_country[index].Format("Ecuador");
	m_peakLatitude[index]		= -0.077;
	m_peakLongitude[index]		= -77.656;
	m_peakHeight[index]			= 3562;
	m_hoursToGMT[index]			= -5; // Ecuador
	m_observatory[index]		= 3;

	++index;
	m_name[index].Format("Sangay");
	m_simpleName[index].Format("sangay");
	m_number[index].Format("1502-09=");
	m_country[index].Format("Ecuador");
	m_peakLatitude[index]		= -2.002;
	m_peakLongitude[index]		= -78.341;
	m_peakHeight[index]			= 5230;
	m_hoursToGMT[index]			= -5; // Ecuador
	m_observatory[index]		= 3;

	++index;
	m_name[index].Format("Cayambe");
	m_simpleName[index].Format("cayambe");
	m_number[index].Format("1502-006");
	m_country[index].Format("Ecuador");
	m_peakLatitude[index]		= 0.029;
	m_peakLongitude[index]		= -77.986;
	m_peakHeight[index]			= 5790;
	m_hoursToGMT[index]			= -5; // Ecuador
	m_observatory[index]		= 3;
	
	++index;
	m_name[index].Format("Guagua Pichincha");
	m_simpleName[index].Format("guagua_pichincha");
	m_number[index].Format("1502-02=");
	m_country[index].Format("Ecuador");
	m_peakLatitude[index]		= -0.171;
	m_peakLongitude[index]		= -78.598;
	m_peakHeight[index]			= 4784;
	m_hoursToGMT[index]			= -5; // Ecuador
	m_observatory[index]		= 3;

	++index;
	m_name[index].Format("Merapi");
	m_simpleName[index].Format("merapi");
	m_number[index].Format("0603-25=");
	m_country[index].Format("Indonesia");
	m_peakLatitude[index]		= -7.542;
	m_peakLongitude[index]		= 110.442;
	m_peakHeight[index]			= 2968;
	m_hoursToGMT[index]			= +11; // Indonesia
	m_observatory[index]		= 1;

	++index;
	m_name[index].Format("Kilauea");
	m_simpleName[index].Format("Kilauea");
	m_number[index].Format("1302-01-");
	m_country[index].Format("USA");
	m_peakLatitude[index]		= 19.421;
	m_peakLongitude[index]		= -155.287;
	m_peakHeight[index]			= 1222;
	m_hoursToGMT[index]			= -10; // USA, Hawaii
	m_observatory[index]		= 18;
	
	++index;
	m_name[index].Format("Chalmers");
	m_simpleName[index].Format("chalmers");
	m_number[index].Format("0000-000");
	m_country[index].Format("Sweden");
	m_peakLatitude[index]		= 0.0;
	m_peakLongitude[index]		= 0.0;
	m_peakHeight[index]			= 0;
	m_hoursToGMT[index]			= +1;
	m_observatory[index]		= 1;

	++index;
	m_name[index].Format("Harestua");
	m_simpleName[index].Format("harestua");
	m_number[index].Format("0000-000");
	m_country[index].Format("Norway");
	m_peakLatitude[index]		= 60.21;
	m_peakLongitude[index]		= 10.75;
	m_peakHeight[index]			= 600;
	m_hoursToGMT[index]			= +1;
	m_observatory[index]		= 1;

	m_volcanoNum				= index+1;

	m_preConfiguredVolcanoNum	  = m_volcanoNum;
}

CVolcanoInfo::~CVolcanoInfo(void)
{
}

int CVolcanoInfo::GetVolcanoIndex(const CString &name){
	static unsigned int lastIndex = 0;
	
	// first try with the same volcano as the last time
	//	this function was called...
	if(Equals(name, g_volcanoes.m_name[lastIndex]) ||
	   Equals(name, g_volcanoes.m_simpleName[lastIndex])){
		return (int)lastIndex;
	}	

	// this was not the same volcano as last time, search for it...
	for(unsigned int k = 0; k < g_volcanoes.m_volcanoNum; ++k){
		if(Equals(name, g_volcanoes.m_name[k]) || 
		   Equals(name, g_volcanoes.m_simpleName[k])){
			lastIndex = k;
			return (int)k;
		}
	}
	return -1; // no volcano found
}

void CVolcanoInfo::GetVolcanoName(unsigned int index, CString &name){
	if(index >= g_volcanoes.m_volcanoNum){
		name.Format("");
	}else{
		name.Format(g_volcanoes.m_name[index]);
	}
}

void CVolcanoInfo::GetVolcanoCode(unsigned int index, CString &code){
	if(index >= g_volcanoes.m_volcanoNum){
		code.Format("");
	}else{
		code.Format(g_volcanoes.m_number[index]);
	}
}

void CVolcanoInfo::GetSimpleVolcanoName(unsigned int index, CString &name){
	if(index >= g_volcanoes.m_volcanoNum){
		name.Format("");
	}else{
		name.Format(g_volcanoes.m_simpleName[index]);
	}
}

