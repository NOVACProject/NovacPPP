#pragma once

#include "CrossSectionData.h"
#include <PPPLib/CString.h>

namespace Evaluation
{
	// the options for the shift and squeeze
	enum SHIFT_TYPE{
		SHIFT_FREE,
		SHIFT_FIX,
		SHIFT_LINK,
		SHIFT_LIMIT
	};

	class CReferenceFile
	{
	public:
		CReferenceFile(void);
		~CReferenceFile(void);

		/** The name of the specie */
		novac::CString m_specieName;

		/** The path to the reference file */
		novac::CString m_path;

		/** assignment operator */
		CReferenceFile &operator=(const CReferenceFile &ref2);

		/** The option for the column value (normally SHIFT_FREE) */
		SHIFT_TYPE m_columnOption;

		/** The value for the column value (only used if m_columnOption is not SHIFT_FREE) */
		double m_columnValue;

		/** if m_columnOption is SHIFT_LIMIT,
			this is the maximum column value allowed
			and m_columnValue is the minimum column value allowed */
		double m_columnMaxValue;

		/** The option for the shift */
		SHIFT_TYPE m_shiftOption;

		/** The value for the shift */
		double m_shiftValue;

		/** if m_shiftOption is SHIFT_LIMIT,
			this is the maximum shift value allowed
			and m_shiftValue is the minimum shift value allowed */
		double m_shiftMaxValue;

		/** The option for the squeeze */
		SHIFT_TYPE m_squeezeOption;

		/** The value for the squeeze */
		double m_squeezeValue;

		/** if m_squeezeOption is SHIFT_LIMIT,
			this is the maximum squeeze value allowed
			and m_squeezeValue is the minimum squeeze value allowed */
		double m_squeezeMaxValue;

		/** this is true if this cross section file is filtered 
				already on disk */
		bool m_isFiltered;
	
		/** The actual data of this reference - file. This is equal to
			NULL if the reference data has not been read in yet, otherwise
			this will contain the data from the reference file. */
		CCrossSectionData *m_data;

		// ------------------------ METHODS ---------------------------

		/** Setting the column.
			if(SHIFT_TYPE) is SHIFT_LIMIT then 'value' is the lower limit
			and value2 is the upper limit 
			otherwise value2 is not used */
		void SetColumn(SHIFT_TYPE option, double value, double value2 = 1e16);

		/** Setting the shift
			if(SHIFT_TYPE) is SHIFT_LIMIT then 'value' is the lower limit
			and value2 is the upper limit 
			otherwise value2 is not used */
		void SetShift(SHIFT_TYPE option, double value, double value2 = 1e16);

		/** Setting the squeeze
			if(SHIFT_TYPE) is SHIFT_LIMIT then 'value' is the lower limit
			and value2 is the upper limit 
			otherwise value2 is not used */
		void SetSqueeze(SHIFT_TYPE option, double value, double value2 = 1e16);
		
		/** Reads the data of this reference file from disk. The file is taken
			from the member variable 'm_path' (which of course must be initialized first)
			and the result is written to 'm_data'. If this fails then m_data will be NULL */
		int ReadCrossSectionDataFromFile();
	};
}