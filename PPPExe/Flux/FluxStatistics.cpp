#include "FluxStatistics.h"

using namespace Flux;

CFluxStatistics::CMeasurementDay::CMeasurementDay(){

}
CFluxStatistics::CMeasurementDay::~CMeasurementDay(){

}
CFluxStatistics::CMeasurementDay &CFluxStatistics::CMeasurementDay::operator =(const CFluxStatistics::CMeasurementDay &m){
	this->day = m.day;
	
	auto p = m.fluxList.GetHeadPosition();
	while(p != nullptr){
		fluxList.AddTail(CFluxResult(m.fluxList.GetNext(p)));
	}
	
	return *this;
}

void CFluxStatistics::CMeasurementDay::GetHeaderLine(novac::CString &str, novac::CList <novac::CString, novac::CString&> &instruments){
	// the statistics of the fluxes
	str.Format("Date\tAverageFlux(kg/s)\tMedianFlux(kg/s)\tStdFlux(kg/s)\tnMeasurements");
	
	// the statistics of the instruments
	auto p = instruments.GetHeadPosition();
	while(p != nullptr){
		str.AppendFormat("\t#MeasFrom_%s", (const char*)instruments.GetNext(p));
	}
	str.AppendFormat("\n");
	
	return;
}

void CFluxStatistics::CMeasurementDay::GetStatistics(novac::CString &str, novac::CList <novac::CString, novac::CString &> &instruments){
	double average, median, std;
	long nMeasurements	= (long)this->fluxList.GetCount();
	double *data		= new double[nMeasurements];
	double *sortedData	= new double[nMeasurements];
	int nMeasurementsFromThisInstrument = 0;
	int k = 0;
	
	// Write the day
	str.Format("%04d.%02d.%02d\t", day.year, day.month, day.day);
	
	// copy the data into the array
	auto p = fluxList.GetHeadPosition();
	while(p != nullptr){
		data[k++] = fluxList.GetNext(p).m_flux;
	}
	
	// get the statistical data
	average = Average(data, nMeasurements);
	std		= Std(data, nMeasurements);
	
	// sort the data to get the median
	FindNLowest(data, nMeasurements, sortedData, nMeasurements);
	if(nMeasurements % 2 == 0)
		median = Average(sortedData + nMeasurements / 2-1, 2);
	else
		median = sortedData[nMeasurements / 2 + 1];
	
	// write what we now know to the string
	str.AppendFormat("%.2lf\t%.2lf\t%.2lf\t%d", average, median, std, nMeasurements);
	
	// go through the data and see how many points we have from each instrument
	auto p2 = instruments.GetHeadPosition();
	while(p2 != nullptr){
		novac::CString &serial = instruments.GetNext(p2);
		nMeasurementsFromThisInstrument = 0;
		auto p = fluxList.GetHeadPosition();
		while(p != nullptr){
			const CFluxResult &flux = fluxList.GetNext(p);
			if(Equals(flux.m_instrument, serial)){
				++nMeasurementsFromThisInstrument;
			}
		}
		
		str.AppendFormat("\t%d", nMeasurementsFromThisInstrument);
	}
	str.AppendFormat("\n");
	
	
	delete [] data;
	delete [] sortedData;
	return;
}


CFluxStatistics::CFluxStatistics(void)
{
	this->Clear();
}

CFluxStatistics::~CFluxStatistics(void)
{
	this->Clear();
}


/** Clears all information here */
void CFluxStatistics::Clear(){
	m_measurements.RemoveAll();
}

/** Attaches the supplied list of flux results to the current set
	of measured data. */
void CFluxStatistics::AttachFluxList(novac::CList <CFluxResult, CFluxResult&>& calculatedFluxes){
	auto p = calculatedFluxes.GetHeadPosition();
	while(p != nullptr){
		AttachFlux(calculatedFluxes.GetNext(p));
	}
}

/** Attaches the given flux result to the current set of 
	measured data */
void CFluxStatistics::AttachFlux(const CFluxResult &result){
	CFluxResult r = result; // make a local copy of the result
	CMeasurementDay measday;
	novac::CDateTime resultDay = novac::CDateTime(result.m_startTime.year, result.m_startTime.month, result.m_startTime.day, 0, 0, 0);

	// find out if we know about this instrument
	bool foundInstrument = false;
	auto instr_p = m_instruments.GetHeadPosition();
	while(instr_p != nullptr){
		if(Equals(m_instruments.GetNext(instr_p), result.m_instrument)){
			foundInstrument = true;
			break;
		}
	}
	if(!foundInstrument){
		m_instruments.AddHead(novac::CString(result.m_instrument));
	}

	// Loop through the list of measurementdays to find the 
	//	right place for the result and insert it.
	auto meas_p = m_measurements.GetHeadPosition();
	while(meas_p != nullptr){
		CMeasurementDay &d = m_measurements.GetAt(meas_p);
		
		if(d.day == resultDay){
			// insert the result on this day.
			d.fluxList.AddTail(r);
			return;
		}else if(resultDay < d.day){
			// insert the result at the position before this day
			measday.day = resultDay;
			measday.fluxList.AddTail(r);
			m_measurements.InsertBefore(meas_p, measday);
			return;
		}
		
		// move on to the next measurement day
		m_measurements.GetNext(meas_p);
	}
	
	// we've passed the whole list without finding anything that's larger than this day
	//	insert the result as a new measurement day in the end of the list
	measday.day = resultDay;
	measday.fluxList.AddTail(r);
	m_measurements.AddTail(measday);
	
	return;
}

/** Calculates statistics on the statistics we have here and writes 
	the results to file. */		
void CFluxStatistics::WriteFluxStat(const novac::CString &fileName){
	novac::CString str;
	FILE *f = nullptr;

	// try to open the file
	if(IsExistingFile(fileName)){
		f = fopen(fileName, "a");
		if(f == nullptr)
			return;
	}else{
		f = fopen(fileName, "w");
		if(f == nullptr)
			return;

		// write the header line
		CMeasurementDay::GetHeaderLine(str, this->m_instruments);
		fprintf(f, "%s", str);
	}

	// For each day in the list of measurement days, calculate the average flux
	//	write the number of measurements made ...
	auto p = m_measurements.GetHeadPosition();
	while(p != nullptr){
		CMeasurementDay &measDay = m_measurements.GetNext(p);
		
		measDay.GetStatistics(str, this->m_instruments);
		
		fprintf(f, "%s", str);			
	}


	// remember to close the file
	fclose(f);
}

