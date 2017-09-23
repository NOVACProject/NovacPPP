#include "StdAfx.h"
#include "postprocessing.h"

// Include synchronization classes
// #include <afxmt.h>

// the PostEvaluationController takes care of the DOAS evaluations
#include "Evaluation/PostEvaluationController.h"

// The FluxCalculator takes care of calculating the fluxes
#include "Flux/FluxCalculator.h"

// The Stratospherecalculator takes care of calculating Stratospheric VCD's
#include "Stratosphere/StratosphereCalculator.h"

// The flux CFluxStatistics takes care of the statistcal part of the fluxes
#include "Flux/FluxStatistics.h"

// This is the configuration of the network
#include "Configuration/NovacPPPConfiguration.h"

// This is the settings for how to do the procesing
#include "Configuration/UserConfiguration.h"

// We also need to read the evaluation-log files
#include "Common/EvaluationLogFileHandler.h"

#include "WindMeasurement/WindSpeedCalculator.h"

#include "Meteorology/XMLWindFileReader.h"
#include "VolcanoInfo.h"

// we want to make some statistics on the processing
#include "PostProcessingStatistics.h"

// we need to be able to download data from the FTP-server
#include "Communication/FTPServerConnection.h"

extern Configuration::CNovacPPPConfiguration        g_setup;	   // <-- The settings
extern Configuration::CUserConfiguration			g_userSettings;// <-- The settings of the user
extern CVolcanoInfo									g_volcanoes;   // <-- A list of all known volcanoes
CPostProcessingStatistics							g_processingStats; // <-- The statistics of the processing itself


// this is the working-thread that takes care of evaluating one scan
UINT EvaluateOneScan(LPVOID pParam);
// this takes care of adding the evaluated log-files to the list in an synchronized way
//	 the parameter passed in a reference to an array of strings holding the names of the 
//	 eval-log files generated
void AddResultToList(const CString &pakFileName, const CString (&evalLog)[MAX_FIT_WINDOWS], const CPlumeInScanProperty &scanProperties);
//	this retrieves the next  .pak-file from the list of files. Called from
//		'EvaluateOneScan'
int GetNextPakFileToProcess(CString &pakFileName);

CPostProcessing::CPostProcessing(void)
{
}

CPostProcessing::~CPostProcessing(void)
{
}

void CPostProcessing::DoPostProcessing_Flux(){
	CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> evalLogFiles;
	CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*> geometryResults;
	CString messageToUser, statFileName, windFileName;

	// --------------- PREPARING FOR THE PROCESSING -----------

	// Checks that the evaluation is ok and that all settings makes sense.
	if(CheckSettings()){
		ShowMessage("Exiting post processing");
		return;
	}

	// Prepare for the evaluation by reading in the reference files
	if(PrepareEvaluation()){
		ShowMessage("Exiting post processing");
		return;
	}
	
	// Prepare for the flux-calculation by reading in the wind-field
	if(ReadWindField()){
		ShowMessage("Exiting post processing");
		return;
	}
	
	// Prepare for the flux-calculation by compiling a set of pausible plume heights
	if(PreparePlumeHeights()){
		ShowMessage("Exiting post processing");
		return;
	}

	// --------------- DOING THE PROCESSING -----------

	// 1. Find all .pak files in the directory.
	CList <CString, CString&> pakFileList;
	if(g_userSettings.m_LocalDirectory.GetLength() > 3){
		CheckForSpectraInDir(g_userSettings.m_LocalDirectory, pakFileList);
	}
	if(g_userSettings.m_FTPDirectory.GetLength() > 9){
		CheckForSpectraOnFTPServer(pakFileList);
	}
	if(pakFileList.GetCount() == 0){
		ShowMessage("No spectrum files found. Exiting");
		return;
	}
	
	// Evaluate the scans. This at the same time generates a list of evaluation-log
	//	files with the evaluated results
	EvaluateScans(pakFileList, evalLogFiles);
	messageToUser.Format("%d evaluation log files accepted", evalLogFiles.GetCount());
	ShowMessage(messageToUser);	
	
	// Sort the evaluation-logs in order of increasing start-time, this to make
	//	the looking for matching files in 'CalculateGeometries' faster
	ShowMessage("Evaluation done. Sorting the evaluation log files");
	SortEvaluationLogs(evalLogFiles);
	ShowMessage("Sort done.");

	// 3. Loop through list with output text files from evaluation and calculate
	//		the geometries
	CalculateGeometries(evalLogFiles, geometryResults);

	// 4. Apply correction from AC-DC model
	while(ApplyACDCCorrections(evalLogFiles, geometryResults)){
		CalculateGeometries(evalLogFiles, geometryResults);
	}
	
	// 4.1 write the calculations to file, for later checking or other uses...
	WriteCalculatedGeometriesToFile(geometryResults);
	
	// 4.2 Insert the calculated geometries into the plume height database
	InsertCalculatedGeometriesIntoDataBase(geometryResults);
	
	// 5. Calculate the wind-speeds from the wind-speed measurements
	//		the plume heights are taken from the database
	CalculateDualBeamWindSpeeds(evalLogFiles);
	
	// 6. Calculate flux from evaluation text files
	CalculateFluxes(evalLogFiles);
	
	// 7. Write the statistics
	statFileName.Format("%s\\ProcessingStatistics.txt", g_userSettings.m_outputDirectory);
	Common::ArchiveFile(statFileName);
	g_processingStats.WriteStatToFile(statFileName);
	
	// 8. Also write the wind field that we have created to file
	windFileName.Format("%s\\GeneratedWindField.wxml", g_userSettings.m_outputDirectory);	
	Common::ArchiveFile(windFileName);
	m_windDataBase.WriteToFile(windFileName);
	
	// 9. Upload the results to the FTP-server
	if(g_userSettings.m_uploadResults){
		UploadResultsToFTP();
	}
	
	// ------------ Clean up -----------
	POSITION p = geometryResults.GetTailPosition();
	while(p != NULL){
		Geometry::CGeometryResult *g = geometryResults.GetAt(p);
		delete g;
		geometryResults.RemoveTail();
		p = geometryResults.GetTailPosition();
	}
}

/** Performs an post processing of the data in order to extract
	good stratospheric data */
void CPostProcessing::DoPostProcessing_Strat(){
	CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> evalLogFiles;
	CString messageToUser, statFileName;
	Stratosphere::CStratosphereCalculator strat;

	// --------------- PREPARING FOR THE PROCESSING -----------

	// Checks that the evaluation is ok and that all settings makes sense.
	if(CheckSettings()){
		ShowMessage("Exiting post processing");
		return;
	}

	// Prepare for the evaluation by reading in the reference files
	if(PrepareEvaluation()){
		ShowMessage("Exiting post processing");
		return;
	}	

	// --------------- DOING THE PROCESSING -----------

	// 1. Find all .pak files in the directory.
	CList <CString, CString&> pakFileList;
	if(g_userSettings.m_LocalDirectory.GetLength() > 3){
		CheckForSpectraInDir(g_userSettings.m_LocalDirectory, pakFileList);
	}
	if(g_userSettings.m_FTPDirectory.GetLength() > 9){
		CheckForSpectraOnFTPServer(pakFileList);
	}
	if(pakFileList.GetCount() == 0){
		ShowMessage("No spectrum files found. Exiting");
		return;
	}
	
	// Evaluate the scans. This at the same time generates a list of evaluation-log
	//	files with the evaluated results
	EvaluateScans(pakFileList, evalLogFiles);
	messageToUser.Format("%d evaluation log files accepted", evalLogFiles.GetCount());
	ShowMessage(messageToUser);	
	
	// Sort the evaluation-logs in order of increasing start-time, this to make
	//	the looking for matching files in 'CalculateGeometries' faster
	ShowMessage("Evaluation done. Sorting the evaluation log files");
	SortEvaluationLogs(evalLogFiles);
	ShowMessage("Sort done.");
	
	ShowMessage("Calculating stratospheric columns");
//	strat.CalculateVCDs(evalLogFiles);
	ShowMessage("Calculation Done");
	
	// 7. Write the statistics
	statFileName.Format("%s\\ProcessingStatistics.txt", g_userSettings.m_outputDirectory);
	Common::ArchiveFile(statFileName);
	g_processingStats.WriteStatToFile(statFileName);	
}

/** Performs an post processing of already evaluated data in order
	to generate plume heights and wind directions */
void CPostProcessing::DoPostProcessing_Geometry(){
	CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> evalLogFiles;
	CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*> geometryResults;
	CString messageToUser, statFileName;

}

void CPostProcessing::CheckForSpectraInDir(const CString &path, CList <CString, CString&> &fileList){
	int channel;
	CDateTime startTime;
	CString serial, fileName, userMessage;
	MEASUREMENT_MODE mode;
	HANDLE hFile;
	WIN32_FIND_DATA FindFileData;
	char fileToFind[MAX_PATH];

	userMessage.Format("Searching for .pak - files in directory %s", path);
	ShowMessage(userMessage);

	// ------------------------------------ OPTION 1 -----------------------------------------
	// ------ If we want to search for sub-directories, then search for all directories... ---
	// ---------------------------------------------------------------------------------------
	if(g_userSettings.m_includeSubDirectories_Local){
		sprintf(fileToFind, "%s\\*", path);

		// Search for the file
		hFile = FindFirstFile(fileToFind, &FindFileData);

		if(hFile != INVALID_HANDLE_VALUE){
			do{
				fileName.Format("%s\\%s", path, FindFileData.cFileName);

				if(Equals(FindFileData.cFileName, ".") || Equals(FindFileData.cFileName, ".."))
					continue;

				// if this is a directory...
				if(FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY){
					CheckForSpectraInDir(fileName, fileList);
				}else{
					// if this is a .pak - file
					if(Equals(fileName.Right(4), ".pak")){
						// if this is a incomplete scan, then don't process it.
						if(strstr(FindFileData.cFileName, "Incomplete"))
							continue;

						if(Equals(FindFileData.cFileName, "Upload.pak"))
							continue; // don't add upload.pak, it's never complete...

						// check that this file is in the time-interval that we should evaluate spectra...
						FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(FindFileData.cFileName, startTime, serial, channel, mode);

						if(startTime < g_userSettings.m_fromDate || g_userSettings.m_toDate < startTime)
							continue;

						// We've passed all the tests for the .pak-file.
						// Append the found file to the list of files to split and evaluate...
						fileList.AddTail(fileName);
					}
				}
			}while(0 != FindNextFile(hFile, &FindFileData));

			FindClose(hFile);
			return;
		}
	}
	
	// ------------------------------------ OPTION 2 -----------------------------------------
	// --------------------- Find all .pak-files in the specified directory ------------------
	// ---------------------------------------------------------------------------------------

	sprintf(fileToFind, "%s\\*.pak", path);

	// Search for the file
	hFile = FindFirstFile(fileToFind, &FindFileData);

	if(hFile == INVALID_HANDLE_VALUE)
		return; // no files found

	do{
		fileName.Format("%s\\%s", path, FindFileData.cFileName);

		// if this is a incomplete scan, then don't process it.
		if(strstr(FindFileData.cFileName, "Incomplete"))
			continue;

		if(Equals(FindFileData.cFileName, "Upload.pak"))
			continue; // don't add upload.pak, it's never complete...

		// check that this file is in the time-interval that we should evaluate spectra...
		FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(FindFileData.cFileName, startTime, serial, channel, mode);

		if(startTime < g_userSettings.m_fromDate || g_userSettings.m_toDate < startTime)
			continue;
		
		// We've passed all the tests for the .pak-file.
		// Append the found file to the list of files to split and evaluate...
		fileList.AddTail(fileName);

	}while(0 != FindNextFile(hFile, &FindFileData));

	FindClose(hFile);

	return;
}

/** Scans through the given FTP-server in search for .pak-files
	The files will be downloaded to the local computer (to the 
	temporary directory) and the returned path's will be pointing
	there.

	@param path - the directory (on the local computer) where 
		to search for files
	@param fileList - will be appended with the path's and 
	file-names of the found .pak-files */
void CPostProcessing::CheckForSpectraOnFTPServer(CList <CString, CString&> &fileList){
	Communication::CFTPServerConnection *serverDownload = new Communication::CFTPServerConnection();
	
	int ret = serverDownload->DownloadDataFromFTP(g_userSettings.m_FTPDirectory, 
		g_userSettings.m_FTPUsername, 
		g_userSettings.m_FTPPassword, 
		fileList);

	delete serverDownload;

	if(ret == 0){
		ShowMessage("Successfully downloaded all data files.");
	}else{
		ShowMessage("Error happened when downloading data from FTP.");
	}
}

CWinThread **evalThreads = NULL;
volatile int nThreadsRunning;
volatile unsigned long s_nFilesProcessed, s_nFilesToProcess;
const CList <CString, CString &> *s_pakFileList;
POSITION s_pakFileListPosition = NULL;
CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> *s_evalLogs;
CCriticalSection s_evalLogFileListCritSect; // synchronization access to the list of eval-log-files
CCriticalSection s_pakListCritSect; // synchronization access to the list of pak-files

/** Runs through the supplied list of .pak-files and evaluates each one
	using the setups found in the global settings. 
	@param pakFiles - the list of pak-files to evaluate.
	@param evalLogs - will on successful return be filled with the path's and
		filenames of each evaluation log file generated.
	*/
void CPostProcessing::EvaluateScans(const CList <CString, CString &> &pakFileList, CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> &evalLogFiles){
	s_nFilesProcessed = 0;
	s_nFilesToProcess = (long)pakFileList.GetCount();
	CString messageToUser;
	
	// share the list of eval-logs with the other functions around here
	s_evalLogs = &evalLogFiles;
	
	// share the list of pak-files with the other functions around here
	s_pakFileList			= &pakFileList;
	s_pakFileListPosition	= s_pakFileList->GetHeadPosition();

	// Keep the user informed about what we're doing
	messageToUser.Format("%ld spectrum files found. Begin evaluation", s_nFilesToProcess);
	ShowMessage(messageToUser);

	// start the threads
	evalThreads			= (CWinThread **)calloc(g_userSettings.m_maxThreadNum, sizeof(CWinThread *));
	nThreadsRunning		= 0;
	for(unsigned int k = 0; k < g_userSettings.m_maxThreadNum; ++k){
		evalThreads[k] = AfxBeginThread(EvaluateOneScan, NULL, THREAD_PRIORITY_NORMAL, 0, 0, NULL);
		Common::SetThreadName(evalThreads[k]->m_nThreadID, "EvaluateOneScan");
	}

	// make sure that all threads have time to finish before we say that we're ready
	while(nThreadsRunning > 0){
		Sleep(500);
	}

	messageToUser.Format("All %ld scans evaluated.", s_nFilesToProcess);
	ShowMessage(messageToUser);
}

UINT EvaluateOneScan(LPVOID pParam){
	// increase the count of how many threads are running
	++nThreadsRunning;

	CString evalLog[MAX_FIT_WINDOWS];
	CString fileName;
	int fitWindowIndex;
	CString messageToUser;
	CPlumeInScanProperty scanProperties[MAX_FIT_WINDOWS];

	// create a new CPostEvaluationController
	Evaluation::CPostEvaluationController *eval = new Evaluation::CPostEvaluationController();

	// while there are more .pak-files
	while(GetNextPakFileToProcess(fileName) == 0){
		// evaluate the .pak-file in all the specified fit-windows and retrieve the name of the 
		//	eval-logs. If any of the fit-windows fails then the scan is not inserted.
		for(fitWindowIndex = 0; fitWindowIndex < g_userSettings.m_nFitWindowsToUse; ++fitWindowIndex){
			if(0 != eval->EvaluateScan(fileName, g_userSettings.m_fitWindowsToUse[fitWindowIndex], &evalLog[fitWindowIndex], &scanProperties[fitWindowIndex])){
				fitWindowIndex = -2; // this is used to signal that we don't want to insert the eval-log into the list
				break;
			}
		}
		
		if(fitWindowIndex < 0){
			continue;
		}

		// If we made it this far then the measurement is ok, insert it into the list!
		AddResultToList(fileName, evalLog, scanProperties[g_userSettings.m_mainFitWindow]);
		
		// Tell the user what is happening
		messageToUser.Format(" + Inserted scan %s into list of evaluation logs", evalLog[g_userSettings.m_mainFitWindow]);
		ShowMessage(messageToUser);
	}
	
	--nThreadsRunning;
	delete eval;
	return 0;
}

// this function takes care of extracting .pak-file names
//	from the list in a synchronized way.
int GetNextPakFileToProcess(CString &pakFileName){
	CString messageToUser;

	// check if we've passed the end of the list
	if(s_pakFileListPosition == NULL)
		return 1;

	// lock access to the list
	CSingleLock singleLock(&s_pakListCritSect);
	singleLock.Lock();
	if(singleLock.IsLocked()){
		pakFileName = s_pakFileList->GetNext(s_pakFileListPosition);

		// Also tell the user the progress of our work
		if(++s_nFilesProcessed % 10 == 0){
			messageToUser.Format("Evaluated %ld scans out of %ld (%.1lf %% done)", s_nFilesProcessed, s_nFilesToProcess, 100.0 * s_nFilesProcessed / (double)s_nFilesToProcess);
			ShowMessage(messageToUser);
		}
	}

	singleLock.Unlock();
	return 0;
}

// this function takes care of adding filenames to the list
//	in a synchronized way so that no two threads access
//	the list at the same time...
void AddResultToList(const CString &pakFileName, const CString (&evalLog)[MAX_FIT_WINDOWS], const CPlumeInScanProperty &scanProperties){
	// these are not used...
	CString serial;
	int channel;
	MEASUREMENT_MODE mode;

	CSingleLock singleLock(&s_evalLogFileListCritSect);
	singleLock.Lock();
	if(singleLock.IsLocked()){
	
		// Create a new Extended scan result and add it to the end of the list
		Evaluation::CExtendedScanResult newResult;
		newResult.m_pakFile.Format(pakFileName);
		for(int fitWindowIndex = 0; fitWindowIndex < g_userSettings.m_nFitWindowsToUse; ++ fitWindowIndex){
			newResult.m_evalLogFile[fitWindowIndex].Format(evalLog[fitWindowIndex]);
			newResult.m_fitWindowName[fitWindowIndex].Format(g_userSettings.m_fitWindowsToUse[fitWindowIndex]);
		}
		FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(evalLog[0], newResult.m_startTime, serial, channel, mode);
		newResult.m_scanProperties = scanProperties;

		// store the name of the evaluation-log file generated
		s_evalLogs->AddTail(newResult);

		// update the statistics
		g_processingStats.InsertAcception(serial);
	}
	singleLock.Unlock();
}

int CPostProcessing::CheckSettings(){
	unsigned int j, k; //iterators
	CString errorMessage;
	CDateTime now;

	// Check that no instrument is duplicated in the list of instruments...
	for(j = 0; j < g_setup.m_instrumentNum; ++j){
		for(k = j + 1; k < g_setup.m_instrumentNum; ++k){
			if(Equals(g_setup.m_instrument[j].m_serial, g_setup.m_instrument[k].m_serial)){
				errorMessage.Format("The instrument %s is defined twice in setup.xml. If the instrument has two locations then define the instrument once but with two locations. Exiting post processsing.", g_setup.m_instrument[k].m_serial);
				MessageBox(NULL, errorMessage, "Error in setup.xml", MB_OK);
				return 1;
			}
		}
	}


	// Check that, for each spectrometer, there's only one fit-window defined
	//	at each instant
	for(j = 0; j < g_setup.m_instrumentNum; ++j){
		if(g_setup.m_instrument[j].m_eval.GetFitWindowNum() == 1){
			continue;
		}else{
			int ret = g_setup.m_instrument[j].m_eval.CheckSettings();
			switch(ret){
				case 0: break; // this is fine
				case 1: 
					errorMessage.Format("No fit window defined for %s. Exiting", g_setup.m_instrument[j].m_serial);
					MessageBox(NULL, errorMessage, "Error in settings", MB_OK);
					return 1;
				case 2:
					errorMessage.Format("Invalid time range found for fit window defined for %s. Exiting", g_setup.m_instrument[j].m_serial);
					MessageBox(NULL, errorMessage, "Error in settings", MB_OK);
					return 1;
				case 3:
					errorMessage.Format("At least two fit windows defined for %s have overlapping time ranges. Exiting", g_setup.m_instrument[j].m_serial);
					MessageBox(NULL, errorMessage, "Error in settings", MB_OK);
					return 1;
			}
		}
	}


	// Check that, for each spectrometer, there's only one location defined
	//	at each instant
	for(j = 0; j < g_setup.m_instrumentNum; ++j){
		if(g_setup.m_instrument[j].m_location.GetLocationNum() == 1){
			continue;
		}else{
			int ret = g_setup.m_instrument[j].m_location.CheckSettings();
			switch(ret){
				case 0: break; // this is fine
				case 1: 
					errorMessage.Format("No location defined for %s. Exiting", g_setup.m_instrument[j].m_serial);
					MessageBox(NULL, errorMessage, "Error in settings", MB_OK);
					return 1;
				case 2:
					errorMessage.Format("Invalid time range found for location defined for %s. Exiting", g_setup.m_instrument[j].m_serial);
					MessageBox(NULL, errorMessage, "Error in settings", MB_OK);
					return 1;
				case 3:
					errorMessage.Format("At least two location defined for %s have overlapping time ranges. Exiting", g_setup.m_instrument[j].m_serial);
					MessageBox(NULL, errorMessage, "Error in settings", MB_OK);
					return 1;
			}
		}
	}



	return 0;	
}

int CPostProcessing::PrepareEvaluation(){
	CDateTime fromTime, toTime; //  these are not used but must be passed onto GetFitWindow...
	CString errorMessage, fileName;

	// this is true if we failed to prepare the evaluation...
	bool failure = false;

	// Loop through each of the configured instruments
	for(unsigned int instrumentIndex = 0; instrumentIndex < g_setup.m_instrumentNum; ++instrumentIndex){
		
		// For each instrument, loop through the fit-windows that are defined
		unsigned long fitWindowNum = g_setup.m_instrument[instrumentIndex].m_eval.GetFitWindowNum();
		for(unsigned int fitWindowIndex = 0; fitWindowIndex < fitWindowNum; ++fitWindowIndex){
			Evaluation::CFitWindow window;
			
			// get the fit window
			g_setup.m_instrument[instrumentIndex].m_eval.GetFitWindow(fitWindowIndex, window, fromTime, toTime);
			
			// For each reference in the fit-window, read it in and make sure that it exists...
			for(int referenceIndex = 0; referenceIndex < window.nRef; ++referenceIndex){
				if(!IsExistingFile(window.ref[referenceIndex].m_path)){
					// the file does not exist, try to change it to include the path of the configuration-directory...
					fileName.Format("%s\\configuration\\%s", m_exePath, window.ref[referenceIndex].m_path);
						if(IsExistingFile(fileName)){
							window.ref[referenceIndex].m_path.Format(fileName);
						}else{
							errorMessage.Format("Cannot read reference file %s", window.ref[referenceIndex].m_path);
							ShowMessage(errorMessage);
							failure = true;
							continue;
						}
				}

				// Read in the cross section
				if(window.ref[referenceIndex].ReadCrossSectionDataFromFile()){
					errorMessage.Format("Failed to read cross section file: %s", window.ref[referenceIndex].m_path);
					ShowMessage(errorMessage);
					failure = true;
					continue;
				}
				
				// If we are supposed to high-pass filter the spectra then
				//	we should also high-pass filter the cross-sections
				if(window.fitType == Evaluation::FIT_HP_DIV || window.fitType == Evaluation::FIT_HP_SUB){
					if(window.ref[referenceIndex].m_isFiltered == false){
						if(Equals(window.ref[referenceIndex].m_specieName, "ring")){
							window.ref[referenceIndex].m_data->HighPassFilter_Ring();
						}else{
							window.ref[referenceIndex].m_data->HighPassFilter();
						}
					}else{
						// The filtered cross sections are scaled to the unit of ppmm
						//	rescale them to molecules/cm2 as all other cross sections
						window.ref[referenceIndex].m_data->Multiply(1.0/2.5e15);
					}
				}// endif
			}
			
			// If the window also contains a fraunhofer-reference then read it too.
			if(window.fraunhoferRef.m_path.GetLength() > 4){
				if(!IsExistingFile(window.fraunhoferRef.m_path)){
					// the file does not exist, try to change it to include the path of the configuration-directory...
					fileName.Format("%s\\configuration\\%s", m_exePath, window.fraunhoferRef.m_path);
						if(IsExistingFile(fileName)){
							window.fraunhoferRef.m_path.Format(fileName);
						}else{
							errorMessage.Format("Cannot read reference file %s", window.fraunhoferRef.m_path);
							ShowMessage(errorMessage);
							failure = true;
							continue;
						}
				}

				if(window.fraunhoferRef.ReadCrossSectionDataFromFile()){
					errorMessage.Format("Failed to read fraunhofer-reference file: %s", window.fraunhoferRef.m_path);
					ShowMessage(errorMessage);
					failure = true;
					continue;
				}
				if(window.fitType == Evaluation::FIT_HP_DIV || window.fitType == Evaluation::FIT_HP_SUB){
					window.fraunhoferRef.m_data->HighPassFilter_Ring();
				}else{
					window.fraunhoferRef.m_data->Log();
				}
			}
			
			// If we've made it this far, then we've managed to read in all the references. Now
			//	store the data in g_setup
			g_setup.m_instrument[instrumentIndex].m_eval.SetFitWindow(fitWindowIndex, window, &fromTime, &toTime);
		}
	}

	if(failure){
		return 1;
	}else{
		return 0;
	}
}

/** Prepares for the flux calculations by reading in the relevant
	wind-field file.
	@return 0 on success, otherwise non-zero */
int CPostProcessing::ReadWindField(){
	CString name1, name2, name3, path1, path2, path3, messageToUser;
	Common common;
	FileHandler::CXMLWindFileReader reader;
	common.GetExePath();

	if(g_userSettings.m_windFieldFileOption == 0){

		// If the user has given a file-name, then try to use that one
		if(g_userSettings.m_windFieldFile.GetLength() > 3){
			messageToUser.Format("Reading wind field from file: %s", g_userSettings.m_windFieldFile);
			ShowMessage(messageToUser);
			
			if(reader.ReadWindFile(g_userSettings.m_windFieldFile, m_windDataBase)){
				messageToUser.Format("Failed to parse wind field file: %s", g_userSettings.m_windFieldFile);
				ShowMessage(messageToUser);
				return 1;
			}else{
				messageToUser.Format("Parsed %s containing %d wind data items", g_userSettings.m_windFieldFile, m_windDataBase.GetDataBaseSize());
				ShowMessage(messageToUser);
			
				name1.Format("%sParsedWindField.wxml", common.m_exePath);
				m_windDataBase.WriteToFile(name1);
				return 0;
			}
		}else{
			// Get the name of the volcano that we are about to process...
			//  there are two options, either the full name or the simple name
			g_volcanoes.GetVolcanoName(g_userSettings.m_volcano, name1);
			g_volcanoes.GetSimpleVolcanoName(g_userSettings.m_volcano, name2);
			g_volcanoes.GetVolcanoCode(g_userSettings.m_volcano, name3);

			// Get the path to the executable, so that we know where to start looking
			path1.Format("%s\\configuration\\%s.wxml", common.m_exePath, name1);
			path2.Format("%s\\configuration\\%s.wxml", common.m_exePath, name2);
			path3.Format("%s\\configuration\\%s.wxml", common.m_exePath, name3);

			// check which of the files exists
			if(IsExistingFile(path1)){
				messageToUser.Format("Reading wind field from file: %s", path1);
				ShowMessage(messageToUser);
				
				if(reader.ReadWindFile(path1, m_windDataBase)){
					messageToUser.Format("Failed to parse wind field file: %s", path1);
					ShowMessage(messageToUser);
					return 1;
				}else{
					return 0;
				}
				
			}else if(IsExistingFile(path2)){
				messageToUser.Format("Reading wind field from file: %s", path2);
				ShowMessage(messageToUser);

				if(reader.ReadWindFile(path2, m_windDataBase)){
					messageToUser.Format("Failed to parse wind field file: %s", path2);
					ShowMessage(messageToUser);
					return 1;
				}else{
					return 0;
				}
			}else if(IsExistingFile(path3)){
				messageToUser.Format("Reading wind field from file: %s", path3);
				ShowMessage(messageToUser);

				if(reader.ReadWindFile(path3, m_windDataBase)){
					messageToUser.Format("Failed to parse wind field file: %s", path3);
					ShowMessage(messageToUser);
					return 1;
				}else{
					return 0;
				}	
			}else{
				messageToUser.Format("Cannot find wind field file. Attempted to read: %s, %s and %s", path1, path2, path3);
				ShowMessage(messageToUser);
				return 1;
			}
		}
		return 1; // fail, should never get here
	}
	
	// If the user has specified a directory of files...
	if(g_userSettings.m_windFieldFileOption == 1){
		if(reader.ReadWindDirectory(g_userSettings.m_windFieldFile, m_windDataBase, &g_userSettings.m_fromDate, &g_userSettings.m_toDate)){
			return 1;
		}
		return 0;
	}
	
	// should never get to this point!
	return 1;
}

/** Prepares for the flux calculation by setting up a reasonable
	set of plume heights. This could also read in a set from file...?
	@return 0 on success, otherwese non-zero */
int CPostProcessing::PreparePlumeHeights(){
	// we need to construct a default plume height to use, if there's nothing else...
	Geometry::CPlumeHeight plumeHeight;
	plumeHeight.m_plumeAltitude			= g_volcanoes.GetPeakAltitude(g_userSettings.m_volcano);
	plumeHeight.m_plumeAltitudeSource	= Meteorology::MET_DEFAULT;
	plumeHeight.m_validFrom				= CDateTime(0,0,0,0,0,0);
	plumeHeight.m_validTo				= CDateTime(9999,12,31,23,59,59);
	
	// the estimated plume height is half of the altitude difference between the highest
	//	instrument for this volcano and the volcano altitude
	double maxInstrumentAltitude = -1e6;
	Configuration::CInstrumentLocation location;
	CString volcanoName;
	g_volcanoes.GetVolcanoName(g_userSettings.m_volcano, volcanoName);
	for(unsigned int k = 0; k < g_setup.m_instrumentNum; ++k){
		unsigned long N = g_setup.m_instrument[k].m_location.GetLocationNum();
		for(unsigned int j = 0; j < N; ++j){
			g_setup.m_instrument[k].m_location.GetLocation(j, location);
			if(Equals(volcanoName, location.m_volcano)){
				maxInstrumentAltitude = max(maxInstrumentAltitude, location.m_altitude);
			}
		}
	}
	if(maxInstrumentAltitude > 0){
		plumeHeight.m_plumeAltitudeError	= fabs(g_volcanoes.GetPeakAltitude(g_userSettings.m_volcano) - maxInstrumentAltitude) / 2.0;
	}else{
		plumeHeight.m_plumeAltitudeError	= g_volcanoes.GetPeakAltitude(g_userSettings.m_volcano) / 2.0;
	}

	m_plumeDataBase.InsertPlumeHeight(plumeHeight);	

	return 0;
}

/** Runs through the supplied list of evaluation - logs and performs
	geometry calculations on the ones which does match. The results
	are returned in the list geometryResults. 
	@param evalLogs - list of CStrings, each holding the full path and filename 
		of an evaluation-log file that should be considered for geometrical 
		calculations
	@param geometryResults - will on successfull return be filled with the
		calculated plume heights and wind-directions.
	*/
void CPostProcessing::CalculateGeometries(const CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult&> &evalLogFiles, CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*> &geometryResults){
	CString serial1, serial2, messageToUser;
	CDateTime startTime1, startTime2;
	MEASUREMENT_MODE measMode1, measMode2;
	int channel;
	unsigned long nFilesChecked1 = 0; // this is for debugging purposes...
	unsigned long nFilesChecked2 = 0; // this is for debugging purposes...
	unsigned long nCalculationsMade = 0; // this is for debugging purposes...
	unsigned long nTooLongdistance = 0; // this is for debugging purposes...
	unsigned long nTooLargeAbsoluteError = 0; // this is for debugging purposes...
	unsigned long nTooLargeRelativeError = 0; // this is for debugging purposes...
	Configuration::CInstrumentLocation location[2];

	// Tell the user what's happening
	ShowMessage("Begin to calculate plume heights from scans");

	// Loop through list with output text files from evaluation and apply geometrical corrections
	POSITION pos1 = evalLogFiles.GetHeadPosition();
	while(pos1 != NULL){
		const CString &evalLog1				= evalLogFiles.GetAt(pos1).m_evalLogFile[g_userSettings.m_mainFitWindow];
		const CPlumeInScanProperty &plume1	= evalLogFiles.GetNext(pos1).m_scanProperties;
		
		++nFilesChecked1; // for debugging...
		
		// if this is the last file in the list, then
		//	quit. There's nothing more to compare to...
		if(pos1 == NULL){
			break;
		}
		
		// if this scan does not see a large enough portion of the plume, then ignore it...
		if(plume1.m_completeness < g_userSettings.m_calcGeometry_CompletenessLimit)
			continue;

		//  Get the information about evaluation log file #1
		FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(evalLog1, startTime1, serial1, channel, measMode1);
		
		// If this is not a flux-measurement, then there's no use in trying to use it...
		if(measMode1 != MODE_FLUX)
			continue;

		// try to combine this evaluation-log file with every other eval-log
		//  use the fact that the list of eval-logs is sorted by increasing start-time
		//  thus we start at the eval-log next after this one and compare with all
		//  eval-logs until the difference in start-time is too big.
		POSITION pos2 = pos1;
		evalLogFiles.GetNext(pos2);
		bool successfullyCombined = false; // this is true if evalLog1 was combined with (at least one) other eval-log to make a geomery calculation.
		while(pos2 != NULL){
			const CString &evalLog2				= evalLogFiles.GetAt(pos2).m_evalLogFile[g_userSettings.m_mainFitWindow];
			const CPlumeInScanProperty &plume2	= evalLogFiles.GetNext(pos2).m_scanProperties;

			++nFilesChecked2; // for debugging...

			// if this scan does not see a large enough portion of the plume, then ignore it...
			if(plume2.m_completeness < g_userSettings.m_calcGeometry_CompletenessLimit)
				continue;

			//  Get the information about evaluation log file # 2
			FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(evalLog2, startTime2, serial2, channel, measMode2);
			
			// The time elapsed between the two measurements must not be more than 
			//	the user defined time-limit (in seconds)
			double timeDifference = fabs(CDateTime::Difference(startTime1, startTime2));
			if(timeDifference > g_userSettings.m_calcGeometry_MaxTimeDifference){
				pos2 = NULL;
				continue;
			}

			// If this is not a flux-measurement, then there's no use in trying to use it...
			if(measMode2 != MODE_FLUX)
				continue;

			// the serials must be different (i.e. the two measurements must be
			//  from two different instruments)
			if(Equals(serial1, serial2)){
				continue;
			}
			
			// Get the locations of the two instruments
			if(g_setup.GetInstrumentLocation(serial1, startTime1, location[0]))
				continue;
			if(g_setup.GetInstrumentLocation(serial2, startTime2, location[1]))
				continue;

			// make sure that the distance between the instruments is not too long....
			double instrumentDistance = Common::GPSDistance(location[0].m_latitude, location[0].m_longitude, location[1].m_latitude, location[1].m_longitude);
			if(instrumentDistance < g_userSettings.m_calcGeometry_MinDistance || instrumentDistance > g_userSettings.m_calcGeometry_MaxDistance){
				++nTooLongdistance;
				continue;
			}

			// count the number of times we calculate a result, for improving the software...
			++nCalculationsMade;

			// If the files have passed these tests then make a geometry-calculation
			Geometry::CGeometryResult *result = new Geometry::CGeometryResult();
			if(Geometry::CGeometryCalculator::CalculateGeometry(plume1, startTime1, plume2, startTime2, location, *result)){

				// Check the quality of the measurement before we insert it...
				if(result->m_plumeAltitudeError > g_userSettings.m_calcGeometry_MaxPlumeAltError){
					++nTooLargeAbsoluteError;
					delete result; // too bad, continue.
				}else if((result->m_plumeAltitudeError > 0.5*result->m_plumeAltitude) || (result->m_windDirectionError > g_userSettings.m_calcGeometry_MaxWindDirectionError)){
					++nTooLargeRelativeError;
					delete result; // too bad, continue.
				}else if(result->m_windDirectionError > g_userSettings.m_calcGeometry_MaxWindDirectionError){
					delete result; // too bad, continue.
				}else{
					// remember which instruments were used
					result->m_instr1.Format(serial1);
					result->m_instr2.Format(serial2);

					geometryResults.AddTail(result);
					
					messageToUser.Format(" + Calculated a plume altitude of %.0lf � %.0lf meters and wind direction of %.0lf � %.0lf degrees by combining measurements from %s and %s",
						result->m_plumeAltitude, result->m_plumeAltitudeError, result->m_windDirection, result->m_windDirectionError, serial1, serial2);
					ShowMessage(messageToUser);
					
					successfullyCombined = true;
				}
			}else{
				// something went wrong... delete the 'info'
				delete result;
			}
		} // end while(pos2 != NULL)

		// if it was not possible to combine this scan with any other to generate an
		//	estimated plume height and wind direction we might still be able to use it to calculate
		//	a wind direction given the plume height at the time of the measurement.
		if(!successfullyCombined){
			Meteorology::CWindField windField;
			Geometry::CPlumeHeight plumeHeight;

			// Get the location of the instrument
			if(g_setup.GetInstrumentLocation(serial1, startTime1, location[0]))
				continue;

			Geometry::CGeometryResult *result = new Geometry::CGeometryResult();
			
			// Get the altitude of the plume at this moment. First look into the
			//	general database. Then have a look in the list of geometry-results
			//	that we just generated to see if there's anything better there...
			m_plumeDataBase.GetPlumeHeight(startTime1, plumeHeight);
			POSITION gp = geometryResults.GetTailPosition();
			while(gp != NULL){
				const Geometry::CGeometryResult *oldResult = geometryResults.GetPrev(gp);
				if(fabs(CDateTime::Difference(oldResult->m_averageStartTime, startTime1)) < g_userSettings.m_calcGeometryValidTime){
					if((oldResult->m_plumeAltitudeError < plumeHeight.m_plumeAltitudeError) && (oldResult->m_plumeAltitude > NOT_A_NUMBER)){
						plumeHeight.m_plumeAltitude			= oldResult->m_plumeAltitude;
						plumeHeight.m_plumeAltitudeError	= oldResult->m_plumeAltitudeError;
						plumeHeight.m_plumeAltitudeSource	= oldResult->m_calculationType;
					}
				}
			}
			
			// Try to calculate the wind-direction
			if(Geometry::CGeometryCalculator::CalculateWindDirection(evalLog1, 0, plumeHeight, location[0], *result)){
				// Success!!
				result->m_instr1.Format(serial1);
				geometryResults.AddTail(result);
				
				// tell the user			
				messageToUser.Format(" + Calculated a wind direction of %.0lf � %.0lf degrees from a scan by instrument %s",
					result->m_windDirection, result->m_windDirectionError, serial1);
				ShowMessage(messageToUser);
			}else{
				delete result;
				continue;
			}
		}
	} // end while(pos1 != NULL)
	
	// Tell the user what we have done
	if(geometryResults.GetCount() == 0){
		ShowMessage("No plume heights could be calculated");
	}else{
		messageToUser.Format("Done calculating geometries. Plume height calculated on %d occasions", geometryResults.GetCount());
		ShowMessage(messageToUser);
	}
	messageToUser.Format("nFilesChecked1 = %ld, nFilesChecked2 = %ld, nCalculationsMade = %ld", nFilesChecked1, nFilesChecked2, nCalculationsMade);
	ShowMessage(messageToUser);	
}

/** Runs through the supplied list of evaluation - logs and performs
	AC-DC corrections on the derived columns. The results are not
	returned, instead the files are re-written with the updated 
	column values. 
	@param evalLogs - list of CStrings, each holding the full path and filename 
		of an evaluation-log file that should be considered for geometrical 
		calculations
	@param geometryResults - list of calculated geometrical results. These will
		be used to apply the radiative corrections to the columns. For scans when 
		no plume-height can be taken from the calculations - a default plume height equal
		to the altitude of the summit of the volcano will be used.
	@return - true if so large changes are made that the geometries would need to
		be re-calculated. Otherwise false.
	*/	
bool CPostProcessing::ApplyACDCCorrections(const CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> &evalLogs, const CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*> &geometryResults){

	ShowMessage("Applying ACDC corrections - This is not yet implemented!!");

	return false;
}

/** Runs through the supplied list of evaluation-logs and 
	calculates the flux for each scan. The resulting fluxes are written
	to a flux-log file in the output directory. 
	@param evalLogs - list of CStrings, each holding the full path and filename 
		of an evaluation-log file that should be considered for geometrical 
		calculations
	@param geometryResults - list of calculated geometrical results. These will
		be used to calculate the fluxes. For scans when no plume-height can be 
		taken from the calculations - a default plume height equal to the 
		altitude of the summit of the volcano will be used.
	*/
void CPostProcessing::CalculateFluxes(const CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> &evalLogFiles){
	CDateTime scanStartTime;
	CString serial, messageToUser;
	Geometry::CPlumeHeight plumeHeight; // the altitude of the plume, in meters above sea level
	MEASUREMENT_MODE measMode;
	int channel;
	Flux::CFluxStatistics stat;

	// we keep the calculated fluxes in a list
	CList <Flux::CFluxResult, Flux::CFluxResult &> calculatedFluxes;

	// Initiate the flux-calculator
	Flux::CFluxCalculator *fluxCalc = new Flux::CFluxCalculator();

	// Loop through the list of evaluation log files. For each of them, find
	//	the best available wind-speed, wind-direction and plume height and
	//	calculate the flux.
	POSITION pos = evalLogFiles.GetHeadPosition();
	while(pos != NULL){
		// Get the name of this eval-log
		const CString &evalLog				= evalLogFiles.GetAt(pos).m_evalLogFile[g_userSettings.m_mainFitWindow];
		const CPlumeInScanProperty &plume	= evalLogFiles.GetNext(pos).m_scanProperties;

		// if the completeness is too low then ignore this scan.
		if(plume.m_completeness < (g_userSettings.m_completenessLimitFlux + 0.01)){
			messageToUser.Format(" - Scan %s has completeness = %.2lf which is less than limit of %.2lf. Rejected!", evalLog, plume.m_completeness, g_userSettings.m_completenessLimitFlux);
			ShowMessage(messageToUser);
			continue;
		}

		// Extract the date and time of day when the measurement was made
		FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(evalLog, scanStartTime, serial, channel, measMode);
		
		// If this is not a flux-measurement, then there's no point in calculating any flux for it
		if(measMode != MODE_FLUX)
			continue;
		
		// Extract a plume height at this time of day
		m_plumeDataBase.GetPlumeHeight(scanStartTime, plumeHeight);

		// tell the user
		messageToUser.Format("Calculating flux for measurement %s", evalLog);
		ShowMessage(messageToUser);

		// Calculate the flux. This also takes care of writing
		//	the results to file
		Flux::CFluxResult fluxResult;
		if(0 == fluxCalc->CalculateFlux(evalLog, m_windDataBase, plumeHeight, fluxResult)){
			calculatedFluxes.AddTail(fluxResult);
		}
	}
	
	// Now we can write the final fluxes to file
	ShowMessage("Writing flux log");
	WriteFluxResult_XML(calculatedFluxes);
	WriteFluxResult_Txt(calculatedFluxes);
	
	// Also write the statistics for the flux
	ShowMessage("Writing flux statistics");
	stat.AttachFluxList(calculatedFluxes);
	stat.WriteFluxStat(_T(g_userSettings.m_outputDirectory + "\\FluxStatistics.txt"));
	
	// clean up...
	delete fluxCalc;
}

void CPostProcessing::WriteFluxResult_XML(const CList <Flux::CFluxResult, Flux::CFluxResult &> &calculatedFluxes){
	CString fluxLogFile, styleFile, wsSrc, wdSrc, phSrc, typeStr;
	CDateTime now;

	// get the current time
	now.SetToNow();
		
	// the name and path of the final flux-log file
	fluxLogFile.Format("%s\\FluxLog.xml", g_userSettings.m_outputDirectory);

	// Check if there's already a file like this, then archive it...
	Common::ArchiveFile(fluxLogFile);

	// Try to open the file
	FILE *f = fopen(fluxLogFile, "w");
	if(f == NULL){
		ShowMessage("Could not open flux log file for writing. Writing of results failed. ");
		return;
	}
	
	// Write the header and the starting comments
	fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
	fprintf(f, "<?xml-stylesheet type=\"text/xsl\" href=\"fluxresult.xsl\"?>\n");
	fprintf(f, "<!-- This is result of the flux calculations the NOVAC Post Processing Program -->\n");
	fprintf(f, "<!-- File generated on %04d.%02d.%02d at %02d:%02d:%02d -->\n\n", now.year, now.month, now.day, now.hour, now.minute, now.second);

	fprintf(f, "<NovacPPPFluxResults>\n");
	POSITION pos = calculatedFluxes.GetHeadPosition();
	while(pos != NULL){
		// Get the next flux result in the list
		const Flux::CFluxResult &fluxResult = calculatedFluxes.GetNext(pos);
		
		// extract the sources of information about wind-speed, wind-direction and plume-height
		fluxResult.m_windField.GetWindSpeedSource(wsSrc);
		fluxResult.m_windField.GetWindDirectionSource(wdSrc);
		Meteorology::MetSourceToString(fluxResult.m_plumeHeight.m_plumeAltitudeSource, phSrc);
		
		// write a <flux> section
		fprintf(f, "\t<flux>\n");
		
		fprintf(f, "\t\t<startTime>%04d.%02d.%02dT%02d:%02d:%02d</startTime>\n", 
			fluxResult.m_startTime.year, fluxResult.m_startTime.month, fluxResult.m_startTime.day, 
			fluxResult.m_startTime.hour, fluxResult.m_startTime.minute,fluxResult.m_startTime.second);
		fprintf(f, "\t\t<stopTime>%04d.%02d.%02dT%02d:%02d:%02d</stopTime>\n", 
			fluxResult.m_stopTime.year, fluxResult.m_stopTime.month, fluxResult.m_stopTime.day, 
			fluxResult.m_stopTime.hour, fluxResult.m_stopTime.minute,fluxResult.m_stopTime.second);

		fprintf(f, "\t\t<serial>%s</serial>\n",							fluxResult.m_instrument);

		// extract the instrument type
		if(fluxResult.m_instrumentType == INSTR_HEIDELBERG){
			typeStr.Format("heidelberg");
		}else{
			typeStr.Format("gothenburg");
		}		
		fprintf(f, "\t\t<instrumentType>%s</instrumentType>\n",			typeStr);

		fprintf(f, "\t\t<value>%.2lf</value>\n",						fluxResult.m_flux);

		// The judged quality of the calculated flux
		if(fluxResult.m_fluxQualityFlag == FLUX_QUALITY_GREEN){
			fprintf(f, "\t\t<Quality>g</Quality>\n");
		}else if(fluxResult.m_fluxQualityFlag == FLUX_QUALITY_YELLOW){
			fprintf(f, "\t\t<Quality>y</Quality>\n");
		}else{
			fprintf(f, "\t\t<Quality>r</Quality>\n");
		}		

		// the errors
		fprintf(f, "\t\t<FluxError_Wind_kgs>%.2lf</FluxError_Wind_kgs>\n",					fluxResult.m_fluxError_Wind);
		fprintf(f, "\t\t<FluxError_PlumeHeight_kgs>%.2lf</FluxError_PlumeHeight_kgs>\n",	fluxResult.m_fluxError_PlumeHeight);

		// the wind speed
		fprintf(f, "\t\t<windspeed>%.2lf</windspeed>\n",				fluxResult.m_windField.GetWindSpeed());
		fprintf(f, "\t\t<windspeedError>%.2lf</windspeedError>\n",		fluxResult.m_windField.GetWindSpeedError());
		fprintf(f, "\t\t<windspeedSource>%s</windspeedSource>\n",		wsSrc);

		// the wind direction
		fprintf(f, "\t\t<winddirection>%.2lf</winddirection>\n",			fluxResult.m_windField.GetWindDirection());
		fprintf(f, "\t\t<winddirectionError>%.2lf</winddirectionError>\n", fluxResult.m_windField.GetWindDirectionError());
		fprintf(f, "\t\t<winddirectionSource>%s</winddirectionSource>\n", wdSrc);

		// the plume height
		fprintf(f, "\t\t<plumeheight>%.2lf</plumeheight>\n",			fluxResult.m_plumeHeight.m_plumeAltitude);
		fprintf(f, "\t\t<plumeheightError>%.2lf</plumeheightError>\n",	fluxResult.m_plumeHeight.m_plumeAltitudeError);
		fprintf(f, "\t\t<plumeheightSource>%s</plumeheightSource>\n",	phSrc);

		// some additional information about the scan
		fprintf(f, "\t\t<Compass>%.1lf<Compass>\n",						fluxResult.m_compass);
		fprintf(f, "\t\t<ConeAngle>%.1lf<ConeAngle>\n",					fluxResult.m_coneAngle);
		fprintf(f, "\t\t<Tilt>%.1lf<Tilt>\n",							fluxResult.m_tilt);
		fprintf(f, "\t\t<nSpectra>%d<nSpectra>\n",						fluxResult.m_numGoodSpectra);
		fprintf(f, "\t\t<PlumeCentre_1>%.1lf<PlumeCentre_1>\n",			fluxResult.m_plumeCentre[0]);
		fprintf(f, "\t\t<PlumeCentre_2>%.1lf<PlumeCentre_2>\n",			fluxResult.m_plumeCentre[1]);
		fprintf(f, "\t\t<PlumeCompleteness>%.2lf<PlumeCompleteness>\n",	fluxResult.m_completeness);
		fprintf(f, "\t\t<ScanOffset>%.1e<ScanOffset>\n",				fluxResult.m_scanOffset);

		fprintf(f, "\t</flux>\n");
	}

	fprintf(f, "</NovacPPPFluxResults>\n");

	// remember to close the file
	fclose(f);
	
	// ------------- we also need an xslt - file to display the output -----------------
	styleFile.Format("%s\\fluxresult.xsl", g_userSettings.m_outputDirectory);

	// Try to open the file
	f = fopen(styleFile, "w");
	if(f == NULL){
		return;
	}
	fprintf(f, "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");

	fprintf(f, "<html xsl:version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\" xmlns=\"http://www.w3.org/1999/xhtml\">\n");
	fprintf(f, "<body style=\"font-family:Arial;font-size:12pt;background-color:#EEEEEE\">\n");
  	fprintf(f, "\t<div style=\"background-color:white;color:black;padding:4px\">\n");
	fprintf(f, "\t\t<span style=\"font-weight:bold\">Result of flux calculation</span>\n");
  	fprintf(f, "\t</div>\n");

	fprintf(f, "\t<xsl:for-each select=\"NovacPPPFluxResults/flux\">\n");
	fprintf(f, "\t<div style=\"background-color:white;color:teal;padding:4px\">\n");
	fprintf(f, "\t\t<span style=\"font-weight:bold\">Measurement from <xsl:value-of select=\"startTime\"/> to <xsl:value-of select=\"stopTime\"/> </span>\n");
	fprintf(f, "\t- <xsl:value-of select=\"value\"/> kg/s\n");
	fprintf(f, "\t</div>\n");
	fprintf(f, "\t<div style=\"margin-left:20px;margin-bottom:1em;font-size:10pt\">\n");

	fprintf(f, "\t\t<xsl:value-of select=\"description\"/>\n");
	fprintf(f, "\t\t<span style=\"font-style:italic\">\n");
	fprintf(f, "\t\t\tMade by <xsl:value-of select=\"serial\"/>\n");
	fprintf(f, "\t\t</span>\n");
	fprintf(f, "\t</div>\n");
	fprintf(f, "\t</xsl:for-each>\n");
	fprintf(f, "</body>\n");
	fprintf(f, "</html>\n");

	fclose(f);	
}

void CPostProcessing::WriteFluxResult_Txt(const CList <Flux::CFluxResult, Flux::CFluxResult &> &calculatedFluxes){
	CString fluxLogFile, wsSrc, wdSrc, phSrc, typeStr;
	CDateTime now;

	// get the current time
	now.SetToNow();
		
	// the name and path of the final flux-log file
	fluxLogFile.Format("%s\\FluxLog.txt", g_userSettings.m_outputDirectory);

	// Try to open the file
	if(IsExistingFile(fluxLogFile)){
		Common::ArchiveFile(fluxLogFile);
	}

	FILE *f = fopen(fluxLogFile, "w");
	if(f == NULL){
		ShowMessage("Could not open flux log file for writing. Writing of results failed. ");
		return;
	}
	
	// Write the header and the starting comments
	fprintf(f, "# This is result of the flux calculations the NOVAC Post Processing Program \n");
	fprintf(f, "#   File generated on %04d.%02d.%02d at %02d:%02d:%02d \n\n", now.year, now.month, now.day, now.hour, now.minute, now.second);

	fprintf(f, "#StartTime\tStopTime\tSerial\tInstrumentType\tFlux_kgs\tFluxQuality\tFluxError_Wind_kgs\tFluxError_PlumeHeight_kgs\tWindSpeed_ms\tWindSpeedErr_ms\tWindSpeedSrc\tWindDir_deg\tWindDirErr_deg\tWindDirSrc\tPlumeHeight_m\tPlumeHeightErr_m\tPlumeHeightSrc\t");
	fprintf(f, "Compass\tConeAngle\tTilt\tnSpectra\tPlumeCentre_1\tPlumeCentre_2\tPlumeCompleteness\tScanOffset\n");

	POSITION pos = calculatedFluxes.GetHeadPosition();
	while(pos != NULL){
		// Get the next flux result in the list
		const Flux::CFluxResult &fluxResult = calculatedFluxes.GetNext(pos);
		
		// extract the instrument type
		if(fluxResult.m_instrumentType == INSTR_HEIDELBERG){
			typeStr.Format("heidelberg");
		}else{
			typeStr.Format("gothenburg");
		}
		
		// extract the sources of information about wind-speed, wind-direction and plume-height
		fluxResult.m_windField.GetWindSpeedSource(wsSrc);
		fluxResult.m_windField.GetWindDirectionSource(wdSrc);
		Meteorology::MetSourceToString(fluxResult.m_plumeHeight.m_plumeAltitudeSource, phSrc);
		
		// write the date and time when the measurement started and ended
		fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t", 
			fluxResult.m_startTime.year, fluxResult.m_startTime.month, fluxResult.m_startTime.day, 
			fluxResult.m_startTime.hour, fluxResult.m_startTime.minute,fluxResult.m_startTime.second);
		fprintf(f, "%04d.%02d.%02dT%02d:%02d:%02d\t", 
			fluxResult.m_stopTime.year, fluxResult.m_stopTime.month, fluxResult.m_stopTime.day, 
			fluxResult.m_stopTime.hour, fluxResult.m_stopTime.minute,fluxResult.m_stopTime.second);

		// the type of instrument and the serial-number
		fprintf(f, "%s\t",		fluxResult.m_instrument);
		fprintf(f, "%s\t",		typeStr);

		// The actual flux!!!
		fprintf(f, "%.2lf\t",	fluxResult.m_flux);

		// The judged quality of the calculated flux
		if(fluxResult.m_fluxQualityFlag == FLUX_QUALITY_GREEN){
			fprintf(f, "g\t");
		}else if(fluxResult.m_fluxQualityFlag == FLUX_QUALITY_YELLOW){
			fprintf(f, "y\t");
		}else{
			fprintf(f, "r\t");
		}
		
		// the errors
		fprintf(f, "%.2lf\t",	fluxResult.m_fluxError_Wind);
		fprintf(f, "%.2lf\t",	fluxResult.m_fluxError_PlumeHeight);
		
		// the wind speed
		fprintf(f, "%.2lf\t",	fluxResult.m_windField.GetWindSpeed());
		fprintf(f, "%.2lf\t",	fluxResult.m_windField.GetWindSpeedError());
		fprintf(f, "%s\t",		wsSrc);

		// the wind direction
		fprintf(f, "%.2lf\t",	fluxResult.m_windField.GetWindDirection());
		fprintf(f, "%.2lf\t",	fluxResult.m_windField.GetWindDirectionError());
		fprintf(f, "%s\t",		wdSrc);

		// the plume height
		fprintf(f, "%.2lf\t",	fluxResult.m_plumeHeight.m_plumeAltitude);
		fprintf(f, "%.2lf\t",	fluxResult.m_plumeHeight.m_plumeAltitudeError);
		fprintf(f, "%s\t",		phSrc);
		
		// write additional information about the scan
		fprintf(f, "%.1lf\t",	fluxResult.m_compass);
		fprintf(f, "%.1lf\t",	fluxResult.m_coneAngle);
		fprintf(f, "%.1lf\t",	fluxResult.m_tilt);
		fprintf(f, "%d\t",		fluxResult.m_numGoodSpectra);
		fprintf(f, "%.1lf\t",	fluxResult.m_plumeCentre[0]);
		fprintf(f, "%.1lf\t",	fluxResult.m_plumeCentre[1]);
		fprintf(f, "%.2lf\t",	fluxResult.m_completeness);
		fprintf(f, "%.1e\n",	fluxResult.m_scanOffset);
		
	}

	// remember to close the file
	fclose(f);
}

void CPostProcessing::WriteCalculatedGeometriesToFile(const CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*> &geometryResults){
	if(geometryResults.GetCount() == 0)
		return; // nothing to write...

	FILE *f = NULL;
	CString geomLogFile;
	geomLogFile.Format("%s\\GeometryLog.txt", g_userSettings.m_outputDirectory);

	if(IsExistingFile(geomLogFile)){
		f = fopen(geomLogFile, "a");
		if(f == NULL){
			ShowMessage("Could not open geometry log file for writing. Writing of results failed. ");
			return;
		}
	}else{
		f = fopen(geomLogFile, "w");
		if(f == NULL){
			ShowMessage("Could not open geometry log file for writing. Writing of results failed. ");
			return;
		}
		fprintf(f, "Date\tTime\tDifferenceInStartTime_minutes\tInstrument1\tInstrument2\tPlumeAltitude_masl\tPlumeHeightError_m\tWindDirection_deg\tWindDirectionError_deg\tPlumeCentre1_deg\tPlumeCentreError1_deg\tPlumeCentre2_deg\tPlumeCentreError2_deg\n");
	}

	POSITION pos = geometryResults.GetHeadPosition();
	while(pos != NULL){
		Geometry::CGeometryResult *result = geometryResults.GetNext(pos);
		// write the file
		if(result->m_calculationType == Meteorology::MET_GEOMETRY_CALCULATION){
			fprintf(f, "%04d.%02d.%02d\t",	result->m_averageStartTime.year, result->m_averageStartTime.month, result->m_averageStartTime.day);
			fprintf(f, "%02d:%02d:%02d\t",	result->m_averageStartTime.hour, result->m_averageStartTime.minute, result->m_averageStartTime.second);
			fprintf(f, "%.1lf\t",			result->m_startTimeDifference / 60.0);
			fprintf(f, "%s\t%s\t",			result->m_instr1, result->m_instr2);
			fprintf(f, "%.0lf\t%.0lf\t",	result->m_plumeAltitude, result->m_plumeAltitudeError);
			fprintf(f, "%.0lf\t%.0lf\t",	result->m_windDirection, result->m_windDirectionError);

			fprintf(f, "%.1f\t%.1f\t",		result->m_plumeCentre1, result->m_plumeCentreError1);
			fprintf(f, "%.1f\t%.1f\n",		result->m_plumeCentre2, result->m_plumeCentreError2);
		}else{
			fprintf(f, "%04d.%02d.%02d\t",	result->m_averageStartTime.year, result->m_averageStartTime.month, result->m_averageStartTime.day);
			fprintf(f, "%02d:%02d:%02d\t",	result->m_averageStartTime.hour, result->m_averageStartTime.minute, result->m_averageStartTime.second);
			fprintf(f, "0\t");
			fprintf(f, "%s\t\t",			result->m_instr1);
			fprintf(f, "%.0lf\t%.0lf\t",	result->m_plumeAltitude, result->m_plumeAltitudeError);
			fprintf(f, "%.0lf\t%.0lf\t",	result->m_windDirection, result->m_windDirectionError);

			fprintf(f, "%.1f\t%.1f\t",		result->m_plumeCentre1, result->m_plumeCentreError1);
			fprintf(f, "0\t0\n");
		}
	}
	fclose(f);
}

// 5. Insert the calculated geometries into the plume height database
void CPostProcessing::InsertCalculatedGeometriesIntoDataBase(const CList <Geometry::CGeometryResult*, Geometry::CGeometryResult*> &geometryResults){
	Meteorology::CWindField windField;
	CDateTime validFrom, validTo;
	Configuration::CInstrumentLocation location;

	POSITION pos = geometryResults.GetHeadPosition();
	while(pos != NULL){
		Geometry::CGeometryResult *result = geometryResults.GetNext(pos);

		if(result->m_plumeAltitude > 0.0){
			// insert the plume height into the plume height database
			this->m_plumeDataBase.InsertPlumeHeight(*result);
		}

		if(result->m_windDirection > NOT_A_NUMBER){
			// get the location of the instrument at the time of the measurement
			g_setup.GetInstrumentLocation(result->m_instr1, result->m_averageStartTime, location);
		
			// get the time-interval that the measurement is valid for
			validFrom = CDateTime(result->m_averageStartTime);
			validFrom.Decrement(g_userSettings.m_calcGeometryValidTime);
			validTo   = CDateTime(result->m_averageStartTime);
			validTo.Increment(g_userSettings.m_calcGeometryValidTime);
		
			// insert the wind-direction into the wind database
			m_windDataBase.InsertWindDirection(validFrom, validTo, result->m_windDirection, result->m_windDirectionError, result->m_calculationType, NULL);
		}
	}
}

/** This calculates the wind speeds from the dual-beam measurements that has been made
	@param evalLogs - list of CStrings, each holding the full path and filename 
		of an evaluation-log file. Only the measurements containing a dual-beam measurement
		will be considered.
	The plume heights are taken from the database 'm_plumeDataBase' and the results are written
	to the database 'm_windDataBase'
*/
void CPostProcessing::CalculateDualBeamWindSpeeds(const CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> &evalLogs){
	CList <CString, CString &> masterList; // list of wind-measurements from the master channel
	CList <CString, CString &> slaveList;  // list of wind-measurements from the slave channel
	CList <CString, CString &> heidelbergList;  // list of wind-measurements from the Heidelbergensis
	CDateTime validFrom, validTo;
	
	CString serial, serial2, fileName, fileName2, nonsenseString;
	CString userMessage, windLogFile;
	CDateTime startTime, startTime2;
	int channel, channel2, nWindMeasFound = 0;
	MEASUREMENT_MODE meas_mode, meas_mode2;
	Configuration::CInstrumentLocation location;
	WindSpeedMeasurement::CWindSpeedCalculator calculator;
	Geometry::CPlumeHeight plumeHeight;
	Meteorology::CWindField windField, oldWindField;

	// -------------------------------- step 1. -------------------------------------
	// search through 'evalLogs' for dual-beam measurements from master and from slave
	POSITION pos = evalLogs.GetHeadPosition();
	while(pos != NULL){
		const CString &fileNameAndPath = evalLogs.GetNext(pos).m_evalLogFile[g_userSettings.m_mainFitWindow];
		
		// to know the start-time of the measurement, we need to 
		//	extract just the file-name, i.e. remove the path
		fileName = CString(fileNameAndPath);
		Common::GetFileName(fileName);
		
		FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(fileName, startTime, serial, channel, meas_mode);
		
		if(meas_mode == MODE_WINDSPEED){
			++nWindMeasFound;
			// first check if this is a heidelberg instrument
			if(g_setup.GetInstrumentLocation(serial, startTime, location))
				continue;
			
			if(location.m_instrumentType == INSTR_HEIDELBERG){
				// this is a heidelberg instrument
				heidelbergList.AddTail(CString(fileNameAndPath));
			}else{
				// this is a gothenburg instrument
				if(channel == 0){
					masterList.AddTail(CString(fileNameAndPath));
				}else if(channel == 1){
					slaveList.AddTail(CString(fileNameAndPath));
				}
			}
		}
	}
	if(nWindMeasFound == 0){
		ShowMessage("No dual-beam wind speed measurements found.");
		return; // if nothing was found...
	}

	userMessage.Format("%d dual-beam wind speed measurements found. Calculating wind-speeds", nWindMeasFound);
	ShowMessage(userMessage);
	
	// Create the dual-beam log-file
	windLogFile.Format("%s\\DualBeamLog.txt", g_userSettings.m_outputDirectory);
	calculator.WriteWindSpeedLogHeader(windLogFile);
	

	// -------------------------------- step 2. -------------------------------------
	// loop through each of the measurements from the heidelberg instruments
	//	and calculate the wind speed for each measurement
	pos = heidelbergList.GetHeadPosition();
	while(pos != NULL){
		const CString &fileNameAndPath = heidelbergList.GetNext(pos);
		
		// to know the start-time of the measurement, we need to 
		//	extract just the file-name, i.e. remove the path
		fileName = CString(fileNameAndPath);
		Common::GetFileName(fileName);
		FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(fileName, startTime, serial, channel, meas_mode);
		
		// Get the plume height at the time of the measurement
		m_plumeDataBase.GetPlumeHeight(startTime, plumeHeight);
		
		// Get the location of the instrument at the time of the measurement
		g_setup.GetInstrumentLocation(serial, startTime, location);
		
		// calculate the speed of the wind at the time of the measurement
		if(0 == calculator.CalculateWindSpeed(fileNameAndPath, nonsenseString, location, plumeHeight, windField)){
			// append the results to file
			calculator.AppendResultToFile(windLogFile, startTime, location, plumeHeight, windField);

			// insert the newly calculated wind-speed into the database
			if(windField.GetWindSpeedError() > g_userSettings.m_dualBeam_MaxWindSpeedError){
				userMessage.Format("-Calculated a wind-speed of %.1lf � %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Error too large, measurement discarded.", windField.GetWindSpeed(), windField.GetWindSpeedError(), 
					startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);
			}else{
				// tell the user...
				userMessage.Format("+Calculated a wind-speed of %.1lf � %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Measurement accepted", windField.GetWindSpeed(), windField.GetWindSpeedError(), 
					startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);

				// get the time-interval that the measurement is valid for
				windField.GetValidTimeFrame(validFrom, validTo);

				// insert the new wind speed into the database
				m_windDataBase.InsertWindSpeed(validFrom, validTo, windField.GetWindSpeed(), windField.GetWindSpeedError(), Meteorology::MET_DUAL_BEAM_MEASUREMENT, NULL);
			}
			ShowMessage(userMessage);
		}else{
			userMessage.Format("Failed to calculate wind speed from measurement: %s", fileName);
			ShowMessage(userMessage);					
		}
	}
	
	// -------------------------------- step 3. -------------------------------------
	// loop through each of the measurements from a master-channel and try to match them with a measurement
	//	from a slave channel...
	pos = masterList.GetHeadPosition();
	while(pos != NULL){
		const CString &fileNameAndPath = masterList.GetNext(pos);
		
		// extract just the file-name, i.e. remove the path
		fileName = CString(fileNameAndPath);
		Common::GetFileName(fileName);
		
		FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(fileName, startTime, serial, channel, meas_mode);
		
		// now check if we can match this one with a file in the slave-channel
		POSITION pos2 = slaveList.GetHeadPosition();
		while(pos2 != NULL){
			const CString &fileNameAndPath2 = slaveList.GetNext(pos2);
			
			// extract just the file-name, i.e. remove the path
			fileName2 = CString(fileNameAndPath2);
			Common::GetFileName(fileName2);
			
			FileHandler::CEvaluationLogFileHandler::GetInfoFromFileName(fileName2, startTime2, serial2, channel2, meas_mode2);
			
			if(Equals(serial, serial2) && (startTime == startTime2)){
				// we have found a match!!!

				// Get the plume height at the time of the measurement
				m_plumeDataBase.GetPlumeHeight(startTime, plumeHeight);
				
				// Get the location of the instrument at the time of the measurement
				g_setup.GetInstrumentLocation(serial, startTime, location);
				
				// calculate the speed of the wind at the time of the measurement
				if(0 == calculator.CalculateWindSpeed(fileNameAndPath, fileNameAndPath2, location, plumeHeight, windField)){
					// append the results to file
					calculator.AppendResultToFile(windLogFile, startTime, location, plumeHeight, windField);

					// insert the newly calculated wind-speed into the database
					if(windField.GetWindSpeedError() > g_userSettings.m_dualBeam_MaxWindSpeedError){
						userMessage.Format("-Calculated a wind-speed of %.1lf � %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Error too large, measurement discarded.", windField.GetWindSpeed(), windField.GetWindSpeedError(), 
							startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);
					}else{
						// tell the user...
						userMessage.Format("+Calculated a wind-speed of %.1lf � %.1lf m/s on %04d.%02d.%02d at %02d:%02d. Measurement accepted", windField.GetWindSpeed(), windField.GetWindSpeedError(), 
							startTime.year, startTime.month, startTime.day, startTime.hour, startTime.minute);

						windField.GetValidTimeFrame(validFrom, validTo);

						// insert the new wind speed into the database
						m_windDataBase.InsertWindSpeed(validFrom, validTo, windField.GetWindSpeed(), windField.GetWindSpeedError(), Meteorology::MET_DUAL_BEAM_MEASUREMENT, NULL);
					}
					ShowMessage(userMessage);

				}else{
					userMessage.Format("Failed to calculate wind speed from measurement: %s", fileName);
					ShowMessage(userMessage);					
				}
			}
		}
	}
}

/** Sorts the evaluation logs in order of increasing time */
void CPostProcessing::SortEvaluationLogs(CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> &evalLogs){
	CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> left;
	CList <Evaluation::CExtendedScanResult, Evaluation::CExtendedScanResult &> right;
	POSITION pos1, pos2;

	// If this list consists of only one element, then we're done
	if(evalLogs.GetCount() <= 1)
		return;

	// Divide the list into two, and sort each one of them
	pos1 = evalLogs.GetHeadPosition();
	int index = 0;
	while(pos1 != NULL){
		Evaluation::CExtendedScanResult &log = evalLogs.GetNext(pos1);
		if(index % 2 == 0)
			left.AddTail(log);
		else
			right.AddTail(log);
		++index;		
	}

	SortEvaluationLogs(left);
	SortEvaluationLogs(right);

	// Merge the two lists into one, sorted list
	evalLogs.RemoveAll();
	pos1 = left.GetHeadPosition();
	pos2 = right.GetHeadPosition();
	while(pos1 != NULL && pos2 != NULL){
		Evaluation::CExtendedScanResult &log1 = left.GetAt(pos1);
		Evaluation::CExtendedScanResult &log2 = right.GetAt(pos2);
		
		if(log2.m_startTime < log1.m_startTime){
			evalLogs.AddTail(log2);
			right.GetNext(pos2);
		}else{
			evalLogs.AddTail(log1);
			left.GetNext(pos1);
		}
	}
	while(pos1 != NULL){
		evalLogs.AddTail(left.GetNext(pos1));
	}
	while(pos2 != NULL){
		evalLogs.AddTail(right.GetNext(pos2));
	}
	

	return;
}


/** Takes care of uploading the result files to the FTP server */
void CPostProcessing::UploadResultsToFTP(){
	Communication::CFTPServerConnection *connection = new Communication::CFTPServerConnection();
	CString fileName;

	// Generate a list with all the files we want to upload.
	CList <CString, CString&> fileList;
	
	// 1. the geometry log file
	fileName.Format("%s\\GeometryLog.txt", g_userSettings.m_outputDirectory);
	fileList.AddTail(fileName);

	// 2. the generated wind field
	fileName.Format("%s\\GeneratedWindField.wxml", g_userSettings.m_outputDirectory);
	fileList.AddTail(fileName);

	// 3. the generated flux log
	fileName.Format("%s\\FluxLog.txt", g_userSettings.m_outputDirectory);
	fileList.AddTail(fileName);

	// upload the files
	connection->UploadResults("129.16.35.206", "novacUser", "iht-1inks.", fileList);
	
	delete connection;
	return;
}