# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPPLIB_FILE_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/File/EvaluationConfigurationParser.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/File/EvaluationLogFileHandler.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/File/Filesystem.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/File/ProcessingFileReader.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/File/SetupFileReader.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/File/XMLFileReader.h
    PARENT_SCOPE)    
    
set(NPPLIB_FILE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationConfigurationParser.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationLogFileHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Filesystem.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ProcessingFileReader.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SetupFileReader.cpp
    ${CMAKE_CURRENT_LIST_DIR}/XMLFileReader.cpp
    PARENT_SCOPE)

    