# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPPLIB_CONFIGURATION_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/CommandLineParser.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/DarkCorrectionConfiguration.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/EvaluationConfiguration.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/InstrumentCalibrationConfiguration.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/InstrumentConfiguration.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/InstrumentType.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/InstrumentLocation.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/LocationConfiguration.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/NovacPPPConfiguration.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/ProcessingMode.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Configuration/UserConfiguration.h
    PARENT_SCOPE)    
    
set(NPPLIB_CONFIGURATION_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/CommandLineParser.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/DarkCorrectionConfiguration.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationConfiguration.cpp
    ${CMAKE_CURRENT_LIST_DIR}/InstrumentConfiguration.cpp
    ${CMAKE_CURRENT_LIST_DIR}/InstrumentLocation.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/LocationConfiguration.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/NovacPPPConfiguration.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/UserConfiguration.cpp 
    PARENT_SCOPE)

    