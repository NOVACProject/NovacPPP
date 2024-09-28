# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)


set(NPPLIB_WINDMEASUREMENT_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/WindMeasurement/WindSpeedCalculator.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/WindMeasurement/WindSpeedMeasSettings.h
    PARENT_SCOPE)
    
    
set(NPPLIB_WINDMEASUREMENT_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/WindSpeedCalculator.cpp
    PARENT_SCOPE)

    