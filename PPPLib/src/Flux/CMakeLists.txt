# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPPLIB_FLUX_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/Flux/FluxCalculator.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Flux/FluxResult.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Flux/FluxStatistics.h
    PARENT_SCOPE)
    
    
set(NPPLIB_FLUX_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/FluxCalculator.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FluxResult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FluxStatistics.cpp
    PARENT_SCOPE)

    