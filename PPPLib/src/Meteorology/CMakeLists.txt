# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPPLIB_METEOROLOGY_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/Meteorology/MeteorologySource.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Meteorology/WindDataBase.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Meteorology/WindField.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Meteorology/XMLWindFileReader.h
    PARENT_SCOPE)
    
    
set(NPPLIB_METEOROLOGY_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/WindDataBase.cpp
    ${CMAKE_CURRENT_LIST_DIR}/WindField.cpp
    ${CMAKE_CURRENT_LIST_DIR}/XMLWindFileReader.cpp
    PARENT_SCOPE)

    