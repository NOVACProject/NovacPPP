# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPPLIB_GEOMETRY_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/Geometry/GeometryCalculator.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Geometry/GeometryResult.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Geometry/PlumeDataBase.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Geometry/PlumeHeight.h
    PARENT_SCOPE)
    
    
set(NPPLIB_GEOMETRY_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/GeometryCalculator.cpp
    ${CMAKE_CURRENT_LIST_DIR}/GeometryResult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PlumeDataBase.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PlumeHeight.cpp
    PARENT_SCOPE)

    