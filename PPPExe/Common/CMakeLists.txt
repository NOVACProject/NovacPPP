# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPP_COMMON_HEADERS
    ${CMAKE_CURRENT_LIST_DIR}/Common.h
    ${CMAKE_CURRENT_LIST_DIR}/Version.h
    PARENT_SCOPE)
    
set(NPP_COMMON_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/Common.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Version.cpp
    PARENT_SCOPE)


