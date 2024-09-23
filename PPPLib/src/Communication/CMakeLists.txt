# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPPLIB_COMMUNICATION_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/Communication/FTPServerConnection.h
    PARENT_SCOPE)
    
    
set(NPPLIB_COMMUNICATION_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/FTPServerConnection.cpp
    PARENT_SCOPE)
