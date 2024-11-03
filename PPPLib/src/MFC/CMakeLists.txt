# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPPLIB_MFC_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/MFC/CCriticalSection.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/MFC/CFileUtils.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/MFC/CFtpUtils.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/MFC/CList.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/MFC/CSingleLock.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/MFC/CStdioFile.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/MFC/CString.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/MFC/CStringTokenizer.h
    PARENT_SCOPE)
    
    
set(NPPLIB_MFC_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/CFileUtils.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/CFtpUtils.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/CStdioFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/CString.cpp
    ${CMAKE_CURRENT_LIST_DIR}/CStringTokenizer.cpp
    PARENT_SCOPE)

    