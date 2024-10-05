# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(NPPLIB_EVALUATION_HEADERS
    ${PppLib_INCLUDE_DIRS}/PPPLib/Evaluation/EvaluationUtils.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Evaluation/ExtendedScanResult.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Evaluation/PostEvaluationController.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Evaluation/PostEvaluationIO.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Evaluation/ScanEvaluation.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Evaluation/ScanResult.h
    PARENT_SCOPE)
    
    
set(NPPLIB_EVALUATION_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationUtils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PostEvaluationController.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PostEvaluationIO.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScanEvaluation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScanResult.cpp
    PARENT_SCOPE)
