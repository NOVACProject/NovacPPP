# This is the main project of Novac Post Processing Program (NovacPPP)

cmake_minimum_required (VERSION 3.6)

set(PppTests_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/include)

# Add the different components
add_executable(PPPTests
    ${CMAKE_CURRENT_LIST_DIR}/include/StdOutLogger.h

    ${CMAKE_CURRENT_LIST_DIR}/src/main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_CommandLineParser.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_CArray.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_CFileUtils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_CFtpUtils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_CList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_CString.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_CStringTokenizer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_EvaluationConfiguration.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_EvaluationConfigurationParser.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_NovacPPPConfiguration.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_PostCalibrationStatistics.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_ProcessingFileReader.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_SetupFileReader.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/UnitTest_XmlWindFileReader.cpp
)

target_link_libraries(PPPTests PRIVATE PPPLib)
    
target_include_directories(PPPTests PRIVATE ${PppTests_INCLUDE_DIRS} ${PppLib_INCLUDE_DIRS})

IF(WIN32)
    target_compile_options(PPPTests PRIVATE /W4 /WX /sdl)
ELSE()
    target_compile_options(PPPTests PRIVATE -Wall -std=c++11)
ENDIF()