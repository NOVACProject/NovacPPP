# This is the main project of Novac Post Processing Program (NovacPPP)

cmake_minimum_required (VERSION 3.6)

set(PppLib_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/include CACHE STRING "PPPLib headers")

# Define the different groups of files
add_subdirectory(src/Calibration)
add_subdirectory(src/Communication)
add_subdirectory(src/Configuration)
add_subdirectory(src/File)
add_subdirectory(src/Meteorology)
add_subdirectory(src/MFC)

# Add the different components
add_library(PPPLib
    ${NPPLIB_MFC_HEADERS} 
    ${NPPLIB_MFC_SOURCES}
    ${NPPLIB_CALIBRATION_HEADERS}
    ${NPPLIB_CALIBRATION_SOURCES}
    ${NPP_COMMUNICATION_HEADERS}
    ${NPP_COMMUNICATION_SOURCES}
    ${NPPLIB_CONFIGURATION_HEADERS}
    ${NPPLIB_CONFIGURATION_SOURCES}
    ${NPPLIB_FILE_HEADERS}
    ${NPPLIB_FILE_SOURCES}
    ${NPP_METEOROLOGY_HEADERS}
    ${NPP_METEOROLOGY_SOURCES}

    ${PppLib_INCLUDE_DIRS}/PPPLib/Definitions.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Gps.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Logging.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Measurement.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Molecule.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/PPPLib.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/PostProcessingStatistics.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/SpectrometerId.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/ThreadUtils.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/VolcanoInfo.h

    ${CMAKE_CURRENT_LIST_DIR}/src/Gps.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/VolcanoInfo.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/Molecule.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/PostProcessingStatistics.cpp
)

target_include_directories(PPPLib PRIVATE 
    ${PppLib_INCLUDE_DIRS}
    ${SPECTRALEVAUATION_INCLUDE_DIRS}
    ${PocoNet_DIR}
    ${PocoFoundation_DIR})

# Some classes are using Poco, though we should be able to remove these eventually
target_link_libraries(PPPLib PUBLIC ${Poco_LIBRARIES} NovacSpectralEvaluation)


IF(WIN32)
    target_compile_options(PPPLib PRIVATE /W4 /WX /sdl)
    target_compile_definitions(PPPLib PRIVATE _CRT_SECURE_NO_WARNINGS)
ELSE()
    target_compile_options(PPPLib PRIVATE -Wall -std=c++14)
ENDIF()

IF(WIN32)
    #Organize the files into folders, for ease of use (Windows and Visual Studio only)
    source_group("Header Files\\MFC"            FILES ${NPPLIB_MFC_HEADERS})
    source_group("Header Files\\Communication"  FILES ${NPP_COMMUNICATION_HEADERS})
    source_group("Header Files\\Configuration"  FILES ${NPPLIB_CONFIGURATION_HEADERS})
    source_group("Header Files\\Calibration"    FILES ${NPPLIB_CALIBRATION_HEADERS})
    source_group("Header Files\\File"           FILES ${NPPLIB_FILE_HEADERS})
    source_group("Header Files\\Meteorology"    FILES ${NPP_METEOROLOGY_HEADERS})

    source_group("Source Files\\MFC"            FILES ${NPPLIB_MFC_SOURCES})
    source_group("Source Files\\Communication"  FILES ${NPP_COMMUNICATION_SOURCES})
    source_group("Source Files\\Configuration"  FILES ${NPPLIB_CONFIGURATION_SOURCES})
    source_group("Source Files\\Calibration"    FILES ${NPPLIB_CALIBRATION_SOURCES})
    source_group("Source Files\\File"           FILES ${NPPLIB_FILE_SOURCES})
    source_group("Source Files\\Meteorology"    FILES ${NPP_METEOROLOGY_SOURCES})
endif()

