# This is the main project of Novac Post Processing Program (NovacPPP)

cmake_minimum_required (VERSION 3.6)

set(PppLib_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/include CACHE STRING "PPPLib headers")

# Define the different groups of files
add_subdirectory(src/Calibration)
add_subdirectory(src/Communication)
add_subdirectory(src/Configuration)
add_subdirectory(src/File)
add_subdirectory(src/Flux)
add_subdirectory(src/Geometry)
add_subdirectory(src/Evaluation)
add_subdirectory(src/Meteorology)
add_subdirectory(src/MFC)
add_subdirectory(src/WindMeasurement)

# Add the different components
add_library(PPPLib
    ${NPPLIB_MFC_HEADERS} 
    ${NPPLIB_MFC_SOURCES}
    ${NPPLIB_CALIBRATION_HEADERS}
    ${NPPLIB_CALIBRATION_SOURCES}
    ${NPPLIB_COMMUNICATION_HEADERS}
    ${NPPLIB_COMMUNICATION_SOURCES}
    ${NPPLIB_EVALUATION_HEADERS}
    ${NPPLIB_EVALUATION_SOURCES}
    ${NPPLIB_FLUX_HEADERS}
    ${NPPLIB_FLUX_SOURCES}
    ${NPPLIB_GEOMETRY_HEADERS}
    ${NPPLIB_GEOMETRY_SOURCES}
    ${NPPLIB_CONFIGURATION_HEADERS}
    ${NPPLIB_CONFIGURATION_SOURCES}
    ${NPPLIB_FILE_HEADERS}
    ${NPPLIB_FILE_SOURCES}
    ${NPPLIB_METEOROLOGY_HEADERS}
    ${NPPLIB_METEOROLOGY_SOURCES}
    ${NPPLIB_WINDMEASUREMENT_HEADERS}
    ${NPPLIB_WINDMEASUREMENT_SOURCES}

    ${PppLib_INCLUDE_DIRS}/PPPLib/ContinuationOfProcessing.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Definitions.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Logging.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/PPPLib.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/PostProcessingStatistics.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/PostProcessingUtils.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/SpectrometerId.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/VolcanoInfo.h

    ${CMAKE_CURRENT_LIST_DIR}/src/ContinuationOfProcessing.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/VolcanoInfo.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/PostProcessingStatistics.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/PostProcessingUtils.cpp
)

target_include_directories(PPPLib PRIVATE 
    ${PppLib_INCLUDE_DIRS}
    ${SPECTRALEVAUATION_INCLUDE_DIRS}
    ${PocoNet_DIR}
    ${PocoFoundation_DIR})

# Some classes are using Poco, though we should be able to remove these eventually
target_link_libraries(PPPLib PUBLIC ${Poco_LIBRARIES} NovacSpectralEvaluation)


IF(WIN32)
    target_compile_options(PPPLib PRIVATE /W4 /WX /w44365 /sdl)
    target_compile_definitions(PPPLib PRIVATE _CRT_SECURE_NO_WARNINGS)
ELSE()
    target_compile_options(PPPLib PRIVATE -Wall -std=c++14)
ENDIF()

IF(WIN32)
    #Organize the files into folders, for ease of use (Windows and Visual Studio only)
    source_group("Header Files\\MFC"                FILES ${NPPLIB_MFC_HEADERS})
    source_group("Header Files\\Communication"      FILES ${NPPLIB_COMMUNICATION_HEADERS})
    source_group("Header Files\\Configuration"      FILES ${NPPLIB_CONFIGURATION_HEADERS})
    source_group("Header Files\\Calibration"        FILES ${NPPLIB_CALIBRATION_HEADERS})
    source_group("Header Files\\Evaluation"         FILES ${NPPLIB_EVALUATION_HEADERS})
    source_group("Header Files\\File"               FILES ${NPPLIB_FILE_HEADERS})
    source_group("Header Files\\Flux"               FILES ${NPPLIB_FLUX_HEADERS})
    source_group("Header Files\\Geometry"           FILES ${NPPLIB_GEOMETRY_HEADERS})
    source_group("Header Files\\Meteorology"        FILES ${NPPLIB_METEOROLOGY_HEADERS})
    source_group("Header Files\\WindMeasurement"    FILES ${NPPLIB_WINDMEASUREMENT_HEADERS})

    source_group("Source Files\\MFC"                FILES ${NPPLIB_MFC_SOURCES})
    source_group("Source Files\\Communication"      FILES ${NPPLIB_COMMUNICATION_SOURCES})
    source_group("Source Files\\Configuration"      FILES ${NPPLIB_CONFIGURATION_SOURCES})
    source_group("Source Files\\Calibration"        FILES ${NPPLIB_CALIBRATION_SOURCES})
    source_group("Source Files\\Evaluation"         FILES ${NPPLIB_EVALUATION_SOURCES})
    source_group("Source Files\\File"               FILES ${NPPLIB_FILE_SOURCES})
    source_group("Source Files\\Flux"               FILES ${NPPLIB_FLUX_SOURCES})
    source_group("Source Files\\Geometry"           FILES ${NPPLIB_GEOMETRY_SOURCES})
    source_group("Source Files\\Meteorology"        FILES ${NPPLIB_METEOROLOGY_SOURCES})
    source_group("Source Files\\WindMeasurement"    FILES ${NPPLIB_WINDMEASUREMENT_SOURCES})
endif()

