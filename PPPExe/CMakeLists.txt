##########################################################################
## This is the main project of Novac Post Processing Program (NovacPPP) ##
##########################################################################

cmake_minimum_required (VERSION 3.6)

# Define the different groups of files
add_subdirectory(Common)
add_subdirectory(Stratosphere)

## --------- Creating NovacPPP ---------

# Add the different components
add_executable(NovacPPP
    ${CMAKE_CURRENT_LIST_DIR}/NovacPPP.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ObservatoryInfo.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ObservatoryInfo.h
    ${CMAKE_CURRENT_LIST_DIR}/PostProcessing.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PostProcessing.h

    ${NPP_COMMON_HEADERS} 
    ${NPP_COMMON_SOURCES}
    ${NPP_COMMON_SPECTRA_HEADERS}
    ${NPP_COMMON_SPECTRA_SOURCES}
    ${NPP_COMMON_SPECTRUM_FORMAT_HEADERS}
    ${NPP_COMMON_SPECTRUM_FORMAT_SOURCES}
    ${NPP_STRATOSPHERE_HEADERS}
    ${NPP_STRATOSPHERE_SOURCES}
)

target_include_directories(NovacPPP PRIVATE 
    ${PppLib_INCLUDE_DIRS}
    ${SPECTRALEVAUATION_INCLUDE_DIRS}
    ${PocoNet_DIR}
    ${PocoFoundation_DIR})

target_link_libraries(NovacPPP PRIVATE PPPLib ${Poco_LIBRARIES})

IF(WIN32)
    # TODO: Enable /w44365 (conversion signed/unsigned) as well here
    target_compile_options(NovacPPP PUBLIC /W4 /WX /sdl)
    target_compile_definitions(NovacPPP PUBLIC _CRT_SECURE_NO_WARNINGS)
ELSE()
    target_compile_options(NovacPPP PUBLIC -Wall -std=c++14)
ENDIF()


IF(WIN32)
    #Organize the files into folders, for ease of use (Windows and Visual Studio only)
    source_group("Header Files\\Common"                 FILES ${NPP_COMMON_HEADERS})
    source_group("Header Files\\Configuration"          FILES ${NPP_CONFIGURATION_HEADERS})
    source_group("Header Files\\Stratosphere"           FILES ${NPP_STRATOSPHERE_HEADERS})

    source_group("Source Files\\Common"                 FILES ${NPP_COMMON_SOURCES})
    source_group("Source Files\\Configuration"          FILES ${NPP_CONFIGURATION_SOURCES})
    source_group("Source Files\\Stratosphere"           FILES ${NPP_STRATOSPHERE_SOURCES})
endif()

