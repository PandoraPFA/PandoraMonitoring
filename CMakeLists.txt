########################################################
# cmake file for building PandoraMonitoring
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2 FATAL_ERROR)
########################################################


# project name
PROJECT( PandoraMonitoring )


# project version
SET( PandoraMonitoring_VERSION_MAJOR 0 )
SET( PandoraMonitoring_VERSION_MINOR 17 )
SET( PandoraMonitoring_VERSION_PATCH 0 )


### DEPENDENCIES ############################################################
FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES QUIET )

IF( ILCUTIL_FOUND )
    INCLUDE( ilcsoft_default_settings )
ELSE()
    INCLUDE( Default_settings )
ENDIF()

FIND_PACKAGE( PandoraSDK 0.17 REQUIRED )
FOREACH( pkg PandoraSDK )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

FIND_PACKAGE( ROOT 5.26.00 REQUIRED COMPONENTS Eve Geom RGL EG ) # minimum required ROOT version
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )


### DOCUMENTATION ############################################################
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()


### SOURCE ##################################################################

#include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )
ADD_DEFINITIONS( "-DMONITORING" )

AUX_SOURCE_DIRECTORY( ./src/ PandoraMonitoring_SRCS )

# require proper C++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

# need long long for int64 for now
ADD_DEFINITIONS( "-Wno-long-long -Wno-sign-compare -fno-strict-aliasing" )


### LIBRARY #################################################################

# add library
ADD_SHARED_LIBRARY( PandoraMonitoring ${PandoraMonitoring_SRCS} )

# install library
INSTALL_SHARED_LIBRARY( PandoraMonitoring DESTINATION lib )

# install header files
INSTALL_DIRECTORY( ./include DESTINATION . )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install cmake configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( PandoraMonitoringConfig.cmake PandoraMonitoringConfigVersion.cmake PandoraMonitoringLibDeps.cmake )
