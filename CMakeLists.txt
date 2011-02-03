########################################################
# cmake file for building PandoraMonitoring
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( PandoraMonitoring )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 1 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( ROOT 5.26.00 REQUIRED ) # minimum required ROOT version
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} )


SET( LIB_EVE LIB_EVE-NOTFOUND )
FIND_FILE( LIB_EVE lib/libEve.so ${ROOT_DIR} )
IF( LIB_EVE )
    MESSAGE( STATUS "Check for ROOT_EVE_LIBRARY: " ${LIB_EVE} " -- ok " )
    LINK_LIBRARIES( ${LIB_EVE} )
    ADD_DEFINITIONS( -DROOT_EVE )
ELSE()
    MESSAGE( STATUS "ROOT_EVE_LIBRARY not found" )
ENDIF()


#FIND_PACKAGE( PandoraPFANew QUIET )

IF( PandoraPFANew_FOUND )
    INCLUDE_DIRECTORIES( ${PandoraPFANew_INCLUDE_DIRS} )
ELSEIF( PandoraPFANew_DIR )
    INCLUDE_DIRECTORIES( "${PandoraPFANew_DIR}/Framework/include" )
ELSE()
    MESSAGE( SEND_ERROR "Could not find PandoraPFANew required to build PandoraMonitoring. Please set PandoraPFANew_DIR" )
ENDIF()



### DOCUMENTATION ###########################################################
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC )

    FIND_PACKAGE( Doxygen )

    IF( DOXYGEN_EXECUTABLE )
        ADD_SUBDIRECTORY( ./doc )
    ELSE()
        MESSAGE( SEND_ERROR "Could not find doxygen required to build documentation" )
        MESSAGE( "Please install doxygen or set INSTALL_DOC to OFF" )
        MESSAGE( "" )
    ENDIF()

ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# Add additional linker flag for Apple to defer linking against ROOT until PandoraPFANew build.
IF( APPLE ) 
    SET( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup" )
    # Add Unix compiler definition.  Already happens automatically on Linux.
    IF( UNIX )
        ADD_DEFINITIONS(-D__unix__)
    ENDIF()
ENDIF()


# include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( ${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake ${PROJECT_NAME}LibDeps.cmake )
