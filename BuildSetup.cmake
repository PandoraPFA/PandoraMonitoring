#############################################################################
# cmake build setup for PandoraMonitoring
#
# For building PandoraMonitoring with cmake type:
# (1) $ mkdir build
# (2) $ cd build
# (3) $ cmake -C ../BuildSetup.cmake ..
# (4) $ make install
#
# @author Jan Engels, DESY
#############################################################################


#############################################################################
# Setup path variables
#############################################################################

# ILC_HOME
SET( ILC_HOME "path_to_ilc_home"
    CACHE PATH "Path to ILC Software" FORCE )

# Path to PandoraPFANew
SET( PandoraPFANew_HOME "${ILC_HOME}/PandoraPFANew"
    CACHE PATH "Path to PandoraPFANew" FORCE )

# Path to LCIO
SET( LCIO_HOME "${ILC_HOME}/lcio/v01-11"
    CACHE PATH "Path to LCIO" FORCE )

# Path to GEAR
SET( GEAR_HOME "${ILC_HOME}/gear/v00-11-01"
    CACHE PATH "Path to GEAR" FORCE )

# Path to ROOT
SET( ROOT_HOME "${ILC_HOME}/root/5.16.00"
     CACHE PATH "Path to ROOT" FORCE )

# CMake Modules Path
SET( CMAKE_MODULE_PATH "${ILC_HOME}/CMakeModules/v01-08"
    CACHE PATH "Path to CMake Modules" FORCE )

###############################################
# Project Options
###############################################

#SET( INSTALL_DOC OFF CACHE BOOL "Set to OFF to skip build/install Documentation" FORCE )

# set cmake build type
# possible options are: None Debug Release RelWithDebInfo MinSizeRel
#SET( CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build" FORCE )

###############################################
# Advanced Options
###############################################

#SET( BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF to build static libraries" FORCE )

# installation path for PandoraPFA
#SET( CMAKE_INSTALL_PREFIX "foo/bar" CACHE STRING "Where to install PandoraMonitoring" FORCE )
