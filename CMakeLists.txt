########################################################
# cmake file for building PandoraMonitoring
########################################################

# PandoraMonitoring requires ROOT
FIND_PACKAGE( ROOT 5.26.00 REQUIRED COMPONENTS Eve Geom RGL EG ) # minimum required ROOT version
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
LINK_LIBRARIES( ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES} )


### LIBRARY #################################################################

INCLUDE_DIRECTORIES( ./include )
INCLUDE_DIRECTORIES( ../Framework/include )

AUX_SOURCE_DIRECTORY( ./src/ Monitoring_SRCS )

# require proper C++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

# need long long for int64 for now
ADD_DEFINITIONS( "-Wno-long-long -Wno-sign-compare -fno-strict-aliasing" )


# add library
ADD_SHARED_LIBRARY( PandoraMonitoring ${Monitoring_SRCS} )
TARGET_LINK_LIBRARIES( PandoraMonitoring PandoraFramework )

# install library
INSTALL_SHARED_LIBRARY( PandoraMonitoring DESTINATION lib )

# install header files
INSTALL_DIRECTORY( ./include DESTINATION . )

