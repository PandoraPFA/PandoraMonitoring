README
======

Summary:
--------
Download, compilation and installation instructions of PandoraMonitoring for PandoraPFANew.


Installation of PandoraMonitoring
==============================================================

Requirements:
-------------
For the default installation of PandoraMonitoring using cmake the following software is needed:
- ROOT (>=5.26, built with ROOT TEve included)
- LCIO 
- CMakeModules (if cmake is used for compilation)
- PandoraPFANew

Getting PandoraMonitoring:
--------------------------
- checkout PandoraMonitoring-source from:
svn co https://svnsrv.desy.de/public/PandoraPFANew/PandoraMonitoring/trunk PandoraMonitoring

Installing PandoraMonitoring:
--------------------------
cd PandoraMonitoring
mkdir build
cd build
cmake -C ILCSoft.cmake ..
make install

(where ILCSoft.cmake contains locations of required packages PandoraPFANew, ROOT and CMakeModules)


Installation of PandoraMonitoring without CMake
==============================================================

If cmake is not available, please edit the supplied "Makefile" to specify details of the project directory 
and paths to the project dependencies. The package is then built with the command "make". A commandline 
switch "BUILD_32BIT_COMPATIBLE=1" can be used.


Recompilation of PandoraPFANew with PandoraMonitoring support:
==============================================================

- tell PandoraPFANew the path to the monitoring directory
cmake -C ILCSoft.cmake -DPandoraMonitoring_DIR=/<path_to_monitoring_directory> ..
make install


Usage of PandoraMonitoring:
==========================

PandoraSettings:
----------------
To visualize pandora objects, the "VisualMonitoring"-algorithm can be used. Just plug the following 
pandora-settings snippet into the PandoraSettings file at the location where the state of the pandora
reconstruction should be visualized. 

    <algorithm type = "VisualMonitoring" description = "display all">
        <!--  draw specific named cluster lists -->
        <ClusterListNames>PrimaryClusterList PhotonClusterList</ClusterListNames> 

        <!--  Use ROOT TEve for visualization -->
        <UseROOTEve>true</UseROOTEve>

        <!--  draw Monte Carlo particles -->
        <ShowMCParticles>true</ShowMCParticles> 

        <!--  draw the current PFOs -->
        <ShowCurrentPfos>true</ShowCurrentPfos>

        <!--  draw the current cluster-list   -->
        <ShowCurrentClusters>true</ShowCurrentClusters>

        <!--  draw the current track-list   -->
        <ShowCurrentTracks>false</ShowCurrentTracks>
        
        <!--  when drawing the current calo-hit-list, draw only calohits which are not clustered yet   -->
        <ShowOnlyAvailable>true</ShowOnlyAvailable>

        <!--  draw the current calohit-list   -->
        <ShowCurrentCaloHits>false</ShowCurrentCaloHits>

        <!--  (re)draw the display   -->
        <Show>true</Show>
    </algorithm>

The visualisation can be called at several places in the pandora settings file. The objects will 
be given to PandoraMonitoring and displayed when "Show" is set to true. Each time the 
visualisation algorithm is called with "Show" set to true a new event display will be created. 


In the algorithm code:
----------------------
For in-depth debugging of algorithms, the visualisation can be fed directly from within the source code
by using the Pandora monitoring API. An example is given in the following:

    PANDORA_MONITORING_API(VisualizeClusters(pClusterList, "currentList", AUTO, false, true  ) );
    PANDORA_MONITORING_API(ViewEvent() );

where pClusterList is a pointer to a cluster list (the current cluster list in the given example), 
"currentList" is the name which will be displayed in TEve. With AUTO the automatic color scheme is 
selected (look into PandoraMonitoringAPI.h for the full list of available colors). The two boolean 
variables define if the arrow indicating the linear fit of the calo hits shall be drawn and if the 
tracks associated to the clusters are drawn. 

With "ViewEvent()" the event display is redrawn. 

The full list of available visualisation commands can be inspected in PandoraMonitoringApi.h.
