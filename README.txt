README
======

Recompilation of PandoraPFANew with PandoraMonitoring support:
==============================================================

If using cmake
cmake -DPANDORA_MONITORING=ON OTHER_OPTIONS ..
make install

If using standalone makefile
make MONITORING=1 OTHER_OPTIONS


Usage of PandoraMonitoring:
==========================

PandoraSettings:
----------------
To visualize pandora objects, the "VisualMonitoring"-algorithm can be used. Just plug the following
pandora-settings snippet into the PandoraSettings file at the location where the state of the pandora
reconstruction should be visualized.

    <algorithm type = "VisualMonitoring" description = "display all">
        <ShowMCParticles>false</ShowMCParticles>
        <ShowCurrentCaloHits>false</ShowCurrentCaloHits>
        <ShowCurrentTracks>false</ShowCurrentTracks>
        <ShowCurrentClusters>true</ShowCurrentClusters>
        <ShowCurrentPfos>true</ShowCurrentPfos>
        <HitColors>pfo</HitColors>
        <ShowOnlyAvailable>false</ShowOnlyAvailable>
        <DisplayEvent>true</DisplayEvent>
        <TransparencyThresholdE>-1.</TransparencyThresholdE>
        <EnergyScaleThresholdE>1.</EnergyScaleThresholdE>
        <BlackBackground>false</BlackBackground>
        <ShowDetector>true</ShowDetector>
    </algorithm>

The visualisation can be called at several places in the pandora settings file. The objects will
be given to PandoraMonitoring and displayed when "DisplayEvent" is set to true. Each time the
visualisation algorithm is called with "DisplayEvent" set to true a new event display will be created.


In the algorithm code:
----------------------
For in-depth debugging of algorithms, the visualisation can be fed directly from within the source code
by using the Pandora monitoring API. An example, for use in an algorithm or a Pandora "process", is:

    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), pClusterList, "listName", AUTO));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

where pClusterList the address of a cluster list (the current cluster list in the given example),
"listName" is the name which will be displayed in TEve. With AUTO the automatic color scheme is
selected (look into PandoraMonitoringAPI.h for the full list of available colors). The two boolean
variables define if the arrow indicating the linear fit of the calo hits shall be drawn and if the
tracks associated to the clusters are drawn. 

With "DisplayEvent()" the event display is redrawn. 

The full list of available visualisation commands can be inspected in PandoraMonitoringApi.h.
