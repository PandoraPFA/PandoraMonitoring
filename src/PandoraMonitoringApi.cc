/**
 *  @file   PandoraMonitoring/src/PandoraMonitoringApi.cc
 * 
 *  @brief  Redirection for pandora monitoring api class to the pandora monitoring implementation.
 * 
 *  $Log: $
 */

#include "PandoraMonitoring.h"
#include "PandoraMonitoringApi.h"

using namespace pandora_monitoring;

template <typename T>
void PandoraMonitoringApi::SetTreeVariable(const pandora::Pandora &pandora, const std::string &treeName, const std::string &variableName, T t)
{
    PandoraMonitoring::GetInstance(pandora)->SetTreeVariable(treeName, variableName, t);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::PrintTree(const pandora::Pandora &pandora, const std::string &treeName)
{
    PandoraMonitoring::GetInstance(pandora)->PrintTree(treeName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::ScanTree(const pandora::Pandora &pandora, const std::string &treeName)
{
    PandoraMonitoring::GetInstance(pandora)->ScanTree(treeName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::FillTree(const pandora::Pandora &pandora, const std::string &treeName)
{
    PandoraMonitoring::GetInstance(pandora)->FillTree(treeName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::SaveTree(const pandora::Pandora &pandora, const std::string &treeName, const std::string &fileName,
    const std::string &fileOptions)
{
    PandoraMonitoring::GetInstance(pandora)->SaveTree(treeName, fileName, fileOptions);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::Pandora &pandora, const T &t)
{
    PandoraMonitoring::GetInstance(pandora)->DrawPandoraHistogram(t, "");
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::Pandora &pandora, const T &t, const std::string &options)
{
    PandoraMonitoring::GetInstance(pandora)->DrawPandoraHistogram(t, options);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::SetEveDisplayParameters(const pandora::Pandora &pandora, const bool showDetectors,
    const DetectorView detectorView, const float transparencyThresholdE, const float energyScaleThresholdE)
{
    PandoraMonitoring::GetInstance(pandora)->SetEveDisplayParameters(showDetectors, detectorView, transparencyThresholdE, energyScaleThresholdE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeMCParticles(const pandora::Pandora &pandora, const pandora::MCParticleList *const pMCParticleList,
    const std::string &name, const Color color, const PdgCodeToEnergyMap *pParticleSuppressionMap)
{
    PandoraMonitoring::GetInstance(pandora)->VisualizeMCParticles(pMCParticleList, name, NULL, color, pParticleSuppressionMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeTracks(const pandora::Pandora &pandora, const pandora::TrackList *const pTrackList,
    const std::string &name, const Color color)
{
    PandoraMonitoring::GetInstance(pandora)->VisualizeTracks(pTrackList, name, NULL, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeCaloHits(const pandora::Pandora &pandora, const pandora::CaloHitList *const pCaloHitList,
    const std::string &name, const Color color)
{
    PandoraMonitoring::GetInstance(pandora)->VisualizeCaloHits(pCaloHitList, name, NULL, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeClusters(const pandora::Pandora &pandora, const pandora::ClusterList *const pClusterList,
    const std::string &name, const Color color, bool showAssociatedTracks)
{
    PandoraMonitoring::GetInstance(pandora)->VisualizeClusters(pClusterList, name, NULL, color, showAssociatedTracks);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeParticleFlowObjects(const pandora::Pandora &pandora, const pandora::PfoList *const pPfoList,
    const std::string &name, const Color color, bool showVertices, bool displayPfoHierarchy)
{
    PandoraMonitoring::GetInstance(pandora)->VisualizeParticleFlowObjects(pPfoList, name, NULL, color, showVertices, displayPfoHierarchy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeVertices(const pandora::Pandora &pandora, const pandora::VertexList *const pVertexList,
    const std::string &name, const Color color)
{
    PandoraMonitoring::GetInstance(pandora)->VisualizeVertices(pVertexList, name, NULL, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::AddMarkerToVisualization(const pandora::Pandora &pandora, const pandora::CartesianVector *const pMarkerPoint,
    const std::string &name, const Color color, const unsigned int markerSize)
{
    PandoraMonitoring::GetInstance(pandora)->AddMarkerToVisualization(pMarkerPoint, name, NULL, color, markerSize);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::ViewEvent(const pandora::Pandora &pandora)
{
    PandoraMonitoring::GetInstance(pandora)->ViewEvent();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Pause(const pandora::Pandora &pandora)
{
    PandoraMonitoring::GetInstance(pandora)->Pause();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Create(const pandora::Pandora &pandora)
{
    (void) PandoraMonitoring::GetInstance(pandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Delete(const pandora::Pandora &pandora)
{
    PandoraMonitoring::GetInstance(pandora)->DeleteInstance(pandora);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template void PandoraMonitoringApi::SetTreeVariable(const pandora::Pandora &, const std::string&, const std::string&, float);
template void PandoraMonitoringApi::SetTreeVariable(const pandora::Pandora &, const std::string&, const std::string&, int);
template void PandoraMonitoringApi::SetTreeVariable(const pandora::Pandora &, const std::string&, const std::string&, double);

template void PandoraMonitoringApi::SetTreeVariable(const pandora::Pandora &, const std::string&, const std::string&, std::vector<float> *);
template void PandoraMonitoringApi::SetTreeVariable(const pandora::Pandora &, const std::string&, const std::string&, std::vector<int> *);
template void PandoraMonitoringApi::SetTreeVariable(const pandora::Pandora &, const std::string&, const std::string&, std::vector<double> *);

template void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::Pandora &, const pandora::Histogram &);
template void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::Pandora &, const pandora::TwoDHistogram &);

template void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::Pandora &, const pandora::Histogram &, const std::string &);
template void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::Pandora &, const pandora::TwoDHistogram &, const std::string &);
