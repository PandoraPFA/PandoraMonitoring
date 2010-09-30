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

void PandoraMonitoringApi::Create1DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp, 
    const std::string xAxisTitle, const std::string yAxisTitle)
{
    PandoraMonitoring::GetInstance()->Create1DHistogram(name, title, nBinsX, xLow, xUp, xAxisTitle, yAxisTitle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Create2DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp,
    int nBinsY, double yLow, double yUp, const std::string xAxisTitle, const std::string yAxisTitle)
{
    PandoraMonitoring::GetInstance()->Create2DHistogram(name, title, nBinsX, xLow, xUp, nBinsY, yLow, yUp, xAxisTitle, yAxisTitle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Fill1DHistogram(const std::string &name, float xValue, float weight)
{
    PandoraMonitoring::GetInstance()->Fill1DHistogram(name, xValue, weight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Fill2DHistogram(const std::string &name, float xValue, float yValue, float weight)
{
    PandoraMonitoring::GetInstance()->Fill2DHistogram(name, xValue, yValue, weight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::AddMultiplyOrDivideHistograms(const std::string &nameHisto0, const std::string &nameHisto1, 
				       double coeff0, double coeff1,
				       bool add, bool multiply )
{
    PandoraMonitoring::GetInstance()->AddMultiplyOrDivideHistograms(nameHisto0, nameHisto1, coeff0, coeff1, add, multiply );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DrawHistogram(const std::string &name)
{
    PandoraMonitoring::GetInstance()->DrawHistogram(name, "");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DrawHistogram(const std::string &name, const std::string &options)
{
    PandoraMonitoring::GetInstance()->DrawHistogram(name, options);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::SaveAndCloseHistogram(const std::string &name, const std::string &fileName, const std::string &fileOptions)
{
    PandoraMonitoring::GetInstance()->SaveAndCloseHistogram(name, fileName, fileOptions);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DeleteHistogram(const std::string &name)
{
    PandoraMonitoring::GetInstance()->DeleteHistogram(name);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void PandoraMonitoringApi::SetTreeVariable(const std::string &treeName, const std::string &variableName, T t)
{
    PandoraMonitoring::GetInstance()->SetTreeVariable(treeName, variableName, t);
}

// instantiations of this template member function for the permitted types
template void PandoraMonitoringApi::SetTreeVariable(const std::string&, const std::string&, float  t);
template void PandoraMonitoringApi::SetTreeVariable(const std::string&, const std::string&, int    t);
template void PandoraMonitoringApi::SetTreeVariable(const std::string&, const std::string&, double t);

template void PandoraMonitoringApi::SetTreeVariable(const std::string&, const std::string&, std::vector<float>*  t);
template void PandoraMonitoringApi::SetTreeVariable(const std::string&, const std::string&, std::vector<int>*    t);
template void PandoraMonitoringApi::SetTreeVariable(const std::string&, const std::string&, std::vector<double>* t);

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::PrintTree(const std::string &treeName)
{
    PandoraMonitoring::GetInstance()->PrintTree(treeName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::ScanTree(const std::string &treeName)
{
    PandoraMonitoring::GetInstance()->ScanTree(treeName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::FillTree(const std::string &treeName)
{
    PandoraMonitoring::GetInstance()->FillTree(treeName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::SaveTree(const std::string &treeName, const std::string &fileName, const std::string &fileOptions)
{
    PandoraMonitoring::GetInstance()->SaveTree(treeName, fileName, fileOptions);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::ViewEvent()
{
    PandoraMonitoring::GetInstance()->ViewEvent();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeMCParticles(const pandora::MCParticleList *const pMCParticleList, std::string name, Color color,
    const PdgCodeToEnergyMap *pParticleSuppressionMap)
{
    PandoraMonitoring::GetInstance()->VisualizeMCParticles(pMCParticleList, name, NULL, color, pParticleSuppressionMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeTracks(const pandora::TrackList *const pTrackList, std::string name, Color color )
{
    PandoraMonitoring::GetInstance()->VisualizeTracks(pTrackList, name, NULL, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeCaloHits(const pandora::OrderedCaloHitList *const pOrderedCaloHitList, std::string name, Color color )
{
    PandoraMonitoring::GetInstance()->VisualizeCaloHits(pOrderedCaloHitList, name, NULL, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeClusters(const pandora::ClusterList *const pClusterList, std::string name, Color color, 
    bool showAssociatedTracks, bool showFit)
{
    PandoraMonitoring::GetInstance()->VisualizeClusters(pClusterList, name, NULL, color, showAssociatedTracks, showFit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeParticleFlowObjects(const pandora::ParticleFlowObjectList *const pPfoList, std::string name,
    Color color, bool showAssociatedTracks, bool showFit)
{
    PandoraMonitoring::GetInstance()->VisualizeParticleFlowObjects(pPfoList, name, NULL, color, showAssociatedTracks, showFit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Delete()
{
    PandoraMonitoring::GetInstance()->DeleteInstance();
}
