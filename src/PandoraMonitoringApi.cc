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

void PandoraMonitoringApi::AddHistograms(const std::string &nameHisto0, const std::string &nameHisto1, double coeff0, double coeff1)
{
    PandoraMonitoring::GetInstance()->AddMultiplyOrDivideHistograms(nameHisto0, nameHisto1, coeff0, coeff1, true, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::MultiplyHistograms(const std::string &nameHisto0, const std::string &nameHisto1, double coeff0, double coeff1)
{
    PandoraMonitoring::GetInstance()->AddMultiplyOrDivideHistograms(nameHisto0, nameHisto1, coeff0, coeff1, false, true);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DivideHistograms(const std::string &nameHisto0, const std::string &nameHisto1, double coeff0, double coeff1)
{
    PandoraMonitoring::GetInstance()->AddMultiplyOrDivideHistograms(nameHisto0, nameHisto1, coeff0, coeff1, false, false);
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

template <typename T>
void PandoraMonitoringApi::DrawPandoraHistogram(const T &t)
{
    PandoraMonitoring::GetInstance()->DrawPandoraHistogram(t, "");
}

// instantiations of this template member function for the permitted types
template void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::Histogram &histogram);
template void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::TwoDHistogram &twoDHistogram);

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void PandoraMonitoringApi::DrawPandoraHistogram(const T &t, const std::string &options)
{
    PandoraMonitoring::GetInstance()->DrawPandoraHistogram(t, options);
}

// instantiations of this template member function for the permitted types
template void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::Histogram &histogram, const std::string &options);
template void PandoraMonitoringApi::DrawPandoraHistogram(const pandora::TwoDHistogram &twoDHistogram, const std::string &options);

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

void PandoraMonitoringApi::SetEveDisplayParameters(const bool blackBackground, const bool showDetectors, const float transparencyThresholdE,
    const float energyScaleThresholdE)
{
    PandoraMonitoring::GetInstance()->SetEveDisplayParameters(blackBackground, showDetectors, transparencyThresholdE, energyScaleThresholdE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeMCParticles(const pandora::MCParticleList *const pMCParticleList, std::string name, Color color,
    const PdgCodeToEnergyMap *pParticleSuppressionMap)
{
    PandoraMonitoring::GetInstance()->VisualizeMCParticles(pMCParticleList, name, NULL, color, pParticleSuppressionMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeTracks(const pandora::TrackList *const pTrackList, std::string name, Color color)
{
    PandoraMonitoring::GetInstance()->VisualizeTracks(pTrackList, name, NULL, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeCaloHits(const pandora::CaloHitList *const pCaloHitList, std::string name, Color color)
{
    PandoraMonitoring::GetInstance()->VisualizeCaloHits(pCaloHitList, name, NULL, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeClusters(const pandora::ClusterList *const pClusterList, std::string name, Color color,
    bool showAssociatedTracks)
{
    PandoraMonitoring::GetInstance()->VisualizeClusters(pClusterList, name, NULL, color, showAssociatedTracks);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::VisualizeParticleFlowObjects(const pandora::PfoList *const pPfoList, std::string name,
    Color color, bool showAssociatedTracks)
{
    PandoraMonitoring::GetInstance()->VisualizeParticleFlowObjects(pPfoList, name, NULL, color, showAssociatedTracks);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::AddMarkerToVisualization(const pandora::CartesianVector *const pMarkerPoint, std::string name, Color color,
    const unsigned int markerSize)
{
    PandoraMonitoring::GetInstance()->AddMarkerToVisualization(pMarkerPoint, name, color, markerSize);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Pause()
{
    PandoraMonitoring::GetInstance()->Pause();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Create()
{
    (void) PandoraMonitoring::GetInstance();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Delete()
{
    PandoraMonitoring::GetInstance()->DeleteInstance();
}
