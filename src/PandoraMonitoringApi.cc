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

void PandoraMonitoringApi::Create1DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp)
{
    PandoraMonitoring::GetInstance()->Create1DHistogram(name, title, nBinsX, xLow, xUp);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Create2DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp,
    int nBinsY, double yLow, double yUp)
{
    PandoraMonitoring::GetInstance()->Create2DHistogram(name, title, nBinsX, xLow, xUp, nBinsY, yLow, yUp);
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

void PandoraMonitoringApi::DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList)
{
    PandoraMonitoring::GetInstance()->DrawEvent(detectorView, pTrackList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DrawEvent(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList)
{
    PandoraMonitoring::GetInstance()->DrawEvent(detectorView, pOrderedCaloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DrawEvent(DetectorView detectorView, const pandora::ClusterList *const pClusterList)
{
    PandoraMonitoring::GetInstance()->DrawEvent(detectorView, pClusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList,
    const pandora::OrderedCaloHitList *const pOrderedCaloHitList)
{
    PandoraMonitoring::GetInstance()->DrawEvent(detectorView, pTrackList, pOrderedCaloHitList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList,
    const pandora::ClusterList *const pClusterList)
{
    PandoraMonitoring::GetInstance()->DrawEvent(detectorView, pTrackList, pClusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DrawDetectorOutline(DetectorView detectorView)
{
    PandoraMonitoring::GetInstance()->DetectorOutlineTest(detectorView);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Test()
{
    PandoraMonitoring::GetInstance()->DrawCanvas();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Test2(const pandora::ClusterList *const pClusterList)
{
    PandoraMonitoring::GetInstance()->LookAtClusters(pClusterList);
}
