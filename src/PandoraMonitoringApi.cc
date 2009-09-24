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

void PandoraMonitoringApi::SaveHistogram(const std::string &name, const std::string &fileName)
{
    PandoraMonitoring::GetInstance()->SaveHistogram(name, fileName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DeleteHistogram(const std::string &name)
{
    PandoraMonitoring::GetInstance()->DeleteHistogram(name);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::DrawDetectorOutline(DetectorView detectorView)
{
    PandoraMonitoring::GetInstance()->DrawDetectorOutline(detectorView);
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
