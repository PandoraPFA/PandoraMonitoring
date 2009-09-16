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

void PandoraMonitoringApi::Test()
{
    PandoraMonitoring::GetInstance()->DrawCanvas();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoringApi::Test2(const pandora::ClusterList *const pClusterList)
{
    PandoraMonitoring::GetInstance()->LookAtClusters(pClusterList);
}
