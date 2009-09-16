/**
 *  @file   PandoraMonitoring/src/PandoraMonitoring.cc
 * 
 *  @brief  Implementation of the pandora monitoring class.
 * 
 *  $Log: $
 */

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/OrderedCaloHitList.h"

#include "PandoraMonitoring.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TPad.h"
#include "TSystem.h"

#include <fcntl.h>

namespace pandora_monitoring
{

bool PandoraMonitoring::m_instanceFlag = false;

PandoraMonitoring* PandoraMonitoring::m_pPandoraMonitoring = NULL;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraMonitoring *PandoraMonitoring::GetInstance()
{
    if(!m_instanceFlag)
    {
        m_pPandoraMonitoring = new PandoraMonitoring();
        m_instanceFlag = true;
    }

    return m_pPandoraMonitoring;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawCanvas()
{
    TCanvas *pCanvas = new TCanvas("MyRecoDisplay", "MyRecoDisplay", 750, 750);
    TH1F* pHist = new TH1F("hist", "hist", 100, 0, 1);

    for (int i = 0; i < 1000; ++i)
    {
        pHist->Fill(static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX));
    }

    pCanvas->Draw();
    pHist->Draw();
    this->Pause();

    delete pCanvas;
    delete pHist;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::LookAtClusters(const pandora::ClusterList *const pClusterList)
{
    for (pandora::ClusterList::const_iterator clusterIter = pClusterList->begin(); clusterIter != pClusterList->end(); ++clusterIter)
    {
        std::cout << "Monitoring, Cluster! " << *clusterIter << std::endl;

        const pandora::OrderedCaloHitList orderedCaloHitList((*clusterIter)->GetOrderedCaloHitList());

        for (pandora::OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
        {
            for (pandora::CaloHitList::const_iterator caloHitIter = iter->second->begin(), caloHitIterEnd = iter->second->end(); caloHitIter != caloHitIterEnd; ++caloHitIter)
            {
                std::cout << "calo hit: " << iter->first << ", " << (*caloHitIter)->GetParentCaloHitAddress() << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::Pause()
{
#ifdef __unix__
    std::cout << "Press return to continue ..." << std::endl;
    int flag = fcntl(1, F_GETFL, 0);

    int key = 0;
    while(true)
    {
        gSystem->ProcessEvents();
        fcntl(1, F_SETFL, flag | O_NONBLOCK);
        key = getchar();

        if((key == '\n') || (key == '\r'))
            break;

        usleep(1000);
    }

    fcntl(1, F_SETFL, flag);
#else
    std::cout << "PandoraMonitoring::Pause() is only implemented for unix operating systems." << std::endl;
#endif
}

} // namespace pandora_monitoring
