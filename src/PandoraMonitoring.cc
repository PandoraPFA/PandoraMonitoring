/**
 *  @file   PandoraMonitoring/src/PandoraMonitoring.cc
 * 
 *  @brief  Implementation of the pandora monitoring class.
 * 
 *  $Log: $
 */

#include "PandoraMonitoring.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TPad.h"
#include "TSystem.h"

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
    gSystem->ProcessEvents();

    this->Pause();

    delete pCanvas;
    delete pHist;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::Pause()
{
    std::cout << "Press return to continue... ";
    getchar();
    std::cout << std::endl;
}

} // namespace pandora_monitoring
