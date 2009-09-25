/**
 *  @file   PandoraMonitoring/src/PandoraMonitoring.cc
 * 
 *  @brief  Implementation of the pandora monitoring class.
 * 
 *  $Log: $
 */

// Pandora include files
#include "Helpers/GeometryHelper.h"

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/OrderedCaloHitList.h"

// ROOT include files
#include "TArc.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPad.h"
#include "TPolyLine.h"
#include "TSystem.h"
#include "TStyle.h"

#include "PandoraMonitoring.h"

#include <cmath>
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
        gStyle->SetPalette(1);
    }

    return m_pPandoraMonitoring;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::Reset()
{
    if(m_instanceFlag)
    {
        delete m_pPandoraMonitoring;
        m_pPandoraMonitoring = NULL;
        m_instanceFlag = false;

        for (HistogramMap::iterator iter = m_histogramMap.begin(), iterEnd = m_histogramMap.end(); iter != iterEnd; ++iter)
            delete iter->second;

        m_histogramMap.clear();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::Create1DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp)
{
    if (m_histogramMap.end() != m_histogramMap.find(name))
    {
        std::cout << "PandoraMonitoring::Create1DHistogram, error: Histogram with this name already exists." << std::endl;
        throw std::exception();
    }

    TH1F* pTH1F = new TH1F(name.c_str(), title.c_str(), nBinsX, xLow, xUp);
    m_histogramMap.insert(HistogramMap::value_type(name, pTH1F));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::Create2DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp, int nBinsY,
    double yLow, double yUp)
{
    if (m_histogramMap.end() != m_histogramMap.find(name))
    {
        std::cout << "PandoraMonitoring::Create2DHistogram, error: Histogram with this name already exists." << std::endl;
        throw std::exception();
    }

    TH2F* pTH2F = new TH2F(name.c_str(), title.c_str(), nBinsX, xLow, xUp, nBinsY, yLow, yUp);
    m_histogramMap.insert(HistogramMap::value_type(name, pTH2F));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::Fill1DHistogram(const std::string &name, float xValue, float weight)
{
    HistogramMap::iterator iter = m_histogramMap.find(name);

    if (m_histogramMap.end() == iter)
    {
        std::cout << "PandoraMonitoring::Fill1DHistogram, error: No histogram with this name exists." << std::endl;
        throw std::exception();
    }

    iter->second->Fill(xValue, weight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::Fill2DHistogram(const std::string &name, float xValue, float yValue, float weight)
{
    HistogramMap::iterator iter = m_histogramMap.find(name);

    if (m_histogramMap.end() == iter)
    {
        std::cout << "PandoraMonitoring::Fill2DHistogram, error: No histogram with this name exists." << std::endl;
        throw std::exception();
    }

    TH2F *pTH2F = static_cast<TH2F *>(iter->second);

    pTH2F->Fill(xValue, yValue, weight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawHistogram(const std::string &name, const std::string &options) const
{
    HistogramMap::const_iterator iter = m_histogramMap.find(name);

    if (m_histogramMap.end() == iter)
    {
        std::cout << "PandoraMonitoring::DrawHistogram, error: No histogram with this name exists." << std::endl;
        throw std::exception();
    }

    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    iter->second->Draw(options.c_str());
    this->Pause();

    delete pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::SaveHistogram(const std::string &name, const std::string &fileName)
{
    HistogramMap::iterator iter = m_histogramMap.find(name);

    if (m_histogramMap.end() == iter)
    {
        std::cout << "PandoraMonitoring::SaveHistogram, error: No histogram with this name exists." << std::endl;
        throw std::exception();
    }

    TFile* pTFile = new TFile(fileName.c_str(), "recreate");

    iter->second->SetDirectory(pTFile);
    iter->second->Write(name.c_str(), TObject::kOverwrite);

    pTFile->Write();
    pTFile->Close();
    m_histogramMap.erase(iter);
    delete pTFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DeleteHistogram(const std::string &name)
{
    HistogramMap::iterator iter = m_histogramMap.find(name);

    if (m_histogramMap.end() == iter)
    {
        std::cout << "PandoraMonitoring::DeleteHistogram, error: No histogram with this name exists." << std::endl;
        throw std::exception();
    }

    delete iter->second;
    m_histogramMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawDetectorOutline(DetectorView detectorView)
{
    if (!m_isOutlineConstructed)
        this->MakeDetectorOutline();

    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    TH2F *pAxesHist = NULL;

    switch(detectorView)
    {
    case DETECTOR_VIEW_XY:
        pAxesHist = new TH2F("xyview", "x-y view", 100, -4200, 4200, 100, -4200, 4200);
        pAxesHist->SetFillColor(kWhite);
        pAxesHist->SetStats(kFALSE);
        pAxesHist->Draw();

        for(unsigned int i = 0; i <m_2DLinesXY.size(); ++i)
        {
            m_2DLinesXY[i]->SetLineWidth(2);
            m_2DLinesXY[i]->Draw();
        }

        for(unsigned int i = 0; i < m_2DCirclesXY.size(); ++i)
        {
            m_2DCirclesXY[i]->SetLineWidth(2);
            m_2DCirclesXY[i]->Draw();
        }
        break;

    case DETECTOR_VIEW_XZ:
        pAxesHist = new TH2F("xzview", "x-z view", 100, -4200, 4200, 100, -4200, 4200);
        pAxesHist->SetFillColor(kWhite);
        pAxesHist->SetStats(kFALSE);
        pAxesHist->Draw();

        for(unsigned int i = 0; i < m_2DLinesXZ.size(); ++i)
        {
            m_2DLinesXZ[i]->SetLineWidth(2);
            m_2DLinesXZ[i]->Draw();
        }
        break;

    default:
        std::cout << "PandoraMonitoring::DrawDetectorOutline, error: request for an unsupported detector view." << std::endl;
        throw std::exception();
    }

    this->Pause();
    delete pCanvas;

    if (NULL != pAxesHist)
        delete pAxesHist;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawCanvas() const
{
    TH1F* pHist = new TH1F("PandoraMonitoringHist", "PandoraMonitoringHist", 100, 0, 1);

    for (int i = 0; i < 1000; ++i)
    {
        pHist->Fill(static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX));
    }

    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    pHist->Draw();
    this->Pause();

    delete pCanvas;
    delete pHist;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::LookAtClusters(const pandora::ClusterList *const pClusterList) const
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

void PandoraMonitoring::Pause() const
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

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeDetectorOutline()
{
    pandora::GeometryHelper *pGeometryHelper = pandora::GeometryHelper::GetInstance();

    // HCal, ECal (Barrel and EndCap)
    const pandora::GeometryHelper::SubDetectorParameters eCalBarrelParameters = pGeometryHelper->GetECalBarrelParameters();
    const pandora::GeometryHelper::SubDetectorParameters eCalEndCapParameters = pGeometryHelper->GetECalEndCapParameters();
    const pandora::GeometryHelper::SubDetectorParameters hCalBarrelParameters = pGeometryHelper->GetHCalBarrelParameters();
    const pandora::GeometryHelper::SubDetectorParameters hCalEndCapParameters = pGeometryHelper->GetHCalEndCapParameters();

    this->MakeXYLayerOutline(eCalBarrelParameters.GetInnerSymmetryOrder(), eCalBarrelParameters.GetInnerPhiCoordinate(), eCalBarrelParameters.GetInnerRCoordinate(), 1, 1);
    this->MakeXYLayerOutline(eCalBarrelParameters.GetOuterSymmetryOrder(), eCalBarrelParameters.GetOuterPhiCoordinate(), eCalBarrelParameters.GetOuterRCoordinate(), 1, 1);
    this->MakeXYLayerOutline(hCalBarrelParameters.GetInnerSymmetryOrder(), hCalBarrelParameters.GetInnerPhiCoordinate(), hCalBarrelParameters.GetInnerRCoordinate(), 1, 1);
    this->MakeXYLayerOutline(hCalBarrelParameters.GetOuterSymmetryOrder(), hCalBarrelParameters.GetOuterPhiCoordinate(), hCalBarrelParameters.GetOuterRCoordinate(), 1, 1);

    this->MakeXZLayerOutline(eCalBarrelParameters.GetInnerRCoordinate(), eCalBarrelParameters.GetOuterRCoordinate(), eCalBarrelParameters.GetInnerZCoordinate(), eCalBarrelParameters.GetOuterZCoordinate(), 1, 1);
    this->MakeXZLayerOutline(hCalBarrelParameters.GetInnerRCoordinate(), hCalBarrelParameters.GetOuterRCoordinate(), hCalBarrelParameters.GetInnerZCoordinate(), hCalBarrelParameters.GetOuterZCoordinate(), 1, 1);
    this->MakeXZLayerOutline(eCalEndCapParameters.GetInnerRCoordinate(), eCalEndCapParameters.GetOuterRCoordinate(), eCalEndCapParameters.GetInnerZCoordinate(), eCalEndCapParameters.GetOuterZCoordinate(), 1, 1);
    this->MakeXZLayerOutline(hCalEndCapParameters.GetInnerRCoordinate(), hCalEndCapParameters.GetOuterRCoordinate(), hCalEndCapParameters.GetInnerZCoordinate(), hCalEndCapParameters.GetOuterZCoordinate(), 1, 1);

    // TPC
    this->MakeXYLayerOutline(0, 0, pGeometryHelper->GetMainTrackerInnerRadius(), 1, 1);
    this->MakeXYLayerOutline(0, 0, pGeometryHelper->GetMainTrackerOuterRadius(), 1, 1);
    this->MakeXZLayerOutline(pGeometryHelper->GetMainTrackerInnerRadius(), pGeometryHelper->GetMainTrackerOuterRadius(), 0, pGeometryHelper->GetMainTrackerZExtent(), 1, 1);

    // Additional subdetectors
    for (pandora::GeometryHelper::SubDetectorParametersList::const_iterator iter = pGeometryHelper->GetAdditionalSubDetectors().begin(),
        iterEnd = pGeometryHelper->GetAdditionalSubDetectors().end(); iter != iterEnd; ++iter)
    {
        this->MakeXYLayerOutline(iter->GetInnerSymmetryOrder(), iter->GetInnerPhiCoordinate(), iter->GetInnerRCoordinate(), 1, 1);
        this->MakeXYLayerOutline(iter->GetOuterSymmetryOrder(), iter->GetOuterPhiCoordinate(), iter->GetOuterRCoordinate(), 1, 1);
        this->MakeXZLayerOutline(iter->GetInnerRCoordinate(), iter->GetOuterRCoordinate(), iter->GetInnerZCoordinate(), iter->GetOuterZCoordinate(), 1, 1);
    }

    m_isOutlineConstructed = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeXYLayerOutline(int symmetryOrder, float phi0, float closestDistanceToIp, int lineWidth, int lineColor)
{
    if (symmetryOrder > 2)
    {
        float *pXCoordinate = new float[symmetryOrder + 1];
        float *pYCoordinate = new float[symmetryOrder + 1];

        const float pi(3.1415927);
        const float x0(-1. * closestDistanceToIp * tan(pi / float(symmetryOrder)));
        const float y0(closestDistanceToIp);

        for (int i = 0; i <= symmetryOrder; ++i)
        {
            const float theta(phi0 + (2 * pi * float(i) / float(symmetryOrder)));

            pXCoordinate[i] = x0 * cos(theta) + y0 * sin(theta);
            pYCoordinate[i] = y0 * cos(theta) - x0 * sin(theta);
        }

        TPolyLine *pTPolyLine = new TPolyLine(symmetryOrder + 1, pXCoordinate, pYCoordinate);
        pTPolyLine->SetLineWidth(lineWidth);
        pTPolyLine->SetLineColor(lineColor);

        m_2DLinesXY.push_back(pTPolyLine);

        delete [] pXCoordinate;
        delete [] pYCoordinate;
    }
    else if (symmetryOrder == 0)
    {
        TArc* pTArc = new TArc(0., 0., closestDistanceToIp);
        pTArc->SetLineWidth(lineWidth);
        pTArc->SetLineColor(lineColor);

        m_2DCirclesXY.push_back(pTArc);
    }
     else
     {
         std::cout << "PandoraMonitoring::MakeXYLayerOutline, symmetryOrder is " << symmetryOrder << std::endl;
     }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeXZLayerOutline(float innerRCoordinate, float outerRCoordinate, float innerZCoordinate, float outerZCoordinate,
    int lineWidth, int lineColor)
{

    for (int i = 0; i< 4; ++i)
    {
        float x[5], z[5];

        switch(i)
        {
        case 0:
            x[0] = innerZCoordinate; x[1] = outerZCoordinate; x[2] = outerZCoordinate; x[3] = innerZCoordinate; x[4] = x[0];
            z[0] = innerRCoordinate; z[1] = innerRCoordinate; z[2] = outerRCoordinate; z[3] = outerRCoordinate; z[4] = z[0];
            break;
        case 1:
            x[0] =  innerZCoordinate; x[1] =  outerZCoordinate; x[2] =  outerZCoordinate; x[3] =  innerZCoordinate; x[4] = x[0];
            z[0] = -innerRCoordinate; z[1] = -innerRCoordinate; z[2] = -outerRCoordinate; z[3] = -outerRCoordinate; z[4] = z[0];
            break;
        case 2:
            x[0] = -innerZCoordinate; x[1] = -outerZCoordinate; x[2] = -outerZCoordinate; x[3] = -innerZCoordinate; x[4] = x[0];
            z[0] =  innerRCoordinate; z[1] =  innerRCoordinate; z[2] =  outerRCoordinate; z[3] =  outerRCoordinate; z[4] = z[0];
            break;
        case 3:
            x[0] = -innerZCoordinate; x[1] = -outerZCoordinate; x[2] = -outerZCoordinate; x[3] = -innerZCoordinate; x[4] = x[0];
            z[0] = -innerRCoordinate; z[1] = -innerRCoordinate; z[2] = -outerRCoordinate; z[3] = -outerRCoordinate; z[4] = z[0];
            break;
        default:
            throw std::exception();
        }

        int linesToDisplay = 5;

        if (0 == innerZCoordinate)
            linesToDisplay = 4;

        TPolyLine *pTPolyLine = new TPolyLine(linesToDisplay, x, z);
        pTPolyLine->SetLineWidth(lineWidth);
        pTPolyLine->SetLineColor(lineColor);

        m_2DLinesXZ.push_back(pTPolyLine);
    }
}

} // namespace pandora_monitoring
