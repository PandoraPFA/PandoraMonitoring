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

        for(unsigned int i = 0; i < m_2DBoxesXZ.size(); ++i)
        {
            m_2DBoxesXZ[i]->SetLineWidth(2);
            m_2DBoxesXZ[i]->SetLineColor(1);
            m_2DBoxesXZ[i]->SetFillStyle(0);
            m_2DBoxesXZ[i]->Draw();
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

    const pandora::GeometryHelper::SubDetectorParameters eCalBarrelParameters = pGeometryHelper->GetECalBarrelParameters();
    const pandora::GeometryHelper::SubDetectorParameters eCalEndCapParameters = pGeometryHelper->GetECalEndCapParameters();
    const pandora::GeometryHelper::SubDetectorParameters hCalBarrelParameters = pGeometryHelper->GetHCalBarrelParameters();
    const pandora::GeometryHelper::SubDetectorParameters hCalEndCapParameters = pGeometryHelper->GetHCalEndCapParameters();

    this->MakeLayerOutline(DETECTOR_VIEW_XY, eCalBarrelParameters.GetInnerSymmetryOrder(), eCalBarrelParameters.GetInnerPhiCoordinate(),
        eCalBarrelParameters.GetInnerRCoordinate(), 1, 1);
    this->MakeLayerOutline(DETECTOR_VIEW_XY, eCalBarrelParameters.GetOuterSymmetryOrder(), eCalBarrelParameters.GetOuterPhiCoordinate(),
        eCalBarrelParameters.GetOuterRCoordinate(), 1, 1);
    this->MakeLayerOutline(DETECTOR_VIEW_XY, hCalBarrelParameters.GetInnerSymmetryOrder(), hCalBarrelParameters.GetInnerPhiCoordinate(),
        hCalBarrelParameters.GetInnerRCoordinate(), 1, 1);
    this->MakeLayerOutline(DETECTOR_VIEW_XY, hCalBarrelParameters.GetOuterSymmetryOrder(), hCalBarrelParameters.GetOuterPhiCoordinate(),
        hCalBarrelParameters.GetOuterRCoordinate(), 1, 1);

    m_isOutlineConstructed = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeLayerOutline(DetectorView detectorView, int symmetryOrder, float phi0, float closestDistanceToIp,
    int lineWidth, int lineColor)
{
    if (symmetryOrder > 2)
    {
        float *pCoordinate1 = new float[symmetryOrder + 1];
        float *pCoordinate2 = new float[symmetryOrder + 1];

        const float pi(3.1415927);
        const float magnitude(2 * closestDistanceToIp * tan(pi / float(symmetryOrder)));

        pCoordinate1[0] = (-0.5 * magnitude); // TODO initial phi0
        pCoordinate2[0] = closestDistanceToIp;

        for (int i = 0; i < symmetryOrder; ++i)
        {
            pCoordinate1[i + 1] = pCoordinate1[i] + (magnitude * cos(phi0 + (2 * pi * float(i) / float(symmetryOrder))));
            pCoordinate2[i + 1] = pCoordinate2[i] - (magnitude * sin(phi0 + (2 * pi * float(i) / float(symmetryOrder))));
        }

        TPolyLine *pTPolyLine = new TPolyLine(symmetryOrder + 1, pCoordinate1, pCoordinate2);
        pTPolyLine->SetLineWidth(lineWidth);
        pTPolyLine->SetLineColor(lineColor);

        if (DETECTOR_VIEW_XY == detectorView)
        {
            m_2DLinesXY.push_back(pTPolyLine);
        }
        else if (DETECTOR_VIEW_XZ == detectorView)
        {
            m_2DLinesXZ.push_back(pTPolyLine);
        }
        else
        {
            std::cout << "PandoraMonitoring::MakeLayerOutline, error: request for an unsupported detector view." << std::endl;
            throw std::exception();
        }

        for (int i = 0; i < symmetryOrder + 1; ++i)
        {
            std::cout << " i " << i << " pCoordinate1[i] " << pCoordinate1[i] << " pCoordinate2[i] " << pCoordinate2[i] << std::endl;
        }

        delete [] pCoordinate1;
        delete [] pCoordinate2;
    }
    else if (symmetryOrder == 0)
    {
        TArc* pTArc = new TArc(0., 0., closestDistanceToIp);
        pTArc->SetLineWidth(lineWidth);
        pTArc->SetLineColor(lineColor);

        if (DETECTOR_VIEW_XY != detectorView)
        {
            std::cout << "PandoraMonitoring::MakeLayerOutline, error: Unsupported detector structure." << std::endl;
            throw std::exception();
        }

        m_2DCirclesXY.push_back(pTArc);
    }
    else
    {
        std::cout << "PandoraMonitoring::MakeLayerOutline, symmetryOrder is " << symmetryOrder << std::endl;
    }
}

} // namespace pandora_monitoring
