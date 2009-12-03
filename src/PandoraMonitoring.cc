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
#include "Objects/CartesianVector.h"
#include "Objects/Cluster.h"
#include "Objects/OrderedCaloHitList.h"
#include "Objects/Track.h"

// LCIO includes
#include "EVENT/Track.h"

// ROOT include files
#include "TArc.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPolyMarker.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TTree.h"

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
        std::cout << "PandoraMonitoring::Create1DHistogram, error: Histogram with name '"<< name <<"' already exists." << std::endl;
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
        std::cout << "PandoraMonitoring::Create2DHistogram, error: Histogram with name '"<< name <<"' already exists." << std::endl;
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
        std::cout << "PandoraMonitoring::Fill1DHistogram, error: No histogram with name '"<< name <<"' exists." << std::endl;
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
        std::cout << "PandoraMonitoring::Fill2DHistogram, error: No histogram with name '"<< name <<"' exists." << std::endl;
        throw std::exception();
    }

    TH2F *pTH2F = dynamic_cast<TH2F *>(iter->second);

    if (NULL == pTH2F)
        throw std::exception();

    pTH2F->Fill(xValue, yValue, weight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawHistogram(const std::string &name, const std::string &options) const
{
    HistogramMap::const_iterator iter = m_histogramMap.find(name);

    if (m_histogramMap.end() == iter)
    {
        std::cout << "PandoraMonitoring::DrawHistogram, error: No histogram with name '"<< name <<"' exists." << std::endl;
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

void PandoraMonitoring::SaveAndCloseHistogram(const std::string &name, const std::string &fileName, const std::string &fileOptions)
{
    HistogramMap::iterator iter = m_histogramMap.find(name);

    if (m_histogramMap.end() == iter)
    {
        std::cout << "PandoraMonitoring::SaveHistogram, error: No histogram with name '" << name << "' exists." << std::endl;
        throw std::exception();
    }

    TFile* pTFile = new TFile(fileName.c_str(), fileOptions.c_str());

    iter->second->SetDirectory(pTFile);
    iter->second->Write(name.c_str(), TObject::kOverwrite);
    delete iter->second;

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
        std::cout << "PandoraMonitoring::DeleteHistogram, error: No histogram with name '"<< name <<"' exists." << std::endl;
        throw std::exception();
    }

    delete iter->second;
    m_histogramMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename VariableType>
void PandoraMonitoring::SetTreeVariable(const std::string &treeName, const std::string &variableName, VariableType variable )
{
    m_treeWrapper.Set( treeName, variableName, variable );
}

// instantiations of this template member function for the permitted types
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, float  variable );
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, int    variable );
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, double variable );

template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, std::vector<float>*  variable );
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, std::vector<int>*    variable );
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, std::vector<double>* variable );


//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::FillTree(const std::string &treeName)
{
    try {
        m_treeWrapper.Fill( treeName );
    }
    catch( TTreeWrapper::TreeNotFoundError& excpt )
    {
        std::cout << "PandoraMonitoring::FillTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
        throw;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::FillTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::PrintTree(const std::string &treeName)
{
    try {
        m_treeWrapper.Print( treeName );
    }
    catch( TTreeWrapper::TreeNotFoundError& excpt )
    {
        std::cout << "PandoraMonitoring::PrintTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
        throw;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::PrintTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
        throw;
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::SaveTree(const std::string &treeName, const std::string &fileName, const std::string &fileOptions)
{
    TTree* tree = NULL;
    try {
        tree = m_treeWrapper.GetTree( treeName );
    }
    catch( TTreeWrapper::TreeNotFoundError& excpt )
    {
        std::cout << "PandoraMonitoring::SaveTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
        throw;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::SaveTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
        throw;
    }

    TFile* pTFile = new TFile(fileName.c_str(), fileOptions.c_str());

    tree->SetDirectory(pTFile);
    tree->Write(treeName.c_str(), TObject::kOverwrite);

//    pTFile->Write();
    pTFile->Close();

    delete pTFile;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList)
{
    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    this->DrawDetectorOutline(detectorView);
    this->DrawTracks(detectorView, pTrackList);
    this->Pause();

    for(unsigned int i = 0; i < m_eventArrows.size(); ++i)
        delete m_eventArrows[i];

    for(unsigned int i = 0; i < m_eventMarkers.size(); ++i)
        delete m_eventMarkers[i];

    m_eventArrows.clear();
    m_eventMarkers.clear();
    delete pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawEvent(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList)
{
    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    this->DrawDetectorOutline(detectorView);
    this->DrawCaloHits(detectorView, pOrderedCaloHitList);
    this->Pause();

    for(unsigned int i = 0; i < m_eventMarkers.size(); ++i)
        delete m_eventMarkers[i];

    m_eventMarkers.clear();
    delete pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawEvent(DetectorView detectorView, const pandora::ClusterList *const pClusterList)
{
    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    this->DrawDetectorOutline(detectorView);
    this->DrawClusters(detectorView, pClusterList);
    this->Pause();

    for(unsigned int i = 0; i < m_eventMarkers.size(); ++i)
        delete m_eventMarkers[i];

    m_eventMarkers.clear();
    delete pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList, const pandora::OrderedCaloHitList *const pOrderedCaloHitList)
{
    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    this->DrawDetectorOutline(detectorView);
    this->DrawTracks(detectorView, pTrackList);
    this->DrawCaloHits(detectorView, pOrderedCaloHitList);
    this->Pause();

    for(unsigned int i = 0; i < m_eventArrows.size(); ++i)
        delete m_eventArrows[i];

    for(unsigned int i = 0; i < m_eventMarkers.size(); ++i)
        delete m_eventMarkers[i];

    m_eventArrows.clear();
    m_eventMarkers.clear();
    delete pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList, const pandora::ClusterList *const pClusterList)
{
    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    this->DrawDetectorOutline(detectorView);
    this->DrawTracks(detectorView, pTrackList);
    this->DrawClusters(detectorView, pClusterList);
    this->Pause();

    for(unsigned int i = 0; i < m_eventArrows.size(); ++i)
        delete m_eventArrows[i];

    for(unsigned int i = 0; i < m_eventMarkers.size(); ++i)
        delete m_eventMarkers[i];

    m_eventArrows.clear();
    m_eventMarkers.clear();
    delete pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DetectorOutlineTest(DetectorView detectorView)
{
    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    this->DrawDetectorOutline(detectorView);
    this->Pause();
    delete pCanvas;
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

void PandoraMonitoring::DrawDetectorOutline(DetectorView detectorView)
{
    if (!m_isOutlineConstructed)
        this->MakeDetectorOutline();

    switch(detectorView)
    {
    case DETECTOR_VIEW_XY:
        m_pXYAxes->Draw();
        for(unsigned int i = 0; i < m_2DObjectsXY.size(); ++i)
        {
            m_2DObjectsXY[i]->Draw("f");
            m_2DObjectsXY[i]->Draw("l");
        }
        break;

    case DETECTOR_VIEW_XZ:
        m_pXZAxes->Draw();
        for(unsigned int i = 0; i < m_2DObjectsXZ.size(); ++i)
        {
            m_2DObjectsXZ[i]->Draw();
        }
        break;

    default:
        std::cout << "PandoraMonitoring::DrawDetectorOutline, error: request for an unsupported detector view." << std::endl;
        throw std::exception();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawTracks(DetectorView detectorView, const pandora::TrackList *const pTrackList)
{
    for (pandora::TrackList::const_iterator iter = pTrackList->begin(), iterEnd = pTrackList->end(); iter != iterEnd; ++iter)
    {
        const EVENT::Track* pTrack = static_cast<const EVENT::Track*>((*iter)->GetParentTrackAddress());

        EVENT::TrackerHitVec trackHitVector = pTrack->getTrackerHits();
        const unsigned int nTrackHits(trackHitVector.size());

        float *pX = new float[nTrackHits];
        float *pY = new float[nTrackHits];
        float *pZ = new float[nTrackHits];

        for(unsigned int i = 0; i < nTrackHits; ++i)
        {
            EVENT::TrackerHit* pTrackHit = trackHitVector[i];
            pX[i] = (float)pTrackHit->getPosition()[0];
            pY[i] = (float)pTrackHit->getPosition()[1];
            pZ[i] = (float)pTrackHit->getPosition()[2];
        }

        TArrow *pTArrow = NULL;
        TPolyMarker *pTPolyMarker = NULL;

        const pandora::CartesianVector &trackSeedPosition((*iter)->GetTrackStateAtECal().GetPosition());
        const pandora::CartesianVector trackSeedDirection((*iter)->GetTrackStateAtECal().GetMomentum().GetUnitVector());

        if(DETECTOR_VIEW_XY == detectorView)
        {
            pTArrow = new TArrow(trackSeedPosition.GetX(), trackSeedPosition.GetY(), trackSeedPosition.GetX() + trackSeedDirection.GetX(),
                trackSeedPosition.GetY() + trackSeedDirection.GetY(), 0.01);
            pTPolyMarker = new TPolyMarker(nTrackHits, pX, pY);
        }
        else if(DETECTOR_VIEW_XZ == detectorView)
        {
            pTArrow = new TArrow(trackSeedPosition.GetZ(), trackSeedPosition.GetX(), trackSeedPosition.GetZ() + trackSeedDirection.GetZ(),
                trackSeedPosition.GetX() + trackSeedDirection.GetX(), 0.01);
            pTPolyMarker = new TPolyMarker(nTrackHits, pZ, pX);
        }

        pTPolyMarker->SetMarkerStyle(20);
        pTPolyMarker->SetMarkerColor(1);
        pTPolyMarker->SetMarkerSize(0.1);
        pTPolyMarker->Draw();

        pTArrow->Draw();

        m_eventArrows.push_back(pTArrow);
        m_eventMarkers.push_back(pTPolyMarker);

        delete [] pX;
        delete [] pY;
        delete [] pZ;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawCaloHits(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList)
{
    for (pandora::OrderedCaloHitList::const_iterator iter = pOrderedCaloHitList->begin(), iterEnd = pOrderedCaloHitList->end(); iter != iterEnd; ++iter)
    {
        const pandora::PseudoLayer pseudoLayer(iter->first);
        const pandora::CaloHitList *const pCaloHitList(iter->second);

        const unsigned int nCaloHits(pCaloHitList->size());

        if (0 == nCaloHits)
            continue;

        int i = 0;
        float *pX = new float[nCaloHits];
        float *pY = new float[nCaloHits];
        float *pZ = new float[nCaloHits];

        for (pandora::CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
        {
//            if ((DETECTOR_VIEW_XY == detectorView) && (pandora::ENDCAP == (*hitIter)->GetDetectorRegion()))
//                continue;

//            if ((DETECTOR_VIEW_XZ == detectorView) && (pandora::BARREL == (*hitIter)->GetDetectorRegion()))
//                continue;

            const pandora::CartesianVector &positionVector((*hitIter)->GetPositionVector());
            pX[i] = positionVector.GetX();
            pY[i] = positionVector.GetY();
            pZ[i] = positionVector.GetZ();
            i++;
        }

        TPolyMarker *pTPolyMarker = NULL;

        if(DETECTOR_VIEW_XY == detectorView)
        {
            pTPolyMarker = new TPolyMarker(i, pX, pY);
        }
        else if(DETECTOR_VIEW_XZ == detectorView)
        {
            pTPolyMarker = new TPolyMarker(i, pZ, pX);
        }

        pTPolyMarker->SetMarkerStyle(20);
        pTPolyMarker->SetMarkerSize(0.5);

        unsigned int layerModN(pseudoLayer % 7);
        unsigned int color = kBlack;

        switch (layerModN)
        {
            case 0: color = kRed;       break;
            case 1: color = kMagenta;   break;
            case 2: color = kBlue;      break;
            case 3: color = kCyan;      break;
            case 4: color = kGreen;     break;
            case 5: color = kYellow;    break;
            default:color = kBlack;     break;
        }

        pTPolyMarker->SetMarkerColor(color);
        pTPolyMarker->Draw();
        m_eventMarkers.push_back(pTPolyMarker);

        delete [] pX;
        delete [] pY;
        delete [] pZ;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawClusters(DetectorView detectorView, const pandora::ClusterList *const pClusterList)
{
    for (pandora::ClusterList::const_iterator clusterIter = pClusterList->begin(), clusterIterEnd = pClusterList->end();
        clusterIter != clusterIterEnd; ++clusterIter)
    {
        const pandora::OrderedCaloHitList orderedCaloHitList((*clusterIter)->GetOrderedCaloHitList());

        unsigned int hitIndex = 0;
        const unsigned int nHits((*clusterIter)->GetNCaloHits());

        float *pX = new float[nHits];
        float *pY = new float[nHits];
        float *pZ = new float[nHits];

        for (pandora::OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end();
            iter != iterEnd; ++iter)
        {
            for (pandora::CaloHitList::const_iterator caloHitIter = iter->second->begin(), caloHitIterEnd = iter->second->end();
                caloHitIter != caloHitIterEnd; ++caloHitIter)
            {
                const pandora::CartesianVector position = (*caloHitIter)->GetPositionVector();

                pX[hitIndex] = position.GetX();
                pY[hitIndex] = position.GetY();
                pZ[hitIndex] = position.GetZ();

                ++hitIndex;
            }
        }

        TPolyMarker *pTPolyMarker = NULL;

        if(DETECTOR_VIEW_XY == detectorView)
        {
            pTPolyMarker = new TPolyMarker(nHits, pX, pY);
        }
        else if(DETECTOR_VIEW_XZ == detectorView)
        {
            pTPolyMarker = new TPolyMarker(nHits, pZ, pX);
        }

        pTPolyMarker->SetMarkerStyle(20);

        if (1 == nHits)
        {
            pTPolyMarker->SetMarkerSize(0.2);
            pTPolyMarker->SetMarkerColor(3);
        }
        else
        {
            int colour(kBlue + static_cast<int>((*clusterIter)->GetElectromagneticEnergy()));
            pTPolyMarker->SetMarkerSize(0.4);
            pTPolyMarker->SetMarkerColor(colour);
        }

        pTPolyMarker->Draw();

        m_eventMarkers.push_back(pTPolyMarker);

        delete [] pX;
        delete [] pY;
        delete [] pZ;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeDetectorOutline()
{
    m_pXYAxes = new TH2F("xyview", "x-y view", 100, -5000, 5000, 100, -5000, 5000);
    m_pXYAxes->SetFillColor(kWhite);
    m_pXYAxes->SetStats(kFALSE);

    m_pXZAxes = new TH2F("xzview", "x-z view", 100, -5000, 5000, 100, -5000, 5000);
    m_pXZAxes->SetFillColor(kWhite);
    m_pXZAxes->SetStats(kFALSE);

    pandora::GeometryHelper *pGeometryHelper = pandora::GeometryHelper::GetInstance();

    typedef std::vector<pandora::GeometryHelper::SubDetectorParameters> SubDetectorParametersList;
    SubDetectorParametersList subDetectorParametersList;

    // Additional sub detectors are stored in a map from name to parameters. We don't care about the names here - just copy into a vector for sorting.
    const pandora::GeometryHelper::SubDetectorParametersMap &subDetectorParametersMap(pGeometryHelper->GetAdditionalSubDetectors());
    for (pandora::GeometryHelper::SubDetectorParametersMap::const_iterator iter = subDetectorParametersMap.begin(); iter != subDetectorParametersMap.end(); ++iter)
    {
        subDetectorParametersList.push_back(iter->second);
    }

    subDetectorParametersList.push_back(pGeometryHelper->GetECalBarrelParameters());
    subDetectorParametersList.push_back(pGeometryHelper->GetECalEndCapParameters());
    subDetectorParametersList.push_back(pGeometryHelper->GetHCalBarrelParameters());
    subDetectorParametersList.push_back(pGeometryHelper->GetHCalEndCapParameters());

    XYOutlineParametersList xyOutlineParametersList;
    xyOutlineParametersList.push_back(XYOutlineParameters(0, 0, pGeometryHelper->GetMainTrackerInnerRadius()));
    xyOutlineParametersList.push_back(XYOutlineParameters(0, 0, pGeometryHelper->GetMainTrackerOuterRadius()));
    this->MakeXZLayerOutline(pGeometryHelper->GetMainTrackerInnerRadius(), pGeometryHelper->GetMainTrackerOuterRadius(), 0,
        pGeometryHelper->GetMainTrackerZExtent());

    xyOutlineParametersList.push_back(XYOutlineParameters(0, 0, pGeometryHelper->GetCoilInnerRadius()));
    xyOutlineParametersList.push_back(XYOutlineParameters(0, 0, pGeometryHelper->GetCoilOuterRadius()));
    this->MakeXZLayerOutline(pGeometryHelper->GetCoilInnerRadius(), pGeometryHelper->GetCoilOuterRadius(), 0,
        pGeometryHelper->GetCoilZExtent());

    for (SubDetectorParametersList::const_iterator iter = subDetectorParametersList.begin(); iter != subDetectorParametersList.end(); ++iter)
    {
        XYOutlineParameters xyInnerOutlineParameters(iter->GetInnerSymmetryOrder(), iter->GetInnerPhiCoordinate(), iter->GetInnerRCoordinate());
        XYOutlineParameters xyOuterOutlineParameters(iter->GetOuterSymmetryOrder(), iter->GetOuterPhiCoordinate(), iter->GetOuterRCoordinate());
        xyOutlineParametersList.push_back(xyInnerOutlineParameters);
        xyOutlineParametersList.push_back(xyOuterOutlineParameters);
        this->MakeXZLayerOutline(iter->GetInnerRCoordinate(), iter->GetOuterRCoordinate(), iter->GetInnerZCoordinate(),
            iter->GetOuterZCoordinate());
    }

    std::sort(xyOutlineParametersList.begin(), xyOutlineParametersList.end(), XYOutlineParameters::Sort);

    for (XYOutlineParametersList::const_iterator iter = xyOutlineParametersList.begin(), iterEnd = xyOutlineParametersList.end();
        iter != iterEnd; ++iter)
    {
        this->MakeXYLayerOutline(*iter);
    }

    m_isOutlineConstructed = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeXYLayerOutline(const XYOutlineParameters &xyOutlineParameters)
{
    const int symmetryOrder(xyOutlineParameters.GetSymmetryOrder());
    const float closestDistanceToIp(xyOutlineParameters.GetClosestDistanceToIp());

    if (symmetryOrder > 2)
    {
        float *pXCoordinate = new float[symmetryOrder + 1];
        float *pYCoordinate = new float[symmetryOrder + 1];

        const float pi(3.1415927);
        const float x0(-1. * closestDistanceToIp * tan(pi / float(symmetryOrder)));
        const float y0(closestDistanceToIp);

        for (int i = 0; i <= symmetryOrder; ++i)
        {
            const float theta(xyOutlineParameters.GetPhi0() + (2 * pi * float(i) / float(symmetryOrder)));

            pXCoordinate[i] = x0 * cos(theta) + y0 * sin(theta);
            pYCoordinate[i] = y0 * cos(theta) - x0 * sin(theta);
        }

        TGraph *pTGraph = new TGraph(symmetryOrder + 1, pXCoordinate, pYCoordinate);
        pTGraph->SetFillStyle(1001);
        pTGraph->SetFillColor(kWhite);

        m_2DObjectsXY.push_back(pTGraph);

        delete [] pXCoordinate;
        delete [] pYCoordinate;
    }
    else
    {
        TArc* pTArc = new TArc(0., 0., closestDistanceToIp);
        pTArc->SetFillStyle(1001);
        pTArc->SetFillColor(kWhite);

        m_2DObjectsXY.push_back(pTArc);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeXZLayerOutline(float innerRCoordinate, float outerRCoordinate, float innerZCoordinate, float outerZCoordinate)
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

        TGraph *pTGraph = new TGraph(linesToDisplay, x, z);
        pTGraph->SetFillStyle(0);

        m_2DObjectsXZ.push_back(pTGraph);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool PandoraMonitoring::XYOutlineParameters::Sort(const XYOutlineParameters &lhs, const XYOutlineParameters &rhs)
{
    return (lhs.m_closestDistanceToIp > rhs.m_closestDistanceToIp);
}



} // namespace pandora_monitoring
