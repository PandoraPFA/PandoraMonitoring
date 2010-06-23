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
#include "Objects/MCParticle.h"

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
#include "TColor.h"


#define MINIMUM_ROOT_VERSION "5.24"
#ifdef ROOT_EVE
#if ( ROOT_VERSION_CODE >= ROOT_VERSION(5,24,0) )
#define USE_ROOT_EVE 1
#endif
#endif




#ifdef USE_ROOT_EVE
#include <TEveManager.h>
//#include <TEveProjectionManager.h>
#include <TEveEventManager.h>
#include <TEveViewer.h>
#include <TEvePointSet.h>
#include <TEveArrow.h>
#include <TEveRGBAPalette.h>

#include <TGeoXtru.h>
#include <TEveGeoShapeExtract.h>
#include <TEveGeoShape.h>

#include <TGeometry.h>
#include <TGeoMaterial.h>
#include <TGeoManager.h>
#include <TEveGeoNode.h>

#include "TEveTrackPropagator.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h" // for TEveRecTrack


#include <TSystem.h>
#include <TGeoManager.h>
#include <TGeoXtru.h>
#include <TGeoMatrix.h>
#include <TEveManager.h>
#include <TEveGeoNode.h>
#include <TGeoTube.h>
#include <TEveBoxSet.h>
#include <TGeoCompositeShape.h>
#endif

#include <vector>
#include <map>
#include <assert.h>
#include <math.h>
#include <iostream>


#include "PandoraMonitoring.h"

#include <cmath>
#include <fcntl.h>

namespace pandora_monitoring
{

bool  PandoraMonitoring::m_instanceFlag   = false;
bool  PandoraMonitoring::m_eveInitialized = false;
bool  PandoraMonitoring::m_openEveEvent   = false;
float PandoraMonitoring::m_scalingFactor = 0.1;

PandoraMonitoring* PandoraMonitoring::m_pPandoraMonitoring = NULL;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraMonitoring::~PandoraMonitoring()
{
    m_treeWrapper.Clear();
    #ifdef USE_ROOT_EVE
    if( m_eveInitialized )
    {
        TEveManager::Terminate();
        gSystem->ProcessEvents();
    }
    #endif
    m_pApplication->Terminate(0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DeleteInstance()
{
    delete GetInstance();
}

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraMonitoring *PandoraMonitoring::GetInstance()
{
    if(!m_instanceFlag)
    {
        m_pPandoraMonitoring = new PandoraMonitoring();
        m_instanceFlag = true;
        TColor::CreateColorWheel(); // create the ROOT color wheel
        gStyle->SetPalette(1);
        gStyle->SetNumberContours(99);
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
void PandoraMonitoring::SetTreeVariable(const std::string &treeName, const std::string &variableName, VariableType variable)
{
    m_treeWrapper.Set(treeName, variableName, variable);
}

// instantiations of this template member function for the permitted types
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, float  variable);
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, int    variable);
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, double variable);

template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, std::vector<float>*  variable);
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, std::vector<int>*    variable);
template void PandoraMonitoring::SetTreeVariable(const std::string&, const std::string&, std::vector<double>* variable);

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::FillTree(const std::string &treeName)
{
    try
    {
        m_treeWrapper.Fill(treeName);
    }
    catch(TTreeWrapper::TreeNotFoundError& excpt)
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
    try
    {
        m_treeWrapper.Print(treeName);
    }
    catch(TTreeWrapper::TreeNotFoundError& excpt)
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

void PandoraMonitoring::ScanTree(const std::string &treeName)
{
    try
    {
        m_treeWrapper.Scan(treeName);
    }
    catch(TTreeWrapper::TreeNotFoundError& excpt)
    {
        std::cout << "PandoraMonitoring::ScanTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
        throw;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::ScanTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::SaveTree(const std::string &treeName, const std::string &fileName, const std::string &fileOptions)
{
    try
    {
        TTree*& tree = m_treeWrapper.GetTree(treeName);

        TFile* pTFile = new TFile(fileName.c_str(), fileOptions.c_str());

        tree->SetDirectory(pTFile);
        tree->Write(treeName.c_str(), TObject::kOverwrite);

        pTFile->Close();

        tree = 0; // pointer not valid any more

        delete pTFile;
    }
    catch(TTreeWrapper::TreeNotFoundError& excpt)
    {
        std::cout << "PandoraMonitoring::SaveTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
        throw;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::SaveTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
        throw;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::AddClusterList(DetectorView detectorView, const pandora::ClusterList *const pClusterList, Color color)
{
    GetCanvas(detectorView);
    this->DrawClusters(detectorView, pClusterList, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::AddTrackList(DetectorView detectorView, const pandora::TrackList *const pTrackList, Color color)
{
    GetCanvas(detectorView);
    this->DrawTracks(detectorView, pTrackList, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::AddCaloHitList(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList, Color color)
{
    GetCanvas(detectorView);
    this->DrawCaloHits(detectorView, pOrderedCaloHitList, color);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TCanvas* PandoraMonitoring::GetCanvas(DetectorView detectorView)
{
    CanvasMap::iterator itCanvas = m_canvasMap.find(detectorView);
    if( itCanvas != m_canvasMap.end() ) // if canvas is already existing in the canvas-map
    {
        TCanvas* pCanvas = itCanvas->second;
        pCanvas->cd();
        return pCanvas;        // return the canvas
    }

    // create the canvas
    std::stringstream sstr;
    sstr << "PandoraMonitoring_";
    sstr << (detectorView==DETECTOR_VIEW_XY?"XY":"XZ");
    
    TCanvas *pCanvas = new TCanvas(sstr.str().c_str(), sstr.str().c_str(), 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    m_canvasMap.insert( std::make_pair(detectorView,pCanvas) );

    this->DrawDetectorOutline(detectorView);

    pCanvas->cd();
    return pCanvas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::ViewEvent()
{
    this->Pause();

    for(unsigned int i = 0; i < m_eventArrows.size(); ++i)
        delete m_eventArrows[i];

    for(unsigned int i = 0; i < m_eventMarkers.size(); ++i)
        delete m_eventMarkers[i];

    m_eventArrows.clear();
    m_eventMarkers.clear();

    for( CanvasMap::iterator itCanvas = m_canvasMap.begin(), itCanvasEnd = m_canvasMap.end(); itCanvas != itCanvasEnd; ++itCanvas )
    {
        TCanvas* pCanvas = itCanvas->second;
        delete pCanvas;
    }
    m_canvasMap.clear();
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

void PandoraMonitoring::DrawTracks(DetectorView detectorView, const pandora::TrackList *const pTrackList, Color color)
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
        pTPolyMarker->SetMarkerColor(GetColor(color));
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

void PandoraMonitoring::DrawCaloHits(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList, Color color)
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

        int col = GetColor(color);
        if( color == AUTO )
        {
            unsigned int layerModN(pseudoLayer % 7);

            switch (layerModN)
            {
                case 0: col = kRed;       break;
                case 1: col = kMagenta;   break;
                case 2: col = kBlue;      break;
                case 3: col = kCyan;      break;
                case 4: col = kGreen;     break;
                case 5: col = kYellow;    break;
                default:col = kBlack;     break;
            }
        }
        pTPolyMarker->SetMarkerColor(col);
        pTPolyMarker->Draw();
        m_eventMarkers.push_back(pTPolyMarker);

        delete [] pX;
        delete [] pY;
        delete [] pZ;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DrawClusters(DetectorView detectorView, const pandora::ClusterList *const pClusterList, Color color )
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
            pTPolyMarker->SetMarkerColor(GetColor(color));
        }
        else
        {
            pTPolyMarker->SetMarkerSize(0.4);
            if( color == AUTO )
            {
                int autoColor = kBlue + static_cast<int>((*clusterIter)->GetElectromagneticEnergy());
                pTPolyMarker->SetMarkerColor(autoColor);
            }
            else
                pTPolyMarker->SetMarkerColor(GetColor(color));
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

EColor PandoraMonitoring::GetColor(Color color)
{
    switch(color)
    {
    case WHITE:
        return kWhite;
    case BLACK:
        return kBlack;
    case RED:
        return kRed;
    case GREEN:
        return kGreen;
    case BLUE:
        return kBlue;
    case MAGENTA:
        return kMagenta;
    case CYAN:
        return kCyan;
    case VIOLET:
        return kViolet;
    case PINK:
        return kPink;
    case ORANGE:
        return kOrange;
    case YELLOW:
        return kYellow;
    case SPRING:
        return kSpring;
    case TEAL:
        return kTeal;
    case AZURE:
        return kAzure;
    case GRAY:
        return kGray;
    case DARKRED:
        return EColor(TColor::GetColorDark(kRed));
    case DARKGREEN:
        return EColor(TColor::GetColorDark(kGreen));
    case DARKBLUE:
        return EColor(TColor::GetColorDark(kBlue));
    case DARKMAGENTA:
        return EColor(TColor::GetColorDark(kMagenta));
    case DARKCYAN:
        return EColor(TColor::GetColorDark(kCyan));
    case DARKVIOLET:
        return EColor(TColor::GetColorDark(kViolet));
    case DARKPINK:
        return EColor(TColor::GetColorDark(kPink));
    case DARKORANGE:
        return EColor(TColor::GetColorDark(kOrange));
    case DARKYELLOW:
        return EColor(TColor::GetColorDark(kYellow));
    case LIGHTGREEN:
        return EColor(TColor::GetColorBright(kGreen));
    case LIGHTBLUE:
        return EColor(TColor::GetColorBright(kBlue));
    case LIGHTMAGENTA:
        return EColor(TColor::GetColorBright(kMagenta));
    case LIGHTCYAN:
        return EColor(TColor::GetColorBright(kCyan));
    case LIGHTVIOLET:
        return EColor(TColor::GetColorBright(kViolet));
    case LIGHTPINK:
        return EColor(TColor::GetColorBright(kPink));
    case LIGHTORANGE:
        return EColor(TColor::GetColorBright(kOrange));
    case LIGHTYELLOW:
        return EColor(TColor::GetColorBright(kYellow));
    default:
        return kBlack;
    }
    return kBlack;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::InitializeEve( Char_t transparency )
{
#ifndef USE_ROOT_EVE
    std::cout << "ERROR: Visualization with ROOT TEve needs ROOT version >= " << MINIMUM_ROOT_VERSION << " !" << std::endl;
#else
    if( PandoraMonitoring::m_eveInitialized )
    {
        if( !m_openEveEvent )
        {
            gEve->AddEvent( new TEveEventManager("Event","Event") );
            m_openEveEvent = true;
        }
        return;
    }

    gSystem->Load("libGeom");
    TGeoManager *geom = new TGeoManager("DetectorGeometry", "detector geometry");


    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);

    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Aluminium",2, matAl);


    //--- make the top container volume
    TGeoVolume *top = geom->MakeBox("Detector", Vacuum, 1000., 1000., 100.);

    geom->SetTopVolume(top);



    pandora::GeometryHelper *pGeometryHelper = pandora::GeometryHelper::GetInstance();

    typedef std::vector<std::pair<pandora::GeometryHelper::SubDetectorParameters, std::string> > SubDetectorParametersList;
    SubDetectorParametersList subDetectorParametersList;

    // Additional sub detectors are stored in a map from name to parameters. We don't care about the names here - just copy into a vector for sorting.
    const pandora::GeometryHelper::SubDetectorParametersMap &subDetectorParametersMap(pGeometryHelper->GetAdditionalSubDetectors());
    int i = 0;
    for (pandora::GeometryHelper::SubDetectorParametersMap::const_iterator iter = subDetectorParametersMap.begin(); iter != subDetectorParametersMap.end(); ++iter)
    { 
        std::stringstream sstr;
        sstr << "subDet_" << i;

        subDetectorParametersList.push_back(std::make_pair(iter->second, sstr.str()) );
        ++i;
    }

    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetECalBarrelParameters(), "ECalBarrel") );
    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetECalEndCapParameters(), "ECalEndCap") );
    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetHCalBarrelParameters(), "HCalBarrel") );
    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetHCalEndCapParameters(), "HCalEndCap") );

    TGeoVolume* mainTracker = NULL;
    mainTracker = MakePolygoneTube( "Tracker", 
                                    0, 0,
                                    pGeometryHelper->GetMainTrackerInnerRadius()*m_scalingFactor, 
                                    pGeometryHelper->GetMainTrackerOuterRadius()*m_scalingFactor, 
                                    0.0,0.0,
                                    pGeometryHelper->GetMainTrackerZExtent()*m_scalingFactor,
                                    Al );

    mainTracker->SetLineColor(kGreen);
    mainTracker->SetTransparency(transparency);
    mainTracker->SetVisibility(kFALSE);
    top->AddNode( mainTracker, 0, new TGeoTranslation(0,0,0));

    TGeoVolume* coil = NULL;
    coil = MakePolygoneTube( "Coil", 
                             0, 0,
                             pGeometryHelper->GetCoilInnerRadius()*m_scalingFactor, 
                             pGeometryHelper->GetCoilOuterRadius()*m_scalingFactor, 
                             0.0,0.0,
                             pGeometryHelper->GetCoilZExtent()*m_scalingFactor,
                             Al );

    coil->SetLineColor(kBlue);
    coil->SetTransparency(transparency);
    coil->SetVisibility(kFALSE);
    top->AddNode( coil, 0, new TGeoTranslation(0,0,0));
    


    Int_t col = 2;
    for (SubDetectorParametersList::const_iterator iter = subDetectorParametersList.begin(); iter != subDetectorParametersList.end(); ++iter)
    {
        bool left = true;
        for( int lr = 0; lr <= 1; ++lr )
        {
            const pandora::GeometryHelper::SubDetectorParameters& detPar = (*iter).first;
            std::string name = (*iter).second;


            std::stringstream sstr;
            sstr << name;
            sstr << (left? "_left" : "_right" );

            TGeoVolume* subDetVol = NULL;

            int sign = (left? -1 : 1 );
            double zMin = detPar.GetInnerZCoordinate()*m_scalingFactor;
            double zMax = detPar.GetOuterZCoordinate()*m_scalingFactor;
            double zThick = zMax-zMin;
            zMin *= sign;
            zMax *= sign;
            double zPosition = zMin+sign*(zThick/2.0);

            subDetVol = MakePolygoneTube( sstr.str().c_str(), 
                                          detPar.GetInnerSymmetryOrder(),
                                          detPar.GetOuterSymmetryOrder(),
                                          detPar.GetInnerRCoordinate()*m_scalingFactor,
                                          detPar.GetOuterRCoordinate()*m_scalingFactor,
                                          detPar.GetInnerPhiCoordinate(),
                                          detPar.GetOuterPhiCoordinate(),
                                          (zThick/2.0), 
                                          Al );

            subDetVol->SetLineColor(GetColor(Color(col)));
            subDetVol->SetFillColor(GetColor(Color(col)));
            subDetVol->SetTransparency(transparency);

            size_t found=name.find("subDet_");
            if (found!=std::string::npos)
                subDetVol->SetVisibility(kFALSE);

            top->AddNode( subDetVol, 0, new TGeoTranslation(0,0, zPosition) );
            left = false;
        }

        ++col;
    }




    //--- close the geometry
    geom->CloseGeometry();


    TEveManager::Create();

    TGeoNode* node = gGeoManager->GetTopNode();
    TEveGeoTopNode* en = new TEveGeoTopNode(gGeoManager, node);
    en->SetVisLevel(1);
    en->GetNode()->GetVolume()->SetVisibility(kFALSE);

    gEve->AddGlobalElement(en);

    gEve->GetDefaultViewer()->SetMainColor( kWhite );
    //    gEve->SetMainColor( kWhite );
    gEve->Redraw3D(kTRUE);

    m_eveInitialized = true;
    m_openEveEvent = true;
#endif
}



//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::ComputePolygonCorners( int symmetryOrder, double closestDistanceToIp, double phi0, std::vector<std::pair<Double_t,Double_t> > &coordinates )
{
    if (symmetryOrder > 2)
    {
        const Double_t pi(3.1415927);
        const Double_t x0(-1. * closestDistanceToIp * tan(pi / Double_t(symmetryOrder)));
        const Double_t y0(closestDistanceToIp);

        for (int i = 0; i < symmetryOrder; ++i)
        {
            Double_t theta = 0.f; 
            theta = phi0 + (2 * pi * Double_t(i) / Double_t(symmetryOrder));

            Double_t x = x0 * cos(theta) + y0 * sin(theta);
            Double_t y = y0 * cos(theta) - x0 * sin(theta);
            coordinates.push_back( std::pair<double,double>( x, y ) );
        }
    }
}



//------------------------------------------------------------------------------------------------------------------------------------------

TGeoVolume* PandoraMonitoring::MakePolygoneTube( std::string name, 
                                                 int innerSymmetryOrder, int outerSymmetryOrder, 
                                                 double innerClosestDistanceToIp, double outerClosestDistanceToIp, 
                                                 double innerPhi0, double outerPhi0, 
                                                 double halfLength, 
                                                 TGeoMedium* medium )
{
#ifndef USE_ROOT_EVE
    std::cout << "ERROR: Visualization with ROOT TEve needs ROOT version >= " << MINIMUM_ROOT_VERSION << " !" << std::endl;
    return NULL;
#else
    TGeoShape* inner = MakePolygoneTube( innerSymmetryOrder, innerClosestDistanceToIp, innerPhi0, halfLength+2 );
    TGeoShape* outer = MakePolygoneTube( outerSymmetryOrder, outerClosestDistanceToIp, outerPhi0, halfLength );

    std::string nameInner = (name+"_I");
    std::string nameOuter = (name+"_O");

    inner->SetName( nameInner.c_str() );
    outer->SetName( nameOuter.c_str() );

    std::string formula = nameOuter+"-"+nameInner;

    TGeoCompositeShape* cs = new TGeoCompositeShape( (name+"_shape").c_str(), formula.c_str() );
    TGeoVolume* csVol = new TGeoVolume( name.c_str(), cs, medium );

    return csVol;
#endif
}


//------------------------------------------------------------------------------------------------------------------------------------------

TGeoShape* PandoraMonitoring::MakePolygoneTube( int symmetryOrder, double closestDistanceToIp, double phi, double halfLength )
{
#ifndef USE_ROOT_EVE
    std::cout << "ERROR: Visualization with ROOT TEve needs ROOT version >= " << MINIMUM_ROOT_VERSION << " !" << std::endl;
    return NULL;
#else
    if( symmetryOrder <= 2 )
    {
        TGeoShape* tube     = new TGeoTube(0, closestDistanceToIp, halfLength );
        return tube;
    }

    std::vector<std::pair<Double_t,Double_t> > vertices;

    ComputePolygonCorners( symmetryOrder, closestDistanceToIp, phi, vertices );
   
    const Int_t nvertices = vertices.size();
   
    Double_t* x = new Double_t[nvertices];
    Double_t* y = new Double_t[nvertices];

    int index = 0;
    for( std::vector< std::pair<Double_t,Double_t> >::iterator itCoord = vertices.begin(), itCoordEnd = vertices.end(); itCoord != itCoordEnd; ++itCoord )
    {
        x[index] = (*itCoord).first;
        y[index] = (*itCoord).second;

        ++index;
    }

    TGeoXtru *xtru = new TGeoXtru(2);

    xtru->DefinePolygon(nvertices,x,y);
    Double_t z0 = -halfLength, x0 = 0, y0 = 0;
    Double_t z1 = halfLength,   x1 = 0, y1 = 0;

    Double_t scale0 = 1.0;
    xtru->DefineSection(0, z0, x0, y0, scale0); // Z position, offset and scale for first section
    xtru->DefineSection(1, z1, x1, y1, scale0); // -''- go forward

    delete[] x;
    delete[] y;
    
    return xtru;
#endif
}




//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::View()
{
#ifndef USE_ROOT_EVE
    std::cout << "ERROR: Visualization with ROOT TEve needs ROOT version >= " << MINIMUM_ROOT_VERSION << " !" << std::endl;
#else
    InitializeEve();

    gEve->Redraw3D();
    
    this->Pause();

    TEveEventManager* current = gEve->GetCurrentEvent();
    if( current )
        current->SetRnrSelfChildren(kFALSE,kFALSE);

    m_openEveEvent = false;

    std::cout << "View done" << std::endl;
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement* PandoraMonitoring::VisualizeClusters(const pandora::ClusterList *const pClusterList, std::string& nameInput, TEveElement* parent, Color color )
{
#ifndef USE_ROOT_EVE
    std::cout << "ERROR: Visualization with ROOT TEve needs ROOT version >= " << MINIMUM_ROOT_VERSION << " !" << std::endl;
    return NULL;
#else

    InitializeEve();

    std::string name = nameInput;
    if( name.empty() )
    {
        std::stringstream sstr;
        sstr << "clusters, #" << pClusterList->size();
        name = sstr.str();
    }
    
    TEveElement* clusterListElement = new TEveElementList();
    clusterListElement->SetElementNameTitle( name.c_str(), name.c_str() );



    for (pandora::ClusterList::const_iterator clusterIter = pClusterList->begin(), clusterIterEnd = pClusterList->end();
        clusterIter != clusterIterEnd; ++clusterIter)
    { 
        pandora::Cluster* pCluster = (*clusterIter);
        if( pCluster->GetNCaloHits() == 0 )
            continue;

        Color clusterColor = color;
        if( color == AUTO )
        {
            const pandora::TrackList& pTrackList = pCluster->GetAssociatedTrackList();
            bool clusterHasTracks = !(pTrackList.empty());
            bool clusterIsPhoton  = pCluster->IsPhoton();
            
            if( clusterIsPhoton )
                clusterColor = DARKYELLOW;
            else if( clusterHasTracks )
                clusterColor = MAGENTA;
            else
                clusterColor = LIGHTBLUE;
        }

        std::stringstream sstr;
        sstr << "Eem=" << pCluster->GetElectromagneticEnergy() << " Ehad=" << pCluster->GetHadronicEnergy() << " NHits=" << pCluster->GetNCaloHits();

//         TEveElement* clusterElement = new TEveElementList();
//         clusterElement->SetElementNameTitle( sstr.str().c_str(), sstr.str().c_str() );

        // show hits
        const pandora::OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        TEveElement* caloHitsElement = VisualizeCaloHits(&orderedCaloHitList, sstr.str().c_str(), clusterListElement, clusterColor );



        try
        {
            const pandora::ClusterHelper::ClusterFitResult& fit = pCluster->GetFitToAllHitsResult();
            const pandora::CartesianVector intercept = fit.GetIntercept()*m_scalingFactor;
            pandora::CartesianVector direction       = fit.GetDirection()*m_scalingFactor;

            const double length = 100;
            pandora::CartesianVector displacedStart = intercept- (direction*(length/2));
            direction *= length;

            TEveArrow* clusterArrow = new TEveArrow(direction.GetX(), direction.GetY(), direction.GetZ(), displacedStart.GetX(),displacedStart.GetY(), displacedStart.GetZ());
            clusterArrow->SetConeR( 0.03 );
            clusterArrow->SetConeL( 0.2 );
            clusterArrow->SetMainColor( GetColor(clusterColor) );
            clusterArrow->SetPickable( kTRUE );
            clusterArrow->SetElementNameTitle( sstr.str().c_str(), sstr.str().c_str() );
            
            caloHitsElement->AddElement( clusterArrow );
        }
        catch(...)
        {
        }


        // show tracks
        const pandora::TrackList& pTrackList = pCluster->GetAssociatedTrackList();
        if( !pTrackList.empty() )
        {
            VisualizeTracks(&pTrackList, "", caloHitsElement, color );
        }

    }

    if( parent )
    {
        parent->AddElement(clusterListElement);
    }
    else
    {
        gEve->AddElement(clusterListElement);
        gEve->Redraw3D();
    }

    return clusterListElement;
#endif
}


//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement* PandoraMonitoring::VisualizeCaloHits(const pandora::OrderedCaloHitList *const pOrderedCaloHitList, std::string name, TEveElement* parent, Color color )
{
#ifndef USE_ROOT_EVE
    std::cout << "ERROR: Visualization with ROOT TEve needs ROOT version >= " << MINIMUM_ROOT_VERSION << " !" << std::endl;
    return NULL;
#else
    InitializeEve();

    if( name.empty() )
        name = "Hits";

    TEveBoxSet* hits = new TEveBoxSet(name.c_str());
    hits->Reset(TEveBoxSet::kBT_FreeBox, kTRUE, 64);
//    hits->SetMainColor( GetColor( color ) );
    //    hits->UseSingleColor();

    TEvePointSet* hitsMarkers = new TEvePointSet((name+"_markers").c_str());
    hits->AddElement(hitsMarkers);

    hits->SetOwnIds(kTRUE);

    Float_t corners[24]; // 8 corners x 3 dimensions

    for (pandora::OrderedCaloHitList::const_iterator iter = pOrderedCaloHitList->begin(), iterEnd = pOrderedCaloHitList->end();
         iter != iterEnd; ++iter)
    {
//        uint layer = iter->first;
        pandora::CaloHitList* pCaloHitList = iter->second;

        uint hitIndex = 0;
        for (pandora::CaloHitList::const_iterator caloHitIter = pCaloHitList->begin(), caloHitIterEnd = pCaloHitList->end();
             caloHitIter != caloHitIterEnd; ++caloHitIter)
        {
            const pandora::CaloHit* caloHit = (*caloHitIter);

            std::stringstream sstr;

            const pandora::CartesianVector position = caloHit->GetPositionVector()*m_scalingFactor;
            const float sizeU = caloHit->GetCellSizeU()*m_scalingFactor;
            const float sizeV = caloHit->GetCellSizeV()*m_scalingFactor;
            const float thickness = caloHit->GetCellThickness()*m_scalingFactor;

            float hitEnergyEm = caloHit->GetElectromagneticEnergy();
            float hitEnergyHad = caloHit->GetHadronicEnergy();

            sstr << "Eem=" << hitEnergyEm << " Eh=" << hitEnergyHad;

            const pandora::MCParticle* pMcParticle = NULL;
            caloHit->GetMCParticle(pMcParticle);

            if( pMcParticle )
                sstr << " MC-PID=" << pMcParticle->GetParticleId();

            // direction vectors (normal and U)
            const pandora::CartesianVector normal   = caloHit->GetNormalVector();

            const pandora::DetectorRegion& detReg = caloHit->GetDetectorRegion();
            pandora::CartesianVector dirU( 0, 0, 1 );
            if( detReg == pandora::ENDCAP )
                dirU.SetValues( 0, 1, 0 );

            // compute the corners of the calohit calorimeter-cell
            MakeCaloHitCell( position, normal, dirU, sizeU, sizeV, thickness, corners );


            // at hit marker
            hitsMarkers->SetNextPoint( position.GetX(), position.GetY(), position.GetZ() );
            //            hitMarkers->SetPointId(new TNamed(Form("Point %d", i), ""));
            hitsMarkers->SetMarkerColor(GetColor(color));
            hitsMarkers->SetMarkerSize(5);
            hitsMarkers->SetMarkerStyle(4);

            // add calorimeter-cell for calo-hit
            hits->AddBox( corners );
            hits->SetPickable( kTRUE );
            
//            int transparency = int(255-hitEnergyEm*1000);
//            hits->DigitColor(GetColor(color),transparency);
            hits->DigitColor(GetColor(color));

            hits->SetElementTitle(sstr.str().c_str());
            //            hits->DigitValue( layer );

//             float userData[1];
//             userData[0] = hitEnergy;
//             hits->DigitUserData( userData );
            ++hitIndex;
        }
    }
    
    if( parent )
    {
        parent->AddElement(hits);
    }
    else
    {
        gEve->AddElement(hits);
        gEve->Redraw3D();
    }

    return hits;
#endif
}


//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeCaloHitCell( const pandora::CartesianVector& position, 
                                         const pandora::CartesianVector& normal, 
                                         const pandora::CartesianVector& directionU, 
                                         const float cellSizeU, const float cellSizeV, const float cellSizeThickness,
                                         Float_t corners[24] )
{
    pandora::CartesianVector norm   = normal.GetUnitVector();
    pandora::CartesianVector dirU   = directionU.GetUnitVector();
    const float u2 = cellSizeU/2.0 ;
    const float v2 = cellSizeV/2.0 ;
    const float t2 = cellSizeThickness/2.0 ;

    pandora::CartesianVector dirV( norm.GetCrossProduct( dirU ));
    const float magnitudeV = dirV.GetMagnitude();
    if( magnitudeV<0.00001 )
    {
        std::cout << "ERROR/MakeCaloHitCell: vector of direction U is equal to normal vector!" << std::endl;
        throw std::exception();
    }
    else
    {
        dirV *= 1/magnitudeV;
    }
        
    pandora::CartesianVector cornerVector[8];

    dirU *= u2;
    dirV *= v2;
    norm *= t2;

    // compute all 8 corners for the box
    cornerVector[0] = position - dirU - dirV -norm;
    cornerVector[1] = position + dirU - dirV -norm;
    cornerVector[2] = position + dirU + dirV -norm;
    cornerVector[3] = position - dirU + dirV -norm;

    cornerVector[4] = position - dirU - dirV +norm;
    cornerVector[5] = position + dirU - dirV +norm;
    cornerVector[6] = position + dirU + dirV +norm;
    cornerVector[7] = position - dirU + dirV +norm;


    //    corners = new Float_t[24];
    corners[0] = cornerVector[0].GetX();
    corners[1] = cornerVector[0].GetY();
    corners[2] = cornerVector[0].GetZ();

    corners[3] = cornerVector[1].GetX();
    corners[4] = cornerVector[1].GetY();
    corners[5] = cornerVector[1].GetZ();

    corners[6] = cornerVector[2].GetX();
    corners[7] = cornerVector[2].GetY();
    corners[8] = cornerVector[2].GetZ();

    corners[9] = cornerVector[3].GetX();
    corners[10] = cornerVector[3].GetY();
    corners[11] = cornerVector[3].GetZ();

    corners[12] = cornerVector[4].GetX();
    corners[13] = cornerVector[4].GetY();
    corners[14] = cornerVector[4].GetZ();

    corners[15] = cornerVector[5].GetX();
    corners[16] = cornerVector[5].GetY();
    corners[17] = cornerVector[5].GetZ();

    corners[18] = cornerVector[6].GetX();
    corners[19] = cornerVector[6].GetY();
    corners[20] = cornerVector[6].GetZ();

    corners[21] = cornerVector[7].GetX();
    corners[22] = cornerVector[7].GetY();
    corners[23] = cornerVector[7].GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement* PandoraMonitoring::VisualizeTracks(const pandora::TrackList *const pTrackList, std::string name, TEveElement* parent, Color color )
{
#ifndef USE_ROOT_EVE
    std::cout << "ERROR: Visualization with ROOT TEve needs ROOT version >= " << MINIMUM_ROOT_VERSION << " !" << std::endl;
    return NULL;
#else

    InitializeEve();

    TEveTrackList *trackList = new TEveTrackList();
    trackList->SetMainColor( GetColor( TEAL ) );

    if( name.empty() )
        name = "Tracks";
    TEveTrackPropagator* propagator = trackList->GetPropagator();

    bool isRungeKutta = false; // if false --> helix fit
    if (isRungeKutta)
        propagator->SetStepper(TEveTrackPropagator::kRungeKutta);

    trackList->SetElementNameTitle( name.c_str(), name.c_str() );


    pandora::GeometryHelper *pGeometryHelper = pandora::GeometryHelper::GetInstance();

    float magneticField = pGeometryHelper->GetBField();
    propagator->SetMagFieldObj(new TEveMagFieldConst(0., 0., magneticField));
//     propagator->SetMaxR( pGeometryHelper->GetMainTrackerOuterRadius()*m_scalingFactor );
//     propagator->SetMaxZ( pGeometryHelper->GetMainTrackerZExtent()*m_scalingFactor );
    propagator->SetMaxR( pGeometryHelper->GetECalBarrelParameters().GetOuterRCoordinate()*m_scalingFactor );
    propagator->SetMaxZ( pGeometryHelper->GetECalEndCapParameters().GetOuterZCoordinate()*m_scalingFactor );
    propagator->SetMaxOrbs( 5 );

    for (pandora::TrackList::const_iterator trackIter = pTrackList->begin(), trackIterEnd = pTrackList->end();
        trackIter != trackIterEnd; ++trackIter)
    { 
        pandora::Track* pandoraTrack = (*trackIter);

        const pandora::TrackState& trackState = pandoraTrack->GetTrackStateAtStart();
        const pandora::CartesianVector position = trackState.GetPosition()*m_scalingFactor;
        const pandora::CartesianVector& momentum = trackState.GetMomentum();

        const pandora::TrackState& trackStateAtECal = pandoraTrack->GetTrackStateAtECal();
        const pandora::CartesianVector positionAtECal = trackStateAtECal.GetPosition()*m_scalingFactor;

        const pandora::TrackState& trackStateAtEnd = pandoraTrack->GetTrackStateAtEnd();
        const pandora::CartesianVector positionAtEnd = trackStateAtEnd.GetPosition()*m_scalingFactor;
        const pandora::CartesianVector& momentumAtEnd = trackStateAtEnd.GetMomentum();

        int chargeSign = pandoraTrack->GetChargeSign();

        Color trackColor = color;
        if( color == AUTO )
        {
            if( chargeSign > 0 )
                trackColor = RED;
            else if( chargeSign < 0 )
                trackColor = GREEN;
            else // should not happen
                trackColor = AZURE;
            
        }
        

        TEveRecTrack *rc = new TEveRecTrack();
        rc->fV.Set(position.GetX(),position.GetY(),position.GetZ());
        rc->fP.Set(momentum.GetX(),momentum.GetY(),momentum.GetZ());
        rc->fSign = -chargeSign; // "-" because of strange convention in ALICE : see http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=9456&p=40325&hilit=teve+histogram#p40325

        TEveTrack* track = new TEveTrack(rc, propagator);
        track->SetName(Form("Charge %d, PID %d, p %f", chargeSign, pandoraTrack->GetParticleId(), momentum.GetMagnitude()));
        
        TEvePathMark* pmEnd = new TEvePathMark(TEvePathMark::kReference);
        pmEnd->fV.Set(positionAtEnd.GetX(),positionAtEnd.GetY(),positionAtEnd.GetZ());
        pmEnd->fP.Set(momentumAtEnd.GetX(),momentumAtEnd.GetY(),momentumAtEnd.GetZ());
        track->AddPathMark(*pmEnd);

        TEvePathMark* pmECal = new TEvePathMark(TEvePathMark::kDecay);
        pmECal->fV.Set(positionAtECal.GetX(),positionAtECal.GetY(),positionAtECal.GetZ());
        track->AddPathMark(*pmECal);

        track->SetLineColor( GetColor(trackColor) );
        track->SetLineWidth( 1 );

        track->SetPickable( kTRUE );

        trackList->AddElement( track );
        track->MakeTrack();
    }

    if( parent )
        parent->AddElement(trackList);
    else
    {
        gEve->AddElement(trackList);
        gEve->Redraw3D();
    }

    return trackList;
#endif
}




//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool PandoraMonitoring::XYOutlineParameters::Sort(const XYOutlineParameters &lhs, const XYOutlineParameters &rhs)
{
    return (lhs.m_closestDistanceToIp > rhs.m_closestDistanceToIp);
}



} // namespace pandora_monitoring
