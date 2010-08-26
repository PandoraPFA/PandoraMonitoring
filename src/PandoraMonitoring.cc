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
#include "Objects/ParticleFlowObject.h"
#include "Objects/OrderedCaloHitList.h"
#include "Objects/Track.h"
#include "Objects/MCParticle.h"

// ROOT include files
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TTree.h"

#include <TEveManager.h>
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

#include "PandoraMonitoring.h"

#include <vector>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <limits>

namespace pandora_monitoring
{

bool  PandoraMonitoring::m_instanceFlag = false;
bool  PandoraMonitoring::m_eveInitialized = false;
bool  PandoraMonitoring::m_openEveEvent = false;
float PandoraMonitoring::m_scalingFactor = 0.1;
int   PandoraMonitoring::m_eventDisplayCounter = 0;

PandoraMonitoring* PandoraMonitoring::m_pPandoraMonitoring = NULL;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraMonitoring::~PandoraMonitoring()
{
    m_treeWrapper.Clear();

    if (m_eveInitialized)
    {
        TEveManager::Terminate();
        gSystem->ProcessEvents();
    }

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
    if (!m_instanceFlag)
    {
        m_pPandoraMonitoring = new PandoraMonitoring();
        m_instanceFlag = true;
        TColor::CreateColorWheel();
        gStyle->SetPalette(1);
        gStyle->SetNumberContours(99);
    }

    return m_pPandoraMonitoring;
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

void PandoraMonitoring::InitializeEve(Char_t transparency)
{
    std::stringstream sstr;
    sstr << "Event Display " << m_eventDisplayCounter;

    if( PandoraMonitoring::m_eveInitialized )
    {
        TEveEventManager* currEv = gEve->GetCurrentEvent();
        if( currEv )
            currEv->SetElementNameTitle(sstr.str().c_str(),sstr.str().c_str());

        if( !m_openEveEvent )
        {
            gEve->AddEvent( new TEveEventManager(sstr.str().c_str(),sstr.str().c_str()) );
            m_openEveEvent = true;
            m_eventDisplayCounter++;
        }
        return;
    }

    gSystem->Load("libGeom");
    TGeoManager *pGeoManager = new TGeoManager("DetectorGeometry", "detector geometry");

    //--- define some materials
    TGeoMaterial *pVacuumMaterial = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *pAluminiumMaterial = new TGeoMaterial("Aluminium", 26.98, 13, 2.7);

    //--- define some media
    TGeoMedium *pVacuum = new TGeoMedium("Vacuum",1, pVacuumMaterial);
    TGeoMedium *pAluminium = new TGeoMedium("Aluminium",2, pAluminiumMaterial);

    //--- make the top container volume
    TGeoVolume *pMainDetectorVolume = pGeoManager->MakeBox("Detector", pVacuum, 1000., 1000., 100.);
    pGeoManager->SetTopVolume(pMainDetectorVolume);

    this->InitializeSubDetectors(pMainDetectorVolume, pAluminium, transparency);
    this->InitializeGaps(pMainDetectorVolume, pVacuum, transparency);

    //--- close the geometry
    pGeoManager->CloseGeometry();

    TEveManager::Create();

    TGeoNode* pGeoNode = gGeoManager->GetTopNode();
    TEveGeoTopNode* pEveGeoTopNode = new TEveGeoTopNode(gGeoManager, pGeoNode);
    pEveGeoTopNode->SetVisLevel(1);
    pEveGeoTopNode->GetNode()->GetVolume()->SetVisibility(kFALSE);

    gEve->AddGlobalElement(pEveGeoTopNode);

    gEve->GetDefaultViewer()->SetMainColor(kWhite);
    gEve->Redraw3D(kTRUE);

    m_eveInitialized = true;
    m_openEveEvent = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::InitializeSubDetectors(TGeoVolume *pMainDetectorVolume, TGeoMedium *pSubDetectorMedium, Char_t transparency)
{
    pandora::GeometryHelper *pGeometryHelper = pandora::GeometryHelper::GetInstance();

    typedef std::vector<std::pair<pandora::GeometryHelper::SubDetectorParameters, std::string> > SubDetectorParametersList;
    SubDetectorParametersList subDetectorParametersList;

    const pandora::GeometryHelper::SubDetectorParametersMap &subDetectorParametersMap(pGeometryHelper->GetAdditionalSubDetectors());

    for (pandora::GeometryHelper::SubDetectorParametersMap::const_iterator iter = subDetectorParametersMap.begin(); iter != subDetectorParametersMap.end(); ++iter)
    {
        subDetectorParametersList.push_back(std::make_pair(iter->second, iter->first));
    }

    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetECalBarrelParameters(), "ECalBarrel") );
    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetECalEndCapParameters(), "ECalEndCap") );
    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetHCalBarrelParameters(), "HCalBarrel") );
    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetHCalEndCapParameters(), "HCalEndCap") );
    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetMuonBarrelParameters(), "MuonBarrel") );
    subDetectorParametersList.push_back(std::make_pair(pGeometryHelper->GetMuonEndCapParameters(), "MuonEndCap") );

    typedef std::set<std::string> StringSet;
    StringSet setInvisible;
    setInvisible.insert("MuonBarrel");
    setInvisible.insert("MuonEndCap");

    TGeoVolume* mainTracker = NULL;
    mainTracker = MakePolygonTube("Tracker", 0, 0,pGeometryHelper->GetMainTrackerInnerRadius() * m_scalingFactor,
        pGeometryHelper->GetMainTrackerOuterRadius()*m_scalingFactor, 0.,0., pGeometryHelper->GetMainTrackerZExtent() * m_scalingFactor, pSubDetectorMedium);

    mainTracker->SetLineColor(kGreen);
    mainTracker->SetTransparency(transparency);
    mainTracker->SetVisibility(kFALSE);
    pMainDetectorVolume->AddNode(mainTracker, 0, new TGeoTranslation(0, 0, 0));

    TGeoVolume* coil = NULL;
    coil = MakePolygonTube("Coil", 0, 0, pGeometryHelper->GetCoilInnerRadius() * m_scalingFactor,
        pGeometryHelper->GetCoilOuterRadius() * m_scalingFactor, 0.,0., pGeometryHelper->GetCoilZExtent() * m_scalingFactor, pSubDetectorMedium);

    coil->SetLineColor(kBlue);
    coil->SetTransparency(transparency);
    coil->SetVisibility(kFALSE);
    pMainDetectorVolume->AddNode(coil, 0, new TGeoTranslation(0,0,0));

    Int_t col = 2;
    for (SubDetectorParametersList::const_iterator iter = subDetectorParametersList.begin(); iter != subDetectorParametersList.end(); ++iter)
    {
        bool left = true;
        for( int lr = 0; lr <= 1; ++lr )
        {
            const pandora::GeometryHelper::SubDetectorParameters& detPar = (*iter).first;
            std::string name = (*iter).second;

            StringSet::iterator itSetInvisible = setInvisible.find(name);
            bool drawInvisible = (itSetInvisible != setInvisible.end() ? true : false );

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

            subDetVol = MakePolygonTube(sstr.str().c_str(), detPar.GetInnerSymmetryOrder(), detPar.GetOuterSymmetryOrder(),
                detPar.GetInnerRCoordinate() * m_scalingFactor, detPar.GetOuterRCoordinate() * m_scalingFactor,
                detPar.GetInnerPhiCoordinate(), detPar.GetOuterPhiCoordinate(), (zThick / 2.), pSubDetectorMedium);

            subDetVol->SetLineColor(GetROOTColor(Color(col)));
            subDetVol->SetFillColor(GetROOTColor(Color(col)));
            subDetVol->SetTransparency(transparency);

            size_t found=name.find("subDet_");
            if (found!=std::string::npos || drawInvisible)
                subDetVol->SetVisibility(kFALSE);

            pMainDetectorVolume->AddNode(subDetVol, 0, new TGeoTranslation(0, 0, zPosition));
            left = false;
        }
        ++col;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::InitializeGaps(TGeoVolume *pMainDetectorVolume, TGeoMedium *pGapMedium, Char_t transparency)
{
    pandora::GeometryHelper *pGeometryHelper = pandora::GeometryHelper::GetInstance();

    const pandora::GeometryHelper::GapList &gapList(pGeometryHelper->GetGapList());
    unsigned int gapCounter(0);

    for (pandora::GeometryHelper::GapList::const_iterator iter = gapList.begin(), iterEnd = gapList.end(); iter != iterEnd; ++iter)
    {
        std::string gapName("gap" + pandora::TypeToString(gapCounter++));

        pandora::GeometryHelper::BoxGap *pBoxGap = NULL;
        pBoxGap = dynamic_cast<pandora::GeometryHelper::BoxGap *>(*iter);

        if (NULL != pBoxGap)
        {
            TGeoShape *pGapShape = new TGeoBBox(gapName.c_str(), 0.5f * pBoxGap->m_side1.GetMagnitude() * m_scalingFactor,
                0.5f * pBoxGap->m_side2.GetMagnitude() * m_scalingFactor, 0.5f * pBoxGap->m_side3.GetMagnitude() * m_scalingFactor);

            TGeoVolume *pGapVol = new TGeoVolume(gapName.c_str(), pGapShape, pGapMedium);

            const float vertexZ(pBoxGap->m_vertex.GetZ());
            static const float hcalEndCapInnerZ(pGeometryHelper->GetHCalEndCapParameters().GetInnerZCoordinate());

            // TODO Remove ILD-specific correction, required for endcap box gaps that do not point back to origin in xy plane.
            //      Pandora gaps are self-describing (four vectors), but this does not map cleanly to TGeoBBox class.
            //      Best solution may be to move to different root TGeoShape.
            static const float pi(std::acos(-1.));
            const float correction((std::fabs(vertexZ) < hcalEndCapInnerZ) ? 0 : ((vertexZ > 0) ? pi / 4.f : -pi / 4.f));
            const float phi(correction + std::atan2(pBoxGap->m_vertex.GetX(), pBoxGap->m_vertex.GetY()));

            const TGeoTranslation trans("trans",
                ( 0.5f * pBoxGap->m_side1.GetMagnitude() * std::cos(phi) + 0.5f * pBoxGap->m_side2.GetMagnitude() * std::sin(phi) + pBoxGap->m_vertex.GetX()) * m_scalingFactor,
                (-0.5f * pBoxGap->m_side1.GetMagnitude() * std::sin(phi) + 0.5f * pBoxGap->m_side2.GetMagnitude() * std::cos(phi) + pBoxGap->m_vertex.GetY()) * m_scalingFactor,
                ( 0.5f * (2.f * pBoxGap->m_vertex.GetZ() + pBoxGap->m_side3.GetZ()) * m_scalingFactor));

            const TGeoRotation rot("rot", -180.f * phi / pi, 0, 0);

            pGapVol->SetLineColor(1);
            pGapVol->SetFillColor(1);
            pGapVol->SetTransparency(transparency + 23);

            pMainDetectorVolume->AddNode(pGapVol, 0, new TGeoCombiTrans(trans, rot));
            continue;
        }

        pandora::GeometryHelper::ConcentricGap *pConcentricGap = NULL;
        pConcentricGap = dynamic_cast<pandora::GeometryHelper::ConcentricGap *>(*iter);

        if (NULL != pConcentricGap)
        {
            const double zMin = pConcentricGap->m_minZCoordinate * m_scalingFactor;
            const double zMax = pConcentricGap->m_maxZCoordinate * m_scalingFactor;
            const double zThick = zMax - zMin;

            TGeoVolume *pGapVol = MakePolygonTube(gapName.c_str(), pConcentricGap->m_innerSymmetryOrder, pConcentricGap->m_outerSymmetryOrder,
                pConcentricGap->m_innerRCoordinate * m_scalingFactor, pConcentricGap->m_outerRCoordinate * m_scalingFactor,
                pConcentricGap->m_innerPhiCoordinate, pConcentricGap->m_outerPhiCoordinate, (zThick / 2.), pGapMedium);

            pGapVol->SetLineColor(1);
            pGapVol->SetFillColor(1);
            pGapVol->SetTransparency(transparency + 23);

            pMainDetectorVolume->AddNode(pGapVol, 0, new TGeoTranslation(0, 0, zMin + (zThick / 2.f)));
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::ViewEvent()
{
    InitializeEve();

    gEve->Redraw3D();
    this->Pause();

    TEveEventManager* current = gEve->GetCurrentEvent();
    if( current )
        current->SetRnrSelfChildren(kFALSE,kFALSE);

    m_openEveEvent = false;
    std::cout << "View done" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeCaloHits(const pandora::OrderedCaloHitList *const pOrderedCaloHitList, std::string name,
    TEveElement *parent, Color color)
{
    InitializeEve();

    const float markerSize(0.1);
    const std::string hitListName(name.empty() ? "Hits" : name);
    TEvePointSet *hitsMarkers = new TEvePointSet((hitListName+"_markers").c_str());

    TEveBoxSet *hits = new TEveBoxSet(name.c_str());
    hits->Reset(TEveBoxSet::kBT_FreeBox, kTRUE, 64);
    hits->AddElement(hitsMarkers);
    hits->SetOwnIds(kTRUE);

    typedef std::map<int,float> EnergyForParticleId;
    EnergyForParticleId energyDepositsForMCParticlesId;

    int numberHits = 0;
    pandora::PseudoLayer firstLayer = 0;
    pandora::PseudoLayer lastLayer = 0;

    firstLayer = pOrderedCaloHitList->begin()->first;

    float minInteractionLengthsFromIp = std::numeric_limits<float>::max();
    float maxInteractionLengthsFromIp = std::numeric_limits<float>::min();

    float energySumElectromagnetic = 0.f;
    float energySumHadronic = 0.f;

    for (pandora::OrderedCaloHitList::const_iterator iter = pOrderedCaloHitList->begin(), iterEnd = pOrderedCaloHitList->end();
         iter != iterEnd; ++iter)
    {
        lastLayer = iter->first;

        for (pandora::CaloHitList::const_iterator caloHitIter = iter->second->begin(), caloHitIterEnd = iter->second->end();
             caloHitIter != caloHitIterEnd; ++caloHitIter)
        {
            const pandora::CaloHit *pCaloHit = (*caloHitIter);
            ++numberHits;

            float interactionLengthsFromIp = pCaloHit->GetNInteractionLengthsFromIp();
            if( interactionLengthsFromIp < minInteractionLengthsFromIp )
                minInteractionLengthsFromIp = interactionLengthsFromIp;

            if( interactionLengthsFromIp > maxInteractionLengthsFromIp )
                maxInteractionLengthsFromIp = interactionLengthsFromIp;

            const pandora::MCParticle *pMCParticle = NULL;
            pCaloHit->GetMCParticle(pMCParticle);

            float hitEnergy = pCaloHit->GetElectromagneticEnergy();
            energySumElectromagnetic += hitEnergy;

            float hitEnergyHadronic = pCaloHit->GetHadronicEnergy();
            energySumHadronic += hitEnergyHadronic;

            int particleId = 0;
            if (pMCParticle)
            {
                particleId = pMCParticle->GetParticleId();
            }
            EnergyForParticleId::iterator it = energyDepositsForMCParticlesId.find( particleId );
            if( it == energyDepositsForMCParticlesId.end() )
                energyDepositsForMCParticlesId.insert( std::make_pair(particleId, hitEnergy) );
            else
            {
                float oldValue = it->second;
                it->second = oldValue + hitEnergy;
            }

            // Compute the corners of calohit calorimeter-cell, 8 corners x 3 dimensions
            Float_t corners[24];
            MakeCaloHitCell(pCaloHit, corners);

            // Supply hit marker details
            const pandora::CartesianVector position = pCaloHit->GetPositionVector() * m_scalingFactor;

            hitsMarkers->SetNextPoint(position.GetX(), position.GetY(), position.GetZ());
            hitsMarkers->SetMarkerColor(GetROOTColor(color));
            hitsMarkers->SetMarkerSize(markerSize);
            hitsMarkers->SetMarkerStyle(4);

            // Add calorimeter-cell for calo-hit
            hits->AddBox(corners);
            hits->SetPickable(kTRUE);
            hits->DigitColor(GetROOTColor(color));
            //            hits->DigitValue(90);
        }
    }


    // Build information string
    std::stringstream sstr;

    if (!name.empty())
        sstr << name << "\n";

    sstr << "--- calo-hits"
         << "\nEem=" << energySumElectromagnetic
         << "\nEhad=" << energySumHadronic
         << "\nfirst pseudo-layer=" << firstLayer
         << "\nlast  pseudo-layer=" << lastLayer
         << "\nmin intLenFromIP=" << minInteractionLengthsFromIp
         << "\nmax intLenFromIP=" << maxInteractionLengthsFromIp;
    for( EnergyForParticleId::iterator itEForId = energyDepositsForMCParticlesId.begin(), itEForIdEnd = energyDepositsForMCParticlesId.end(); 
         itEForId != itEForIdEnd; ++itEForId )
    {
        int mcPDG = itEForId->first;
        float energy = itEForId->second;
        if( mcPDG == 0 )
            sstr << "\nCaloHits w/o MC particle = " << energy << " GeV";
        else
            sstr << "\nfrom MC with PDG " << mcPDG << " = " << energy << " GeV";
    }

    hits->SetElementTitle(sstr.str().c_str());


    if (parent)
    {
        parent->AddElement(hits);
    }
    else
    {
        gEve->AddElement(hits);
        gEve->Redraw3D();
    }

    return hits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeTracks(const pandora::TrackList *const pTrackList, std::string name, TEveElement *parent, Color color)
{
    InitializeEve();

    TEveTrackList *trackList = new TEveTrackList();
    const std::string trackListName(name.empty() ? "Tracks" : name);
    trackList->SetElementNameTitle( trackListName.c_str(), trackListName.c_str() );
    trackList->SetMainColor(GetROOTColor(TEAL));

    TEveTrackPropagator *propagator = trackList->GetPropagator();
    // propagator->SetStepper(TEveTrackPropagator::kRungeKutta);

    pandora::GeometryHelper *pGeometryHelper = pandora::GeometryHelper::GetInstance();
    const float magneticField(pGeometryHelper->GetBField());

    propagator->SetMagFieldObj(new TEveMagFieldConst(0., 0., magneticField));
    propagator->SetMaxR(pGeometryHelper->GetECalBarrelParameters().GetOuterRCoordinate() * m_scalingFactor);
    propagator->SetMaxZ(pGeometryHelper->GetECalEndCapParameters().GetOuterZCoordinate() * m_scalingFactor);
    propagator->SetMaxOrbs(5);

    for (pandora::TrackList::const_iterator trackIter = pTrackList->begin(), trackIterEnd = pTrackList->end();
        trackIter != trackIterEnd; ++trackIter)
    { 
        pandora::Track *pandoraTrack = (*trackIter);

        // Extract pandora track states
        const pandora::TrackState &trackState = pandoraTrack->GetTrackStateAtStart();
        const pandora::CartesianVector &momentum = trackState.GetMomentum();
        const pandora::CartesianVector position = trackState.GetPosition() * m_scalingFactor;

        const pandora::TrackState &trackStateAtEnd = pandoraTrack->GetTrackStateAtEnd();
        const pandora::CartesianVector &momentumAtEnd = trackStateAtEnd.GetMomentum();
        const pandora::CartesianVector positionAtEnd = trackStateAtEnd.GetPosition()*m_scalingFactor;

        const pandora::TrackState &trackStateAtECal = pandoraTrack->GetTrackStateAtECal();
        const pandora::CartesianVector positionAtECal = trackStateAtECal.GetPosition()*m_scalingFactor;

        // Color assignment
        int charge = pandoraTrack->GetCharge();

        Color trackColor = color;
        if (color == AUTO)
        {
            trackColor = AZURE;

            if (charge > 0)
            {
                trackColor = RED;
            }
            else if(charge < 0)
            {
                trackColor = GREEN;
            }
        }

        // Build information string
        std::stringstream sstr;

        if (!name.empty())
            sstr << name << "\n";

        sstr << "--- track"
             << "\np=" << momentum.GetMagnitude()
             << "\nCharge=" << charge
             << "\nPDG=" << pandoraTrack->GetParticleId();

        // Create track path, note strange ALICE charge sign convention,
        // see http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=9456&p=40325&hilit=teve+histogram#p40325
        TEveRecTrack *rc = new TEveRecTrack();
        rc->fV.Set(position.GetX(),position.GetY(),position.GetZ());
        rc->fP.Set(momentum.GetX(),momentum.GetY(),momentum.GetZ());
        rc->fSign = -charge;

        TEveTrack *track = new TEveTrack(rc, propagator);
        track->SetName(sstr.str().c_str());
        track->SetTitle(sstr.str().c_str());
        track->SetLineColor(GetROOTColor(trackColor));
        track->SetLineWidth(1);
        track->SetPickable(kTRUE);

        // Create mark at track end
        TEvePathMark* pmEnd = new TEvePathMark(TEvePathMark::kReference);
        pmEnd->fV.Set(positionAtEnd.GetX(),positionAtEnd.GetY(),positionAtEnd.GetZ());
        pmEnd->fP.Set(momentumAtEnd.GetX(),momentumAtEnd.GetY(),momentumAtEnd.GetZ());
        track->AddPathMark(*pmEnd);

        // Create mark at track projection to ecal
        TEvePathMark *pmECal = new TEvePathMark(TEvePathMark::kDecay);
        pmECal->fV.Set(positionAtECal.GetX(),positionAtECal.GetY(),positionAtECal.GetZ());
        track->AddPathMark(*pmECal);

        trackList->AddElement(track);
        track->MakeTrack();
    }

    if (parent)
    {
        parent->AddElement(trackList);
    }
    else
    {
        gEve->AddElement(trackList);
        gEve->Redraw3D();
    }

    return trackList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeParticleFlowObjects(const pandora::ParticleFlowObjectList *const pPfoList, std::string name,
    TEveElement *parent, Color color, bool showAssociatedTracks, bool showFit)
{
    InitializeEve();

    TEveElement* pfoListElement = new TEveElementList();
    const std::string pfoListName(name.empty() ? "Pfos" : name);
    pfoListElement->SetElementNameTitle(pfoListName.c_str(), pfoListName.c_str());

    for (pandora::ParticleFlowObjectList::const_iterator pfoIter = pPfoList->begin(), pfoIterEnd = pPfoList->end();
        pfoIter != pfoIterEnd; ++pfoIter)
    { 
        pandora::ParticleFlowObject *pPfo = (*pfoIter);

        // Build information string
        std::stringstream sstr;
        sstr << "--- PFO"
             << "\nE=" << pPfo->GetEnergy() 
             << "\nm=" << pPfo->GetMass()
             << "\nPDG=" << pPfo->GetParticleId();

        // Default color assignment
        Color pfoColor = color;

        if (color == AUTO)
        {
            pfoColor = GetColorForPdgCode(pPfo->GetParticleId());
        }

        // show clusters and tracks
        const pandora::ClusterList clusterList(pPfo->GetClusterList());
        const pandora::TrackList trackList(pPfo->GetTrackList());

        if (clusterList.empty())
        {
            VisualizeTracks(&trackList, sstr.str().c_str(), pfoListElement, pfoColor);
        }
        else
        {
            VisualizeClusters(&clusterList, sstr.str().c_str(), pfoListElement, pfoColor, showAssociatedTracks, showFit);
        }
    }

    if (parent)
    {
        parent->AddElement(pfoListElement);
    }
    else
    {
        gEve->AddElement(pfoListElement);
        gEve->Redraw3D();
    }

    return pfoListElement;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeClusters(const pandora::ClusterList *const pClusterList, std::string name, TEveElement *parent,
    Color color, bool showAssociatedTracks, bool showFit)
{
    InitializeEve();

    TEveElement* clusterListElement = new TEveElementList();
    const std::string clusterListName(name.empty() ? "Clusters" : name);
    clusterListElement->SetElementNameTitle( name.c_str(), name.c_str());

    for (pandora::ClusterList::const_iterator clusterIter = pClusterList->begin(), clusterIterEnd = pClusterList->end();
        clusterIter != clusterIterEnd; ++clusterIter)
    {
        pandora::Cluster *pCluster = (*clusterIter);

        if (pCluster->GetNCaloHits() == 0)
            continue;

        // Color assignment
        Color clusterColor = color;

        if (color == AUTO)
        {
            const pandora::TrackList& pTrackList = pCluster->GetAssociatedTrackList();
            bool clusterHasTracks = !(pTrackList.empty());
            bool clusterIsPhoton  = pCluster->IsPhotonFast();

            if (clusterIsPhoton)
            {
                clusterColor = DARKYELLOW;
            }
            else if (clusterHasTracks)
            {
                clusterColor = MAGENTA;
            }
            else
            {
                clusterColor = LIGHTBLUE;
            }
        }

        // Build information string
        std::stringstream sstr;

        if (!name.empty())
            sstr << name << "\n";

        sstr << "--- cluster\nEem=" << pCluster->GetElectromagneticEnergy() 
             << "\nEhad=" << pCluster->GetHadronicEnergy() 
             << "\nNHits=" << pCluster->GetNCaloHits();

        // Display constituent calo hits
        const pandora::OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        TEveElement *caloHitsElement = VisualizeCaloHits(&orderedCaloHitList, sstr.str().c_str(), clusterListElement, clusterColor);

        const pandora::ClusterHelper::ClusterFitResult &fit = pCluster->GetFitToAllHitsResult();

        if (showFit && fit.IsFitSuccessful())
        {
            const pandora::CartesianVector intercept = fit.GetIntercept() * m_scalingFactor;
            pandora::CartesianVector direction = fit.GetDirection() * m_scalingFactor;

            const double length = 100;
            pandora::CartesianVector displacedStart = intercept - (direction * (length / 2));
            direction *= length;

            TEveArrow *clusterArrow = new TEveArrow(direction.GetX(), direction.GetY(), direction.GetZ(), displacedStart.GetX(),
                displacedStart.GetY(), displacedStart.GetZ());

            clusterArrow->SetConeR(0.03);
            clusterArrow->SetConeL(0.2);
            clusterArrow->SetMainColor(GetROOTColor(clusterColor));
            clusterArrow->SetPickable(kTRUE);
            clusterArrow->SetElementNameTitle(sstr.str().c_str(), sstr.str().c_str());

            caloHitsElement->AddElement(clusterArrow);
        }

        // Show tracks associated with clusters
        if (showAssociatedTracks)
        {
            const pandora::TrackList& pTrackList = pCluster->GetAssociatedTrackList();

            if (!pTrackList.empty())
            {
                VisualizeTracks(&pTrackList, sstr.str().c_str(), caloHitsElement, color);
            }
        }
    }

    if (parent)
    {
        parent->AddElement(clusterListElement);
    }
    else
    {
        gEve->AddElement(clusterListElement);
        gEve->Redraw3D();
    }

    return clusterListElement;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeCaloHitCell(const pandora::CaloHit *const pCaloHit, Float_t corners[24])
{
    pandora::CartesianVector dirU((pandora::ENDCAP == pCaloHit->GetDetectorRegion()) ? pandora::CartesianVector(0, 1, 0) : pandora::CartesianVector(0, 0, 1));
    pandora::CartesianVector normal(pCaloHit->GetNormalVector());

    pandora::CartesianVector dirV(normal.GetCrossProduct(dirU));
    const float magnitudeV(dirV.GetMagnitude());


    const pandora::CartesianVector position(pCaloHit->GetPositionVector() * m_scalingFactor);
    float u2(pCaloHit->GetCellSizeU() * m_scalingFactor / 2.0);
    float v2(pCaloHit->GetCellSizeV() * m_scalingFactor / 2.0);
    float t2(pCaloHit->GetCellThickness() * m_scalingFactor / 2.0);

    if (magnitudeV < 0.00001)
    {
        const pandora::DetectorRegion dr(pCaloHit->GetDetectorRegion() );
        std::string detectorRegion( (dr==pandora::ENDCAP ? "ENDCAP":(dr==pandora::BARREL?"BARREL":"UNKNOWN")) );
        std::cout << "PandoraMonitoring::MakeCaloHitCell, ERROR: direction vectors U and V are parallel. "
                  << "Normal("
                  << normal.GetX() << ", " << normal.GetY() << ", " << normal.GetZ() << ") "
                  << "U("
                  << dirU.GetX() << ", " << dirU.GetY() << ", " << dirU.GetZ() << ") "
                  << "V("
                  << dirV.GetX() << ", " << dirV.GetY() << ", " << dirV.GetZ() << ") "
                  << " in Detector-region: " << detectorRegion << "  caloHit-coordinates: " 
                  << position.GetX() << ", " << position.GetY() << ", " << position.GetZ()
                  << std::endl;
        normal.SetValues(1,0, 0 );
        dirU.SetValues( 0, 1, 0 );
        dirV.SetValues( 0, 0, 1 );
        
        const float size = 0.1;
        u2 = size*m_scalingFactor/2.0;
        v2 = size*m_scalingFactor/2.0;
        t2 = size*m_scalingFactor/2.0;
        //        throw std::exception();
    }
    else
    {
        dirV *= 1. / magnitudeV;
    }

    dirU *= u2;
    dirV *= v2;
    normal *= t2;

    // Compute all 8 corners for the box
    pandora::CartesianVector cornerVector[8];

    cornerVector[0] = position - dirU - dirV - normal;
    cornerVector[1] = position + dirU - dirV - normal;
    cornerVector[2] = position + dirU + dirV - normal;
    cornerVector[3] = position - dirU + dirV - normal;
    cornerVector[4] = position - dirU - dirV + normal;
    cornerVector[5] = position + dirU - dirV + normal;
    cornerVector[6] = position + dirU + dirV + normal;
    cornerVector[7] = position - dirU + dirV + normal;

    // Assign corner positions to Float_t array
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

TGeoVolume* PandoraMonitoring::MakePolygonTube(std::string name, int innerSymmetryOrder, int outerSymmetryOrder, double innerClosestDistanceToIp,
    double outerClosestDistanceToIp, double innerPhi0, double outerPhi0, double halfLength, TGeoMedium* medium)
{
    TGeoShape* inner = MakePolygonTube(innerSymmetryOrder, innerClosestDistanceToIp, innerPhi0, halfLength + 2);
    TGeoShape* outer = MakePolygonTube(outerSymmetryOrder, outerClosestDistanceToIp, outerPhi0, halfLength);

    std::string nameInner = (name + "_I");
    std::string nameOuter = (name + "_O");

    inner->SetName(nameInner.c_str());
    outer->SetName(nameOuter.c_str());

    std::string formula = nameOuter+"-"+nameInner;

    TGeoCompositeShape* cs = new TGeoCompositeShape( (name+"_shape").c_str(), formula.c_str() );
    TGeoVolume* csVol = new TGeoVolume( name.c_str(), cs, medium );

    return csVol;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TGeoShape *PandoraMonitoring::MakePolygonTube(int symmetryOrder, double closestDistanceToIp, double phi, double halfLength)
{
    if( symmetryOrder <= 2)
    {
        TGeoShape* tube = new TGeoTube(0, closestDistanceToIp, halfLength);
        return tube;
    }

    DoublePairVector vertices;
    ComputePolygonCorners(symmetryOrder, closestDistanceToIp, phi, vertices);

    const Int_t nvertices = vertices.size();
    Double_t *x = new Double_t[nvertices];
    Double_t *y = new Double_t[nvertices];

    int index = 0;
    for (DoublePairVector::iterator itCoord = vertices.begin(), itCoordEnd = vertices.end(); itCoord != itCoordEnd; ++itCoord)
    {
        x[index] = (*itCoord).first;
        y[index] = (*itCoord).second;
        ++index;
    }

    TGeoXtru *xtru = new TGeoXtru(2);

    xtru->DefinePolygon(nvertices,x,y);
    Double_t z0 = -halfLength, x0 = 0, y0 = 0;
    Double_t z1 = halfLength, x1 = 0, y1 = 0;

    Double_t scale0 = 1.0;
    xtru->DefineSection(0, z0, x0, y0, scale0); // Z position, offset and scale for first section
    xtru->DefineSection(1, z1, x1, y1, scale0); // -''- go forward

    delete[] x;
    delete[] y;

    return xtru;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::ComputePolygonCorners( int symmetryOrder, double closestDistanceToIp, double phi0, DoublePairVector &coordinates)
{
    if (symmetryOrder > 2)
    {
        static const Double_t pi(std::acos(-1.));
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

} // namespace pandora_monitoring
