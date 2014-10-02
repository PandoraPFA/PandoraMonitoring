/**
 *  @file   PandoraMonitoring/src/PandoraMonitoring.cc
 * 
 *  @brief  Implementation of the pandora monitoring class.
 * 
 *  $Log: $
 */

// Pandora include files
#include "Helpers/MCParticleHelper.h"

#include "Managers/GeometryManager.h"
#include "Managers/PluginManager.h"

#include "Objects/CaloHit.h"
#include "Objects/CartesianVector.h"
#include "Objects/Cluster.h"
#include "Objects/DetectorGap.h"
#include "Objects/Histograms.h"
#include "Objects/MCParticle.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/SubDetector.h"
#include "Objects/Track.h"
#include "Objects/Vertex.h"

#include "Pandora/Algorithm.h"

#include "Plugins/BFieldPlugin.h"

// ROOT include files
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGLViewer.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TParticlePDG.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "TEveBoxSet.h"
#include "TEveElement.h"
#include "TEveEventManager.h"
#include "TEveGeoNode.h"
#include "TEveManager.h"
#include "TEvePathMark.h"
#include "TEvePointSet.h"
#include "TEveScene.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"

#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TGeoXtru.h"

#include "PandoraMonitoring.h"

#include <algorithm>
#include <cmath>
#include <fcntl.h>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <vector>

using namespace pandora;

namespace pandora_monitoring
{

PandoraMonitoring::MonitoringInstanceMap PandoraMonitoring::m_monitoringInstanceMap;

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraMonitoring *PandoraMonitoring::GetInstance(const Pandora &pandora)
{
    MonitoringInstanceMap::const_iterator iter = m_monitoringInstanceMap.find(&pandora);

    if (m_monitoringInstanceMap.end() != iter)
        return iter->second;

    PandoraMonitoring *pPandoraMonitoring = new PandoraMonitoring(pandora);

    if (m_monitoringInstanceMap.empty())
    {
        TColor::CreateColorWheel();
        gStyle->SetPalette(1);
        gStyle->SetNumberContours(99);
    }

    if (!m_monitoringInstanceMap.insert(MonitoringInstanceMap::value_type(&pandora, pPandoraMonitoring)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return pPandoraMonitoring;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::DeleteInstance(const Pandora &pandora)
{
    MonitoringInstanceMap::iterator iter = m_monitoringInstanceMap.find(&pandora);

    if (m_monitoringInstanceMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    PandoraMonitoring *pPandoraMonitoring = iter->second;
    m_monitoringInstanceMap.erase(iter);
    delete pPandoraMonitoring;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename VariableType>
void PandoraMonitoring::SetTreeVariable(const std::string &treeName, const std::string &variableName, VariableType variable)
{
    m_pTreeWrapper->Set(treeName, variableName, variable);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::FillTree(const std::string &treeName)
{
    try
    {
        m_pTreeWrapper->Fill(treeName);
    }
    catch(TTreeWrapper::TreeNotFoundError& excpt)
    {
        std::cout << "PandoraMonitoring::FillTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::FillTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::PrintTree(const std::string &treeName)
{
    try
    {
        m_pTreeWrapper->Print(treeName);
    }
    catch(TTreeWrapper::TreeNotFoundError& excpt)
    {
        std::cout << "PandoraMonitoring::PrintTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::PrintTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::ScanTree(const std::string &treeName)
{
    try
    {
        m_pTreeWrapper->Scan(treeName);
    }
    catch(TTreeWrapper::TreeNotFoundError& excpt)
    {
        std::cout << "PandoraMonitoring::ScanTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::ScanTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::SaveTree(const std::string &treeName, const std::string &fileName, const std::string &fileOptions)
{
    try
    {
       TTree *&pTTree = m_pTreeWrapper->GetTree(treeName);

       TFile *pTFile = new TFile(fileName.c_str(), fileOptions.c_str());

       pTTree->SetDirectory(pTFile);
       pTTree->Write(TString(treeName.c_str()), TObject::kOverwrite);

       pTFile->Close();
       pTTree = NULL; // pointer not valid any more
       delete pTFile;
    }
    catch(TTreeWrapper::TreeNotFoundError& excpt)
    {
        std::cout << "PandoraMonitoring::SaveTree, error: No tree with name '" << treeName <<"' exists." << std::endl;
    }
    catch(...)
    {
        std::cout << "PandoraMonitoring::SaveTree, unknown error for tree with name '" << treeName <<"'." << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <>
void PandoraMonitoring::DrawPandoraHistogram(const pandora::Histogram &histogram, const std::string &options)
{
    const std::string histogramName("PandoraHistogram");
    TH1F *pTH1F = new TH1F(histogramName.c_str(), histogramName.c_str(), histogram.GetNBinsX(), histogram.GetXLow(), histogram.GetXHigh());

    for (int xBin = -1, xBinEnd = histogram.GetNBinsX(); xBin <= xBinEnd; ++xBin)
    {
        // Beware bin offset between ROOT and Pandora histograms
        pTH1F->SetBinContent(xBin + 1, histogram.GetBinContent(xBin));
    }

    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    pTH1F->Draw(options.c_str());
    this->Pause();

    delete pCanvas;
    delete pTH1F;
}

template <>
void PandoraMonitoring::DrawPandoraHistogram(const pandora::TwoDHistogram &twoDHistogram, const std::string &options)
{
    const std::string histogramName("PandoraHistogram");
    TH2F *pTH2F = new TH2F(histogramName.c_str(), histogramName.c_str(), twoDHistogram.GetNBinsX(), twoDHistogram.GetXLow(),
        twoDHistogram.GetXHigh(), twoDHistogram.GetNBinsY(), twoDHistogram.GetYLow(), twoDHistogram.GetYHigh());

    for (int xBin = -1, xBinEnd = twoDHistogram.GetNBinsX(); xBin <= xBinEnd; ++xBin)
    {
        for (int yBin = -1, yBinEnd = twoDHistogram.GetNBinsY(); yBin <= yBinEnd; ++yBin)
        {
            // Beware bin offset between ROOT and Pandora histograms
            pTH2F->SetBinContent(xBin + 1, yBin + 1, twoDHistogram.GetBinContent(xBin, yBin));
        }
    }

    TCanvas *pCanvas = new TCanvas("PandoraMonitoring", "PandoraMonitoring", 750, 750);
    pCanvas->SetFillColor(kWhite);
    pCanvas->SetHighLightColor(kWhite);
    pCanvas->Draw();

    pTH2F->Draw(options.c_str());
    this->Pause();

    delete pCanvas;
    delete pTH2F;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::SetEveDisplayParameters(const bool showDetectors, const DetectorView detectorView, const float transparencyThresholdE,
    const float energyScaleThresholdE)
{
    m_showDetectors = showDetectors;
    m_detectorView = detectorView;
    m_transparencyThresholdE = transparencyThresholdE;
    m_energyScaleThresholdE = energyScaleThresholdE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeMCParticles(const MCParticleList *const pMCParticleList, const std::string &name,
    TEveElement *parent, const Color color, const PandoraMonitoringApi::PdgCodeToEnergyMap *pParticleSuppressionMap)
{
    this->InitializeEve();
    MCParticleVector mcParticleVector(pMCParticleList->begin(), pMCParticleList->end());
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), PandoraMonitoring::SortMCParticlesByEnergy);

    const std::string starter("---");
    const std::string mcParticleListTitle(name.empty() ? "MCParticles" : name);
    std::string mcParticleListName(mcParticleListTitle);
    std::replace_if(mcParticleListName.begin(), mcParticleListName.end(), std::bind2nd(std::equal_to<char>(),'\n'), '/');

    TEveTrackList *pTEveTrackList = new TEveTrackList();
    pTEveTrackList->SetElementNameTitle(mcParticleListName.c_str(), mcParticleListTitle.c_str() );
    pTEveTrackList->SetMainColor(GetROOTColor(TEAL));

    float bFieldZ(0.f);
    try {bFieldZ = m_pPandora->GetPlugins()->GetBFieldPlugin()->GetBField(CartesianVector(0., 0., 0.));}
    catch (StatusCodeException &) {}

    // Initialize magnetic field for particle propagation, note strange ALICE charge sign convention,
    // see http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=9456&p=40325&hilit=teve+histogram#p40325
    TEveTrackPropagator *pTEveTrackPropagator = pTEveTrackList->GetPropagator();
    pTEveTrackPropagator->SetMagFieldObj(new TEveMagFieldConst(0., 0., -bFieldZ));
    pTEveTrackPropagator->SetMaxOrbs(5);

    try {pTEveTrackPropagator->SetMaxR(m_pPandora->GetGeometry()->GetSubDetector(HCAL_BARREL).GetOuterRCoordinate() * m_scalingFactor);}
    catch (StatusCodeException &) {}

    try {pTEveTrackPropagator->SetMaxZ(std::fabs(m_pPandora->GetGeometry()->GetSubDetector(HCAL_ENDCAP).GetOuterZCoordinate()) * m_scalingFactor);}
    catch (StatusCodeException &) {}

    for (MCParticleVector::const_iterator iter = mcParticleVector.begin(), iterEnd = mcParticleVector.end(); iter != iterEnd; ++iter)
    { 
        MCParticle *pPandoraMCParticle(*iter);

        // Does particle pass suppression conditions?
        const int particleId(pPandoraMCParticle->GetParticleId());
        const float energy(pPandoraMCParticle->GetEnergy());

        if (pParticleSuppressionMap)
        {
            PandoraMonitoringApi::PdgCodeToEnergyMap::const_iterator suppressIter = pParticleSuppressionMap->find(particleId);

            if ((suppressIter != pParticleSuppressionMap->end()) && (suppressIter->second > energy))
                continue;
        }

        TEveMCTrack *pTEveRecTrack = new TEveMCTrack();
        const CartesianVector position(pPandoraMCParticle->GetVertex() * m_scalingFactor);
        const CartesianVector &momentum(pPandoraMCParticle->GetMomentum());
        pTEveRecTrack->SetProductionVertex(position.GetX(), position.GetY(), position.GetZ(), 0.f);
        pTEveRecTrack->SetMomentum(momentum.GetX(), momentum.GetY(), momentum.GetZ(), energy);
        pTEveRecTrack->SetPdgCode(particleId);

        const int charge((pTEveRecTrack->GetPDG()) ? (static_cast<int>(pTEveRecTrack->GetPDG()->Charge()) / 3) : 0);
        const Color mcParticleColor((color < AUTO) ? color : (charge == 0) ? CYAN : (charge > 0) ? LIGHTRED : LIGHTGREEN);
        const float innerRadius(pPandoraMCParticle->GetInnerRadius());
        const float outerRadius(pPandoraMCParticle->GetOuterRadius());

        std::stringstream sstr, sstrName;
        sstr << starter << "MC particle" << "\nPDG=" << particleId << "\np=" << momentum.GetMagnitude() << "\nE=" << energy
            << "\nCharge=" << charge << "\nr_inner=" << innerRadius << "\nr_outer=" << outerRadius;
        sstrName << "MC" << "/PDG=" << particleId << "/p=" << momentum.GetMagnitude() << "/E=" << energy
            << "/Charge=" << charge << "/r_inner=" << innerRadius << "/r_outer=" << outerRadius;

        // Create particle path
        TEveTrack *pTEveTrack = new TEveTrack(pTEveRecTrack, pTEveTrackPropagator);
        pTEveTrack->SetName(sstrName.str().c_str());
        pTEveTrack->SetTitle(sstr.str().c_str());
        pTEveTrack->SetLineColor(GetROOTColor(mcParticleColor));
        pTEveTrack->SetLineWidth(1);
        pTEveTrack->SetLineStyle(2);
        pTEveTrack->SetPickable(kTRUE);

        // Create mark at end position
        TEvePathMark endPositionMark(TEvePathMark::kDecay);
        const CartesianVector positionAtEnd(pPandoraMCParticle->GetEndpoint() * m_scalingFactor);
        endPositionMark.fV.Set(positionAtEnd.GetX(), positionAtEnd.GetY(), positionAtEnd.GetZ());
        pTEveTrack->AddPathMark(endPositionMark);

        pTEveTrackList->AddElement(pTEveTrack);
        pTEveTrack->MakeTrack();
    }

    if (parent)
    {
        parent->AddElement(pTEveTrackList);
    }
    else
    {
        m_pEveManager->GetCurrentEvent()->AddElement(pTEveTrackList);
        m_pEveManager->Redraw3D();
    }

    return pTEveTrackList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeTracks(const TrackList *const pTrackList, const std::string &name, TEveElement *parent, const Color color)
{
    this->InitializeEve();
    TrackVector trackVector(pTrackList->begin(), pTrackList->end());
    std::sort(trackVector.begin(), trackVector.end(), PandoraMonitoring::SortTracksByMomentum);

    const std::string starter("---");
    const std::string trackListTitle(name.empty() ? "Tracks" : name);
    std::string trackListName(trackListTitle);
    if (trackListName.find(starter) != std::string::npos)
        trackListName.replace(trackListName.find(starter), starter.length(), "Tracks/");
    std::replace_if(trackListName.begin(), trackListName.end(), std::bind2nd(std::equal_to<char>(),'\n'), '/');

    TEveTrackList *pTEveTrackList = new TEveTrackList();
    pTEveTrackList->SetElementNameTitle(trackListName.c_str(), trackListTitle.c_str() );
    pTEveTrackList->SetMainColor(GetROOTColor(TEAL));

    float bFieldZ(0.f);
    try {bFieldZ = m_pPandora->GetPlugins()->GetBFieldPlugin()->GetBField(CartesianVector(0., 0., 0.));}
    catch (StatusCodeException &) {}

    // Initialize magnetic field for particle propagation, note strange ALICE charge sign convention,
    // see http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=9456&p=40325&hilit=teve+histogram#p40325
    TEveTrackPropagator *pTEveTrackPropagator = pTEveTrackList->GetPropagator();
    pTEveTrackPropagator->SetMagFieldObj(new TEveMagFieldConst(0., 0., -bFieldZ));
    pTEveTrackPropagator->SetMaxOrbs(5);

    try {pTEveTrackPropagator->SetMaxR(m_pPandora->GetGeometry()->GetSubDetector(ECAL_BARREL).GetOuterRCoordinate() * m_scalingFactor);}
    catch (StatusCodeException &) {}

    try {pTEveTrackPropagator->SetMaxZ(std::fabs(m_pPandora->GetGeometry()->GetSubDetector(ECAL_ENDCAP).GetOuterZCoordinate()) * m_scalingFactor);}
    catch (StatusCodeException &) {}

    for (TrackVector::const_iterator iter = trackVector.begin(), iterEnd = trackVector.end(); iter != iterEnd; ++iter)
    { 
        Track *pPandoraTrack(*iter);

        const CartesianVector &momentum(pPandoraTrack->GetTrackStateAtStart().GetMomentum());
        const int charge(pPandoraTrack->GetCharge());
        const Color trackColor((color < AUTO) ? color : (charge == 0) ? AZURE : (charge > 0) ? RED : GREEN);

        std::stringstream sstr, sstrName;
        sstr << starter << "Track" << "\np=" << momentum.GetMagnitude() << "\nCharge=" << charge << "\nPDG=" << pPandoraTrack->GetParticleId();
        sstrName << "Track" << "/p=" << momentum.GetMagnitude() << "/Charge=" << charge << "/PDG=" << pPandoraTrack->GetParticleId();

        try
        {
            const MCParticle *pMCParticle(pPandoraTrack->GetMainMCParticle());
            const int mcPdg(pMCParticle->GetParticleId());
            sstr << "\nPDG_MC=" << mcPdg;
            sstrName << "/PDG_MC=" << mcPdg;
        }
        catch (StatusCodeException &)
        {
        }

        // Create track path
        TEveRecTrack *pTEveRecTrack = new TEveRecTrack();
        const CartesianVector position(pPandoraTrack->GetTrackStateAtStart().GetPosition() * m_scalingFactor);
        pTEveRecTrack->fV.Set(position.GetX(), position.GetY(), position.GetZ());
        pTEveRecTrack->fP.Set(momentum.GetX(), momentum.GetY(), momentum.GetZ());
        pTEveRecTrack->fSign = charge;

        TEveTrack *pTEveTrack = new TEveTrack(pTEveRecTrack, pTEveTrackPropagator);
        pTEveTrack->SetName(sstrName.str().c_str());
        pTEveTrack->SetTitle(sstr.str().c_str());
        pTEveTrack->SetLineColor(GetROOTColor(trackColor));
        pTEveTrack->SetLineWidth(1);
        pTEveTrack->SetPickable(kTRUE);

        // Create mark at track end
        TEvePathMark endPositionMark(TEvePathMark::kReference);
        const CartesianVector &momentumAtEnd(pPandoraTrack->GetTrackStateAtEnd().GetMomentum());
        const CartesianVector positionAtEnd(pPandoraTrack->GetTrackStateAtEnd().GetPosition() * m_scalingFactor);
        endPositionMark.fV.Set(positionAtEnd.GetX(), positionAtEnd.GetY(), positionAtEnd.GetZ());
        endPositionMark.fP.Set(momentumAtEnd.GetX(), momentumAtEnd.GetY(), momentumAtEnd.GetZ());
        pTEveTrack->AddPathMark(endPositionMark);

        // Create mark at track projection to calorimeter
        TEvePathMark calorimeterPositionMark(TEvePathMark::kDecay);
        const CartesianVector positionAtCalorimeter(pPandoraTrack->GetTrackStateAtCalorimeter().GetPosition() * m_scalingFactor);
        calorimeterPositionMark.fV.Set(positionAtCalorimeter.GetX(), positionAtCalorimeter.GetY(), positionAtCalorimeter.GetZ());
        pTEveTrack->AddPathMark(calorimeterPositionMark);

        pTEveTrackList->AddElement(pTEveTrack);
        pTEveTrack->MakeTrack();
    }

    if (parent)
    {
        parent->AddElement(pTEveTrackList);
    }
    else
    {
        m_pEveManager->GetCurrentEvent()->AddElement(pTEveTrackList);
        m_pEveManager->Redraw3D();
    }

    return pTEveTrackList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeCaloHits(const CaloHitList *const pCaloHitList, const std::string &name, TEveElement *parent, const Color color)
{
    this->InitializeEve();

    const std::string starter("---");
    const std::string caloHitListTitle(name.empty() ? "CaloHits" : name);
    std::string caloHitListName(caloHitListTitle);
    if (caloHitListName.find(starter) != std::string::npos)
        caloHitListName.replace(caloHitListName.find(starter), starter.length(), "CaloHits/");
    std::replace_if(caloHitListName.begin(), caloHitListName.end(), std::bind2nd(std::equal_to<char>(),'\n'), '/');

    TEveElement *pCaloHitVectorElement = new TEveElementList();
    pCaloHitVectorElement->SetElementNameTitle(caloHitListName.c_str(), caloHitListTitle.c_str());

    TEveBoxSet *hits = new TEveBoxSet(name.c_str());
    hits->Reset(TEveBoxSet::kBT_FreeBox, kTRUE, 64);
    hits->SetOwnIds(kTRUE);
    hits->SetPickable(kTRUE);
    hits->SetMainTransparency(5);
    hits->SetAntiFlick(kTRUE);
    pCaloHitVectorElement->AddElement(hits);

    static Double_t r[] = {1., 1.}, g[] = {1., 0.}, b[] = {0., 0.}, stop[] = {0., 1.};
    static const Int_t firstIndex(TColor::CreateGradientColorTable(2, stop, r, g, b, 256));

    Int_t customPalette[256];
    for (int index = 0; index < 256; ++index)
        customPalette[index] = firstIndex + index;

    PandoraMonitoringApi::PdgCodeToEnergyMap pdgCodeToEnergyMap;
    unsigned int firstLayer = std::numeric_limits<unsigned int>::max();
    unsigned int lastLayer = 0;

    float energySumElectromagnetic = 0.f;
    float energySumHadronic = 0.f;

    static int colorIter = RED;
    if (++colorIter >= AUTO)
        colorIter = RED;

    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        const CaloHit *pCaloHit(*iter);

        // Determing extremal pseudolayers
        const unsigned int pseudoLayer(pCaloHit->GetPseudoLayer());

        if (pseudoLayer > lastLayer)
            lastLayer = pseudoLayer;

        if (pseudoLayer < firstLayer)
            firstLayer = pseudoLayer;

        // Energy properties
        const float hitEnergy(pCaloHit->GetElectromagneticEnergy());
        energySumElectromagnetic += hitEnergy;

        const float hitEnergyHadronic(pCaloHit->GetHadronicEnergy());
        energySumHadronic += hitEnergyHadronic;

        // MC particle id
        int particleId = 0;

        try
        {
           const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
           particleId = pMCParticle->GetParticleId();
        }
        catch (StatusCodeException &)
        {
        }

        PandoraMonitoringApi::PdgCodeToEnergyMap::iterator it = pdgCodeToEnergyMap.find(particleId);

        if (pdgCodeToEnergyMap.end() == it)
        {
            pdgCodeToEnergyMap.insert(std::make_pair(particleId, hitEnergy));
        }
        else
        {
            float oldValue = it->second;
            it->second = oldValue + hitEnergy;
        }

        // Supply hit marker details
        EColor hitColor = GetROOTColor(color);

        if (color == AUTOID)
        {
            hitColor = GetROOTColor(GetColorForPdgCode(particleId));
        }
        else if (color == AUTOITER)
        {
            hitColor = GetROOTColor(Color(colorIter));
        }
        else if (color == AUTOENERGY)
        {
            unsigned int customColorIndex = 0;

            if (m_energyScaleThresholdE > 0.f)
                customColorIndex = std::min(255, static_cast<int>(255.f * (hitEnergy / m_energyScaleThresholdE)));

            hitColor = EColor(customPalette[customColorIndex]);
        }

        // Compute the corners of calohit calorimeter-cell, 8 corners x 3 dimensions
        float corners[24];
        this->MakeCaloHitCell(pCaloHit, corners);
        hits->AddBox(corners);

        char transparency = 0;
        if (m_transparencyThresholdE > 0.f)
            transparency = static_cast<char>(std::max(0, 255 - static_cast<int>(255.f * (hitEnergy / m_transparencyThresholdE))));

        hits->DigitColor(hitColor, transparency);
    }

    std::stringstream sstr, sstrName;
    sstr << starter << "CaloHits" << "\nEem=" << energySumElectromagnetic << "\nEhad=" << energySumHadronic << "\nInnerLayer=" << firstLayer
        << "\nOuterLayer=" << lastLayer;
    sstrName << "CaloHits" << "/Eem=" << energySumElectromagnetic << "/Ehad=" << energySumHadronic << "/InnerLayer=" << firstLayer
        << "/OuterLayer=" << lastLayer;

    for (PandoraMonitoringApi::PdgCodeToEnergyMap::const_iterator iter = pdgCodeToEnergyMap.begin(), iterEnd = pdgCodeToEnergyMap.end(); iter != iterEnd; ++iter)
    {
        const int mcPDG(iter->first);
        const float energy(iter->second);

        if (0 == mcPDG)
        {
            sstr << "\nNoMC=" << energy << "GeV";
            sstrName << "/NoMC=" << energy << "GeV";
        }
        else
        {
            sstr << "\nFromMC(PDG " << mcPDG << ")=" << energy << "GeV";
            sstrName << "/FromMC(PDG " << mcPDG << ")=" << energy << "GeV";
        }
    }

    hits->SetElementName(sstrName.str().c_str());
    hits->SetElementTitle(sstr.str().c_str());

    if (parent)
    {
        parent->AddElement(pCaloHitVectorElement);
    }
    else
    {
        m_pEveManager->GetCurrentEvent()->AddElement(pCaloHitVectorElement);
        m_pEveManager->Redraw3D();
    }

    return pCaloHitVectorElement;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeClusters(const ClusterList *const pClusterList, const std::string &name, TEveElement *parent,
    const Color color, bool showAssociatedTracks)
{
    this->InitializeEve();
    ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), PandoraMonitoring::SortClustersByHadronicEnergy);

    const std::string starter("---");
    const std::string clusterListTitle(name.empty() ? "Clusters" : name);
    std::string clusterListName(clusterListTitle);
    if (clusterListName.find(starter) != std::string::npos)
        clusterListName.replace(clusterListName.find(starter), starter.length(), "Clusters/");
    std::replace_if(clusterListName.begin(), clusterListName.end(), std::bind2nd(std::equal_to<char>(),'\n'), '/');

    TEveElement *pClusterVectorElement = new TEveElementList();
    pClusterVectorElement->SetElementNameTitle(clusterListName.c_str(), clusterListTitle.c_str());

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);

        if (pCluster->GetNCaloHits() == 0)
            continue;

        std::stringstream sstr;
        sstr << starter << "Cluster\nEem(corr)=" << pCluster->GetElectromagneticEnergy() << "\nEhad(corr)=" << pCluster->GetHadronicEnergy()
            << "\nNHits=" << pCluster->GetNCaloHits() << "\nInnerHitType=" << this->GetHitTypeString(pCluster->GetInnerLayerHitType())
            << "\nOuterHitType=" << this->GetHitTypeString(pCluster->GetOuterLayerHitType());

        const Color clusterColor((color != AUTO) ? color : (pCluster->GetAssociatedTrackList().empty()) ? LIGHTBLUE : MAGENTA);

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        TEveElement *pCaloHitsElement = VisualizeCaloHits(&caloHitList, sstr.str().c_str(), pClusterVectorElement, clusterColor);

        if (showAssociatedTracks && !pCluster->GetAssociatedTrackList().empty())
        {
            TEveElement *pTrackParentElement(parent ? pClusterVectorElement : pCaloHitsElement);
            (void) VisualizeTracks(&(pCluster->GetAssociatedTrackList()), sstr.str().c_str(), pTrackParentElement, clusterColor);
        }
    }

    if (parent)
    {
        parent->AddElement(pClusterVectorElement);
    }
    else
    {
        m_pEveManager->GetCurrentEvent()->AddElement(pClusterVectorElement);
        m_pEveManager->Redraw3D();
    }

    return pClusterVectorElement;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeParticleFlowObjects(const PfoList *const pPfoList, const std::string &name, TEveElement *parent,
    const Color color, bool showVertices, bool displayPfoHierarchy)
{
    this->InitializeEve();
    PfoVector pfoVector(pPfoList->begin(), pPfoList->end());
    std::sort(pfoVector.begin(), pfoVector.end(), PandoraMonitoring::SortPfosByEnergy);

    const std::string starter("---");
    const std::string pfoListTitle(name.empty() ? "Pfos" : name);
    std::string pfoListName(pfoListTitle);
    std::replace_if(pfoListName.begin(), pfoListName.end(), std::bind2nd(std::equal_to<char>(),'\n'), '/');

    TEveElement *pPfoVectorElement = new TEveElementList();
    pPfoVectorElement->SetElementNameTitle(pfoListName.c_str(), pfoListTitle.c_str());

    for (PfoVector::const_iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo(*iter);

        if (displayPfoHierarchy && !parent && (pPfo->GetNParentPfos() != 0))
            continue;

        std::stringstream sstr, sstrName;
        sstr << starter << "PFO" << "\nE=" << pPfo->GetEnergy() << "\nm=" << pPfo->GetMass() << "\nPDG=" << pPfo->GetParticleId();
        sstrName << "PFO" << "/E=" << pPfo->GetEnergy() << "/m=" << pPfo->GetMass() << "/PDG=" << pPfo->GetParticleId();

        TEveElement *pPfoElement = new TEveElementList();
        pPfoElement->SetElementNameTitle(sstrName.str().c_str(), sstrName.str().c_str());
        pPfoVectorElement->AddElement(pPfoElement);

        const Color pfoColor((color != AUTO) ? color : this->GetColorForPdgCode(pPfo->GetParticleId()));

        const ClusterList &clusterList(pPfo->GetClusterList());
        const TrackList &trackList(pPfo->GetTrackList());

        if (!clusterList.empty())
            (void) VisualizeClusters(&clusterList, sstr.str().c_str(), pPfoElement, pfoColor, false);

        if (!trackList.empty())
            (void) VisualizeTracks(&trackList, sstr.str().c_str(), pPfoElement, pfoColor);

        if (showVertices && !pPfo->GetVertexList().empty())
            (void) VisualizeVertices(&(pPfo->GetVertexList()), "VertexList", pPfoElement, pfoColor);

        if (displayPfoHierarchy && !pPfo->GetDaughterPfoList().empty())
            (void) VisualizeParticleFlowObjects(&(pPfo->GetDaughterPfoList()), "DaughterPfos", pPfoElement, pfoColor, showVertices, displayPfoHierarchy);
    }

    if (parent)
    {
        parent->AddElement(pPfoVectorElement);
    }
    else
    {
        m_pEveManager->GetCurrentEvent()->AddElement(pPfoVectorElement);
        m_pEveManager->Redraw3D();
    }

    return pPfoVectorElement;
}


//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::VisualizeVertices(const pandora::VertexList *const pVertexList, const std::string &name, TEveElement *parent,
    const Color color)
{
    this->InitializeEve();

    const std::string starter("---");
    TEvePointSet *pTEvePointSet = new TEvePointSet(name.c_str(), 1);
    pTEvePointSet->SetOwnIds(kTRUE);

    for (VertexList::const_iterator iter = pVertexList->begin(), iterEnd = pVertexList->end(); iter != iterEnd; ++iter)
    {
        Vertex *pVertex = *iter;
        const CartesianVector &vertexPosition(pVertex->GetPosition());

        std::stringstream sstr;
        sstr << starter << "Vertex\n" << vertexPosition;

        (void) AddMarkerToVisualization(&vertexPosition, sstr.str(), pTEvePointSet, color, 1);
    }

    if (parent)
    {
        parent->AddElement(pTEvePointSet);
    }
    else
    {
        m_pEveManager->GetCurrentEvent()->AddElement(pTEvePointSet);
        m_pEveManager->Redraw3D();
    }

    return pTEvePointSet;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TEveElement *PandoraMonitoring::AddMarkerToVisualization(const CartesianVector *const pMarkerPoint, const std::string &name, TEveElement *parent,
    const Color color, const unsigned int markerSize)
{
    this->InitializeEve();

    const std::string markerTitle(name.empty() ? "Marker" : name);

    TEvePointSet *pTEvePointSet = new TEvePointSet(markerTitle.c_str(), 1);
    pTEvePointSet->SetOwnIds(kTRUE);
    pTEvePointSet->SetPoint(0, pMarkerPoint->GetX() * m_scalingFactor, pMarkerPoint->GetY() * m_scalingFactor, pMarkerPoint->GetZ() * m_scalingFactor);

    const Color chosenColor((color < AUTO) ? color : ORANGE);
    pTEvePointSet->SetMarkerColor(GetROOTColor(chosenColor));
    pTEvePointSet->SetMarkerSize(markerSize);
    pTEvePointSet->SetMarkerStyle(20);

    if (parent)
    {
        parent->AddElement(pTEvePointSet);
    }
    else
    {
        m_pEveManager->GetCurrentEvent()->AddElement(pTEvePointSet);
        m_pEveManager->Redraw3D();
    }

    return pTEvePointSet;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::ViewEvent()
{
    this->InitializeEve();

    m_pEveManager->Redraw3D();
    this->Pause();

    m_pEveManager->GetCurrentEvent()->SetRnrSelfChildren(kFALSE,kFALSE);

    m_openEveEvent = false;
    std::cout << "View done" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::Pause() const
{
    std::cout << "Press return to continue ..." << std::endl;
    int flag = fcntl(1, F_GETFL, 0);

    int key = 0;
    while(true)
    {
        gSystem->ProcessEvents();
        (void) fcntl(1, F_SETFL, flag | O_NONBLOCK);
        key = getchar();

        if((key == '\n') || (key == '\r'))
            break;

        usleep(1000);
    }

    (void) fcntl(1, F_SETFL, flag);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PandoraMonitoring::SortClustersByHadronicEnergy(const Cluster *const pLhs, const Cluster *const pRhs)
{
    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PandoraMonitoring::SortMCParticlesByEnergy(const MCParticle *const pLhs, const MCParticle *const pRhs)
{
    return (pLhs->GetEnergy() > pRhs->GetEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PandoraMonitoring::SortTracksByMomentum(const Track *const pLhs, const Track *const pRhs)
{
    return (pLhs->GetMomentumAtDca().GetMagnitudeSquared() > pRhs->GetMomentumAtDca().GetMagnitudeSquared());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PandoraMonitoring::SortPfosByEnergy(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    return (pLhs->GetEnergy() > pRhs->GetEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraMonitoring::PandoraMonitoring(const Pandora &pandora) :
    m_pPandora(&pandora),
    m_pApplication(NULL),
    m_pEveManager(NULL),
    m_pTreeWrapper(NULL),
    m_scalingFactor(0.1f),
    m_openEveEvent(false),
    m_eventDisplayCounter(0.f),
    m_transparencyThresholdE(-1.f),
    m_energyScaleThresholdE(-1.f),
    m_showDetectors(false),
    m_detectorView(DETECTOR_VIEW_DEFAULT)
{
    if (gApplication)
    {
        m_pApplication = gApplication;

        if (m_pApplication->TestBit(TApplication::kDefaultApplication))
            std::cout << "PandoraMonitoring, only able to use default TApplication (limited functionality)." << std::endl;
    }
    else
    {
        int argc = 0;
        char *argv = (char *)"";
        m_pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
        m_pApplication->SetReturnFromRun(kTRUE);
    }
   
    m_pTreeWrapper = new TTreeWrapper;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PandoraMonitoring::~PandoraMonitoring()
{
    delete m_pTreeWrapper;

    if (NULL != m_pEveManager)
    {
        delete m_pEveManager;
        gSystem->ProcessEvents();
    }

    if (m_monitoringInstanceMap.empty())
        m_pApplication->Terminate(0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::MakeCaloHitCell(const CaloHit *const pCaloHit, float corners[24])
{
    CartesianPointList cartesianPointList;
    pCaloHit->GetCellCorners(cartesianPointList);

    if (8 != cartesianPointList.size())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    unsigned int counter(0);

    for (CartesianPointList::iterator iter = cartesianPointList.begin(), iterEnd = cartesianPointList.end(); iter != iterEnd; ++iter)
    {
        CartesianVector &corner = *iter;
        corner *= m_scalingFactor;
        corners[counter++] = corner.GetX();
        corners[counter++] = corner.GetY();
        corners[counter++] = corner.GetZ();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::ComputePolygonCorners(int symmetryOrder, double closestDistanceToIp, double phi0, DoublePairVector &coordinates)
{
    if (symmetryOrder > 2)
    {
        static const double pi(std::acos(-1.));
        const double x0(-1. * closestDistanceToIp * tan(pi / double(symmetryOrder)));
        const double y0(closestDistanceToIp);

        for (int i = 0; i < symmetryOrder; ++i)
        {
            double theta = 0.f; 
            theta = phi0 + (2 * pi * double(i) / double(symmetryOrder));

            double x = x0 * cos(theta) + y0 * sin(theta);
            double y = y0 * cos(theta) - x0 * sin(theta);
            coordinates.push_back(std::pair<double, double>(x, y));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

TGeoVolume *PandoraMonitoring::MakePolygonTube(std::string name, int innerSymmetryOrder, int outerSymmetryOrder, double innerClosestDistanceToIp,
    double outerClosestDistanceToIp, double innerPhi0, double outerPhi0, double halfLength, TGeoMedium* medium)
{
    TGeoShape *pInnerTGeoShape = MakePolygonTube(innerSymmetryOrder, innerClosestDistanceToIp, innerPhi0, halfLength + 2);
    TGeoShape *pOuterTGeoShape = MakePolygonTube(outerSymmetryOrder, outerClosestDistanceToIp, outerPhi0, halfLength);

    std::string nameInner = (name + "_I");
    std::string nameOuter = (name + "_O");

    pInnerTGeoShape->SetName(nameInner.c_str());
    pOuterTGeoShape->SetName(nameOuter.c_str());

    std::string formula = nameOuter + "-" + nameInner;

    TGeoCompositeShape *pTGeoCompositeShape = new TGeoCompositeShape((name + "_shape").c_str(), formula.c_str());
    TGeoVolume *pTGeoVolume = new TGeoVolume(name.c_str(), pTGeoCompositeShape, medium);

    return pTGeoVolume;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TGeoShape *PandoraMonitoring::MakePolygonTube(int symmetryOrder, double closestDistanceToIp, double phi, double halfLength)
{
    if (symmetryOrder <= 2)
    {
        TGeoShape *pTGeoShape = new TGeoTube(0, closestDistanceToIp, halfLength);
        return pTGeoShape;
    }

    DoublePairVector vertices;
    ComputePolygonCorners(symmetryOrder, closestDistanceToIp, phi, vertices);

    const int nvertices(vertices.size());
    double *x = new double[nvertices];
    double *y = new double[nvertices];

    int index = 0;
    for (DoublePairVector::iterator itCoord = vertices.begin(), itCoordEnd = vertices.end(); itCoord != itCoordEnd; ++itCoord)
    {
        x[index] = (*itCoord).first;
        y[index] = (*itCoord).second;
        ++index;
    }

    TGeoXtru *pTGeoXtru = new TGeoXtru(2);

    pTGeoXtru->DefinePolygon(nvertices,x,y);
    double z0 = -halfLength, x0 = 0, y0 = 0;
    double z1 = halfLength, x1 = 0, y1 = 0;

    double scale0 = 1.0;
    pTGeoXtru->DefineSection(0, z0, x0, y0, scale0); // Z position, offset and scale for first section
    pTGeoXtru->DefineSection(1, z1, x1, y1, scale0); // -''- go forward

    delete[] x;
    delete[] y;

    return pTGeoXtru;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string PandoraMonitoring::GetHitTypeString(const pandora::HitType hitType) const
{
    switch (hitType)
    {
        case TRACKER : return "TRACKER";
        case ECAL : return "ECAL";
        case HCAL : return "HCAL";
        case MUON : return "MUON";
        case TPC_VIEW_U : return "TPC_VIEW_U";
        case TPC_VIEW_V : return "TPC_VIEW_V";
        case TPC_VIEW_W : return "TPC_VIEW_W";
        case TPC_3D : return "TPC_3D";
        case HIT_CUSTOM : return "HIT_CUSTOM";
        default : return "Unknown";
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::InitializeEve(Char_t transparency)
{
    if (NULL == m_pEveManager)
    {
        try
        {
            std::cout << "PandoraMonitoring::InitializeEve(): ";
            const char *pDisplay(::getenv("DISPLAY"));

            if (NULL == pDisplay)
            {
                std::cout << "DISPLAY environment not set" << std::endl;
            }
            else
            {
                std::cout << "DISPLAY environment set to " << pDisplay << std::endl;
            }

            m_pEveManager = TEveManager::Create();
        }
        catch (TEveException &tEveException1)
        {
            std::cout << "PandoraMonitoring::InitializeEve(): Caught TEveException: " << tEveException1.what() << std::endl;

            try
            {
                std::cout << "PandoraMonitoring::InitializeEve(): Attempt to release ROOT from batch mode." << std::endl;
                gROOT->SetBatch(kFALSE);
                m_pEveManager = TEveManager::Create();
            }
            catch (TEveException &tEveException2)
            {
                std::cout << "PandoraMonitoring::InitializeEve(): Caught TEveException: " << tEveException2.what() << std::endl;
                throw tEveException2;
            }
        }

        if (NULL == m_pEveManager)
        {
            std::cout << "PandoraMonitoring, unable to create TEveManager. Bailing." << std::endl;
            throw std::exception();
        }

        if (m_showDetectors && !m_pPandora->GetGeometry()->GetSubDetectorMap().empty())
        {
            TGeoManager *pGeoManager = (NULL != gGeoManager) ? gGeoManager : new TGeoManager("DetectorGeometry", "detector geometry");

            TGeoMaterial *pVacuumMaterial = new TGeoMaterial("Vacuum", 0, 0, 0); // dummy material
            TGeoMaterial *pAluminiumMaterial = new TGeoMaterial("Aluminium", 26.98, 13, 2.7); // dummy material
            TGeoMedium *pVacuum = new TGeoMedium("Vacuum", 1, pVacuumMaterial);
            TGeoMedium *pAluminium = new TGeoMedium("Aluminium", 2, pAluminiumMaterial);

            const bool topVolumeExists(NULL != pGeoManager->GetTopVolume());
            TGeoVolume *pMainDetectorVolume = topVolumeExists ? pGeoManager->GetTopVolume() : pGeoManager->MakeBox("Detector", pVacuum, 1000., 1000., 100.);

            this->InitializeSubDetectors(pMainDetectorVolume, pAluminium, transparency);
            this->InitializeGaps(pMainDetectorVolume, pVacuum, transparency);

            if (!topVolumeExists)
            {
                pGeoManager->SetTopVolume(pMainDetectorVolume);
                TGeoNode *pTGeoNode = pGeoManager->GetTopNode();
                TEveGeoTopNode *pTEveGeoTopNode = new TEveGeoTopNode(pGeoManager, pTGeoNode);
                pTEveGeoTopNode->SetVisLevel(1);
                pTEveGeoTopNode->GetNode()->GetVolume()->SetVisibility(kFALSE);
                m_pEveManager->AddGlobalElement(pTEveGeoTopNode);
            }

            m_pEveManager->FullRedraw3D(kTRUE);
        }

        TGLViewer *pTGLViewer = m_pEveManager->GetDefaultGLViewer();
        pTGLViewer->ColorSet().Background().SetColor(kWhite);
        
        if (DETECTOR_VIEW_XZ == m_detectorView)
        {
            pTGLViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
        }
        else if (DETECTOR_VIEW_XY == m_detectorView)
        {
            pTGLViewer->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
        }
    }
    
    if (!m_openEveEvent)
    {
        std::stringstream sstr;
        sstr << "Event Display " << m_eventDisplayCounter++;
        m_pEveManager->AddEvent(new TEveEventManager(sstr.str().c_str(),sstr.str().c_str()));
        m_openEveEvent = true;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::InitializeSubDetectors(TGeoVolume *pMainDetectorVolume, TGeoMedium *pSubDetectorMedium, Char_t transparency)
{
    const SubDetectorMap &subDetectorMap(m_pPandora->GetGeometry()->GetSubDetectorMap());

    typedef std::set<std::string> StringSet;
    StringSet setInvisible;

    setInvisible.insert("Coil");
    setInvisible.insert("Tracker");
    setInvisible.insert("MuonBarrel");

    // Set muon endcap invisible if barrel is initialized. In case of a test beam set up without barrel the muon endcap (=tail catcher) is drawn
    if (subDetectorMap.count("MuonBarrel"))
        setInvisible.insert("MuonEndCap");

    int color(2);
    for (SubDetectorMap::const_iterator iter = subDetectorMap.begin(); iter != subDetectorMap.end(); ++iter)
    {
        bool isLeft = true;
        for (int lr = 0; lr <= 1; ++lr)
        {
            const SubDetector *const pSubDetector(iter->second);

            if (isLeft && !pSubDetector->IsMirroredInZ())
            {
                isLeft = false;
                continue;
            }

            const std::string name(iter->first);
            StringSet::iterator itSetInvisible = setInvisible.find(name);
            const bool drawInvisible = (itSetInvisible != setInvisible.end() ? true : false);

            std::stringstream sstr;
            sstr << name;
            sstr << (isLeft? "_left" : "_right");

            TGeoVolume *pSubDetVol = NULL;

            const int sign(isLeft ? -1 : 1);
            double zMin = pSubDetector->GetInnerZCoordinate() * m_scalingFactor;
            double zMax = pSubDetector->GetOuterZCoordinate() * m_scalingFactor;
            double zThick = zMax - zMin;
            zMin *= sign;
            zMax *= sign;

            pSubDetVol = MakePolygonTube(sstr.str().c_str(), pSubDetector->GetInnerSymmetryOrder(), pSubDetector->GetOuterSymmetryOrder(),
                pSubDetector->GetInnerRCoordinate() * m_scalingFactor, pSubDetector->GetOuterRCoordinate() * m_scalingFactor,
                pSubDetector->GetInnerPhiCoordinate(), pSubDetector->GetOuterPhiCoordinate(), (zThick / 2.f), pSubDetectorMedium);

            pSubDetVol->SetLineColor(GetROOTColor(Color(color)));
            pSubDetVol->SetFillColor(GetROOTColor(Color(color)));
            pSubDetVol->SetTransparency(transparency);

            if (drawInvisible)
                pSubDetVol->SetVisibility(kFALSE);

            const double zPosition(zMin + sign * (zThick / 2.f));
            pMainDetectorVolume->AddNode(pSubDetVol, 0, new TGeoTranslation(0, 0, zPosition));
            isLeft = false;
        }
        ++color;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PandoraMonitoring::InitializeGaps(TGeoVolume *pMainDetectorVolume, TGeoMedium *pGapMedium, Char_t transparency)
{
    const DetectorGapList &detectorGapList(m_pPandora->GetGeometry()->GetDetectorGapList());
    unsigned int gapCounter(0);

    for (DetectorGapList::const_iterator iter = detectorGapList.begin(), iterEnd = detectorGapList.end(); iter != iterEnd; ++iter)
    {
        std::string gapName("gap" + TypeToString(gapCounter++));

        BoxGap *pBoxGap = NULL;
        pBoxGap = dynamic_cast<BoxGap *>(*iter);

        if (NULL != pBoxGap)
        {
            TGeoShape *pGapShape = new TGeoBBox(gapName.c_str(), 0.5f * pBoxGap->GetSide1().GetMagnitude() * m_scalingFactor,
                0.5f * pBoxGap->GetSide2().GetMagnitude() * m_scalingFactor, 0.5f * pBoxGap->GetSide3().GetMagnitude() * m_scalingFactor);

            TGeoVolume *pGapVol = new TGeoVolume(gapName.c_str(), pGapShape, pGapMedium);

            static const float pi(std::acos(-1.));
            float correction(0.f);

            try
            {
                // TODO Remove ILD-specific correction, required for endcap box gaps that do not point back to origin in xy plane.
                //      Pandora gaps are self-describing (four vectors), but this does not map cleanly to TGeoBBox class.
                //      Best solution may be to move to different root TGeoShape.
                const float vertexZ(pBoxGap->GetVertex().GetZ());
                const float hcalEndCapInnerZ(std::fabs(m_pPandora->GetGeometry()->GetSubDetector(HCAL_ENDCAP).GetInnerZCoordinate()));
                correction = ((std::fabs(vertexZ) < hcalEndCapInnerZ) ? 0 : ((vertexZ > 0) ? pi / 4.f : -pi / 4.f));
            }
            catch (StatusCodeException &)
            {
            }

            const float phi(correction + std::atan2(pBoxGap->GetVertex().GetX(), pBoxGap->GetVertex().GetY()));

            const TGeoTranslation trans("trans",
                ( 0.5f * pBoxGap->GetSide1().GetMagnitude() * std::cos(phi) + 0.5f * pBoxGap->GetSide2().GetMagnitude() * std::sin(phi) + pBoxGap->GetVertex().GetX()) * m_scalingFactor,
                (-0.5f * pBoxGap->GetSide1().GetMagnitude() * std::sin(phi) + 0.5f * pBoxGap->GetSide2().GetMagnitude() * std::cos(phi) + pBoxGap->GetVertex().GetY()) * m_scalingFactor,
                ( 0.5f * (2.f * pBoxGap->GetVertex().GetZ() + pBoxGap->GetSide3().GetZ()) * m_scalingFactor));

            const TGeoRotation rot("rot", -180.f * phi / pi, 0, 0);

            pGapVol->SetLineColor(1);
            pGapVol->SetFillColor(1);
            pGapVol->SetTransparency(transparency + 23);

            pMainDetectorVolume->AddNode(pGapVol, 0, new TGeoCombiTrans(trans, rot));
            continue;
        }

        ConcentricGap *pConcentricGap = NULL;
        pConcentricGap = dynamic_cast<ConcentricGap *>(*iter);

        if (NULL != pConcentricGap)
        {
            const double zMin = pConcentricGap->GetMinZCoordinate() * m_scalingFactor;
            const double zMax = pConcentricGap->GetMaxZCoordinate() * m_scalingFactor;
            const double zThick = zMax - zMin;

            TGeoVolume *pGapVol = MakePolygonTube(gapName.c_str(), pConcentricGap->GetInnerSymmetryOrder(), pConcentricGap->GetOuterSymmetryOrder(),
                pConcentricGap->GetInnerRCoordinate() * m_scalingFactor, pConcentricGap->GetOuterRCoordinate() * m_scalingFactor,
                pConcentricGap->GetInnerPhiCoordinate(), pConcentricGap->GetOuterPhiCoordinate(), (zThick / 2.), pGapMedium);

            pGapVol->SetLineColor(1);
            pGapVol->SetFillColor(1);
            pGapVol->SetTransparency(transparency + 23);

            pMainDetectorVolume->AddNode(pGapVol, 0, new TGeoTranslation(0, 0, zMin + (zThick / 2.f)));
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template void PandoraMonitoring::SetTreeVariable(const std::string &, const std::string &, float);
template void PandoraMonitoring::SetTreeVariable(const std::string &, const std::string &, int);
template void PandoraMonitoring::SetTreeVariable(const std::string &, const std::string &, double);

template void PandoraMonitoring::SetTreeVariable(const std::string &, const std::string &, std::vector<float> *);
template void PandoraMonitoring::SetTreeVariable(const std::string &, const std::string &, std::vector<int> *);
template void PandoraMonitoring::SetTreeVariable(const std::string &, const std::string &, std::vector<double> *);

} // namespace pandora_monitoring
