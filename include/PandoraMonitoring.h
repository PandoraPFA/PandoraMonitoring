/**
 *  @file   PandoraMonitoring/include/PandoraMonitoring.h
 *
 *  @brief  Header file for the pandora monitoring class.
 *
 *  $Log: $
 */
#ifndef PANDORA_MONITORING_H
#define PANDORA_MONITORING_H 1

#include "TColor.h"
#include "TGLViewer.h"

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/PandoraInternal.h"

#include "PandoraMonitoringApi.h"
#include "TTreeWrapper.h"

class TApplication;
class TEveElement;
class TEveManager;
class TEveScene;
class TEveViewer;
class TGeoShape;
class TGeoVolume;
class TGeoMedium;

namespace pandora { class CartesianVector; class LineGap; class BoxGap; class ConcentricGap; class Pandora; }

//------------------------------------------------------------------------------------------------------------------------------------------

namespace pandora_monitoring
{

/**
 *  @brief  PandoraMonitoring singleton class
 */
class PandoraMonitoring
{
public:
    /**
     *  @brief  Get the relevant pandora monitoring instance
     * 
     *  @param  pandora the calling pandora instance
     */
    static PandoraMonitoring *GetInstance(const pandora::Pandora &pandora);

    /**
     *  @brief Delete instance of PandoraMonitoring
     * 
     *  @param  pandora the calling pandora instance
     */
    static void DeleteInstance(const pandora::Pandora &pandora);

    /**
     *  @brief  Set a variable in a tree (create the tree and the branch if not yet existing)
     * 
     *  @param  treeName name of the tree (is created if it does not exist yet)
     *  @param  variableName name of the branch in the tree (the branch is created if it does not exist yet)
     *  @param  t sets value of the variable (permitted types are float/double/int and std::vector<float>*,std::vector<double>*,std::vector<int>*
     */
    template <typename T>
    void SetTreeVariable(const std::string &treeName, const std::string &variableName, T t);

    /**
     *  @brief  Fill the tree with the variables which have been set before with SetTreeVariable
     * 
     *  @param  treeName name of the tree to be filled
     */
    void FillTree(const std::string &treeName);

    /**
     *  @brief  Print the tree
     * 
     *  @param  treeName name of the tree to be printed
     */
    void PrintTree(const std::string &treeName);

    /**
     *  @brief  Scan the tree (print the values of all branches)
     * 
     *  @param  treeName name of the tree to be scanned
     */
    void ScanTree(const std::string &treeName);

    /**
     *  @brief  Save the tree to a file
     *  @param  fileName the file name under which to save the histogram
     *  @param  fileOptions the options associated with opening/recreating a file
     * 
     *  @param  treeName name of the tree to be written to a file
     */
    void SaveTree(const std::string &treeName, const std::string &fileName, const std::string &fileOptions);

    /**
     *  @brief  Draw a pandora histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  options the drawing options
     */
    template <typename T>
    void DrawPandoraHistogram(const T &t, const std::string &options);

    /**
     *  @brief  Set TEve display parameters
     * 
     *  @param  showDetectors turns the visibility of the detector geometry on or off
     *  @param  detectorView the detector view
     *  @param  transparencyThresholdE cell energy for which transparency is saturated (0%, fully opaque)
     *  @param  energyScaleThresholdE cell energy for which color is at top end of continous color palette
     *  @param  scalingFactor TEve works with [cm], Pandora works with [mm] (unless user has decided to use alternative units)
     */
    void SetEveDisplayParameters(const bool showDetectors, const DetectorView detectorView, const float transparencyThresholdE,
        const float energyScaleThresholdE, const float scalingFactor);

    /**
     *  @brief  Add MC particles to the Eve event-display
     * 
     *  @param  pMCParticleList list of tracks to be added to the event display
     *  @param  name of the MC particle list
     *  @param  parent pointer to the parent TEveElement. If NULL, the track will be parent element
     *  @param  color The color the track elements are drawn with
     *  @param  pParticleSuppressionMap map from pdg-codes to energy for suppression of particles types below specific energies
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeMCParticles(const pandora::MCParticleList *const pMCParticleList, const std::string &name, TEveElement *parent,
        const Color color, const PandoraMonitoringApi::PdgCodeToEnergyMap *pParticleSuppressionMap);

    /**
     *  @brief  Add Tracks to the Eve event-display
     * 
     *  @param  pTrackList list of tracks to be added to the event display
     *  @param  name of the track list
     *  @param  parent pointer to the parent TEveElement. If NULL, the track will be parent element
     *  @param  color The color the track elements are drawn with
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeTracks(const pandora::TrackList *const pTrackList, const std::string &name, TEveElement *parent,
        const Color color);

    /**
     *  @brief  Add CaloHits to the Eve event-display
     * 
     *  @param  pCaloHitList list of calohits to be added to the event display
     *  @param  parent name of the calohitlist
     *  @param  parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param  color The color the cluster elements are drawn with
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeCaloHits(const pandora::CaloHitList *const pCaloHitList, const std::string &name, TEveElement *parent,
        const Color color);

    /**
     *  @brief  Add Clusters to the Eve event-display
     * 
     *  @param  pClusterList list of clusters to be added to the event display
     *  @param  name of the cluster list
     *  @param  parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param  color The color the cluster elements are drawn with
     *  @param  showAssociatedTracks draw the tracks associated to the cluster
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeClusters(const pandora::ClusterList *const pClusterList, const std::string &name, TEveElement *parent,
        const Color color, bool showAssociatedTracks);

    /**
     *  @brief  Add Particle flow objects to the Eve event-display
     * 
     *  @param  pPfoList list of particle flow objects to be added to the event display
     *  @param  name of the pfo list
     *  @param  parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param  color The color the cluster elements are drawn with
     *  @param  showAssociatedTracks draw the tracks associated to the cluster
     *  @param  showVertices whether to draw markers to represent the pfo vertices
     *  @param  displayPfoHierarchy whether to draw daughter particles (only) within parent elements
     */
    TEveElement *VisualizeParticleFlowObjects(const pandora::PfoList *const pPfoList, const std::string &name, TEveElement *parent,
        const Color color, bool showVertices, bool displayPfoHierarchy);

    /**
     *  @brief  Add vertices to the Eve event-display
     * 
     *  @param  pVertexList list of vertices to be added to the event display
     *  @param  parent name of the vertexlist
     *  @param  parent pointer to the parent TEveElement. If NULL, the vertex will be parent element
     *  @param  color The color the vertex elements are drawn with
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeVertices(const pandora::VertexList *const pVertexList, const std::string &name, TEveElement *parent,
        const Color color);

    /**
     *  @brief  Add marker to visualization
     * 
     *  @param  pMarkerPoint address of the marker point
     *  @param  name of the marker
     *  @param  parent pointer to the parent TEveElement. If NULL, the marker will be parent element
     *  @param  color the marker color
     *  @param  markerSize the marker size
     */
    TEveElement *AddMarkerToVisualization(const pandora::CartesianVector *const pMarkerPoint, const std::string &name, TEveElement *parent,
        const Color color, const unsigned int markerSize);

    /**
     *  @brief  Add line to visualization
     * 
     *  @param  pStartPoint, pEndPoint address start and end of line
     *  @param  name of the line
     *  @param  parent pointer to the parent TEveElement. If NULL, the marker will be parent element
     *  @param  color the line color
     *  @param  lineWidth the line width
     *  @param  lineStyle the line style
     */
    TEveElement *AddLineToVisualization(const pandora::CartesianVector *const pStartPoint, const pandora::CartesianVector *const pEndPoint,
        const std::string &name, TEveElement *parent, const Color color, const unsigned int lineWidth, const unsigned int lineStyle);

    /**
     *  @brief  Show the Eve Event-display and pause.
     */
    void ViewEvent();

    /**
     *  @brief  Pause thread until user enters 'return'
     */
    void Pause() const;

    /**
     *  @brief  Sort clusters by descending hadronic energy
     * 
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortClustersByHadronicEnergy(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort MCParticles by descending energy
     * 
     *  @param  pLhs address of first MCParticle
     *  @param  pRhs address of second MCParticle
     */
    static bool SortMCParticlesByEnergy(const pandora::MCParticle *const pLhs, const pandora::MCParticle *const pRhs);

    /**
     *  @brief  Sort tracks by descending momentum
     * 
     *  @param  pLhs address of first track
     *  @param  pRhs address of second track
     */
    static bool SortTracksByMomentum(const pandora::Track *const pLhs, const pandora::Track *const pRhs);

    /**
     *  @brief  Sort pfos by descending energy 
     * 
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortPfosByEnergy(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

private:
    /**
     *  @brief  Default constructor
     * 
     *  @param  pandora the associated pandora instance
     */
    PandoraMonitoring(const pandora::Pandora &pandora);

    /**
     *  @brief  Destructor
     */
    ~PandoraMonitoring();

    /**
     *  @brief  Computes the corners of a box in 3D
     * 
     *  @param  pCaloHit address of the calo hit
     *  @param  corners will be filled with the x,y and z-coordinates of all 8 corners of the box
     */
    void MakeCaloHitCell(const pandora::CaloHit *const pCaloHit, float corners[24]);

    typedef std::vector< std::pair<double, double> > DoublePairVector;

    /**
     *  @brief  Compute the polygon corners for the detector outline
     * 
     *  @param  symmetryOrder is the number of polygon corners
     *  @param  closestDistanceToIp is the distance to the closest points on the polygon
     *  @param  phi0 reference angle where to start the polygon
     *  @param  coordinates vector of double,double pairs which is filled with the x and y coordinates of the polygon corners
     */
    void ComputePolygonCorners(int symmetryOrder, double closestDistanceToIp, double phi0, DoublePairVector &coordinates);

    /**
     *  @brief  Creates a "tube" volume with the given symmetry inside and outside. If a symmetryOrder <= 2 is chosen, a circle is
     *          used instead of a polygon
     * 
     *  @param  name of the volume
     *  @param  innerSymmetryOrder symmetry order of the inner polygon (circle)
     *  @param  outerSymmetryOrder symmetry order of the outer polygon (circle)
     *  @param  innerClosestDistanceToIp closest distance between IP and polygon for the inner part of the tube
     *  @param  outerClosestDistanceToIp closest distance between IP and polygon for the outer part of the tube
     *  @param  innerPhi0 starting angle of the inner polygon
     *  @param  outerPhi0 starting angle of the outer polygon
     *  @param  halfLength half length (z) of the tube
     *  @param  medium TGeoMedium of the volume
     */
    TGeoVolume *MakePolygonTube(std::string name, int innerSymmetryOrder, int outerSymmetryOrder, double innerClosestDistanceToIp,
        double outerClosestDistanceToIp, double innerPhi0, double outerPhi0, double halfLength, TGeoMedium *medium = NULL);

    /**
     *  @brief  Creates a extruded polygonal (or cylindrical) shape with the given symmetry. If a symmetryOrder <= 2 is chosen, a
     *          circle is used instead of a polygon
     * 
     *  @param  symmetryOrder symmetry order of the polygon (circle if <=2)
     *  @param  closestDistanceToIp closest distance between IP and polygon (circle radius)
     *  @param  phi starting angle of the polygon
     *  @param  halfLength half length (z) of the tube
     */
    TGeoShape *MakePolygonTube(int symmetryOrder, double closestDistanceToIp, double phi, double halfLength);

    /**
     *  @brief  Get a string representing a member of the pandora hit type enummeration
     *
     *  @return the hit type as a string
     */
    std::string GetHitTypeString(const pandora::HitType hitType) const;

    /**
     *  @brief  Transform a Pandora monitoring API color enum into a ROOT color enum
     * 
     *  @param  color in Pandora monitoring API enum
     */
    EColor GetROOTColor(Color color) const;

    /**
     *  @brief  Get a color for a PDG code
     * 
     *  @param  particleId pdgCode of the particle
     */
    Color GetColorForPdgCode(int particleId) const;

    /**
     *  @brief  Initialize eve elements
     * 
     *  @param  transparency the transparency
     */
    void InitializeEve(Char_t transparency = 70);

    /**
     *  @brief  Initialize the viewers
     */
    void InitializeViews();

    /**
     *  @brief  Add event and geometry scences to a specifed view and specify camera and axis type for view
     *
     *  @param  pTEveViewer viewer 
     *  @param  pTEveEventScene event scene 
     *  @param  pTEveGeometryScene geometry scene
     *  @param  camera orthogonal/perspective view choice
     *  @param  axisType for view
     */
    void AddScenes(TEveViewer *pTEveViewer, TEveScene *pTEveEventScene, TEveScene *pTEveGeometryScene, TGLViewer::ECameraType camera, int axisType);

    /**
     *  @brief  Initialize subdetector eve elements
     * 
     *  @param  pMainDetectorVolume address of the main detector volume
     *  @param  pSubDetectorMedium address of the medium to be used for the subdetectors
     *  @param  transparency the transparency
     */
    void InitializeSubDetectors(TGeoVolume *pMainDetectorVolume, TGeoMedium *pSubDetectorMedium, Char_t transparency);

    /**
     *  @brief  Initialize subdetector eve elements
     * 
     *  @param  pMainDetectorVolume address of the main detector volume
     *  @param  pSubDetectorMedium address of the medium to be used for the subdetectors
     *  @param  transparency the transparency
     */
    void InitializeLArTPCs(TGeoVolume *pMainDetectorVolume, TGeoMedium *pSubDetectorMedium, Char_t transparency);

    /**
     *  @brief  Initialize detector gap eve elements
     * 
     *  @param  pMainDetectorVolume address of the main detector volume
     *  @param  pSubDetectorMedium address of the medium to be used for the gaps
     *  @param  transparency the transparency
     */
    void InitializeGaps(TGeoVolume *pMainDetectorVolume, TGeoMedium *pGapMedium, Char_t transparency);

    /**
     *  @brief  Whether lhs object should appear before rhs object in a sorted list
     *
     *  @param  pLhs first object for comparison
     *  @param  pRhs second object for comparison
     * 
     *  @return boolean
     */
    static bool SortLineGaps(const pandora::LineGap *const pLhs, const pandora::LineGap *const pRhs);

    /**
     *  @brief  Whether lhs object should appear before rhs object in a sorted list
     *
     *  @param  pLhs first object for comparison
     *  @param  pRhs second object for comparison
     * 
     *  @return boolean
     */
    static bool SortBoxGaps(const pandora::BoxGap *const pLhs, const pandora::BoxGap *const pRhs);

    /**
     *  @brief  Whether lhs object should appear before rhs object in a sorted list
     *
     *  @param  pLhs first object for comparison
     *  @param  pRhs second object for comparison
     * 
     *  @return boolean
     */
    static bool SortConcentricGaps(const pandora::ConcentricGap *const pLhs, const pandora::ConcentricGap *const pRhs);

    typedef std::vector<const pandora::LineGap*> LineGapVector;
    typedef std::vector<const pandora::BoxGap*> BoxGapVector;
    typedef std::vector<const pandora::ConcentricGap*> ConcentricGapVector;

    typedef std::unordered_map<const pandora::Pandora*, PandoraMonitoring*> MonitoringInstanceMap;
    static MonitoringInstanceMap    m_monitoringInstanceMap;    ///< The monitoring instance map

    const pandora::Pandora *const   m_pPandora;                 ///< The associated pandora instance
    TApplication                   *m_pApplication;             ///< The root application
    TEveManager                    *m_pEveManager;              ///< The eve manager
    TTreeWrapper                   *m_pTreeWrapper;             ///< wrapper around TTree functionality

    TEveViewer                     *m_p3DView;                  ///< Viewer for the 3D view
    TEveViewer                     *m_p2DUView;                 ///< Viewer for the U 2D view
    TEveViewer                     *m_p2DVView;                 ///< Viewer for the V 2D view
    TEveViewer                     *m_p2DWView;                 ///< Viewer for the W 2D view

    TEveScene                      *m_p3DEventScene;            ///< Scene containing the 3D event
    TEveScene                      *m_p2DUEventScene;           ///< Scene containing the U 2D event
    TEveScene                      *m_p2DVEventScene;           ///< Scene containing the V 2D event
    TEveScene                      *m_p2DWEventScene;           ///< Scene containing the W 2D event

    TEveScene                      *m_p3DGeometryScene;         ///< Scene containing the 3D geometry
    TEveScene                      *m_p2DUGeometryScene;        ///< Scene containing the U 2D geometry
    TEveScene                      *m_p2DVGeometryScene;        ///< Scene containing the V 2D geometry
    TEveScene                      *m_p2DWGeometryScene;        ///< Scene containing the W 2D geometry

    float                           m_scalingFactor;            ///< TEve works with [cm], Pandora works with [mm]
    bool                            m_openEveEvent;             ///< is set if an Event is open to store objects (hits, clusters,...) in it.
    int                             m_eventDisplayCounter;      ///< counter for the event displays
    double                          m_minXLArTPC;               ///< Minimum x position to which any LAr TPC in the event extends
    double                          m_maxXLArTPC;               ///< Maximum x position to which any LAr TPC in the event extends
    double                          m_minYLArTPC;               ///< Minimum y position to which any LAr TPC in the event extends
    double                          m_maxYLArTPC;               ///< Maximum y position to which any LAr TPC in the event extends
    double                          m_minZLArTPC;               ///< Minimum z position to which any LAr TPC in the event extends
    double                          m_maxZLArTPC;               ///< Maximum z position to which any LAr TPC in the event extends
    float                           m_transparencyThresholdE;   ///< Cell energy for which transparency is saturated (0%, fully opaque)
    float                           m_energyScaleThresholdE;    ///< Cell energy for which color is at top end of continous color palette
    bool                            m_showDetectors;            ///< Turns the visibility of the detector geometry on or off
    DetectorView                    m_detectorView;             ///< The detector view
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

// Specify (color, ROOT-color)
#define COLOR_TABLE(d)                                                      \
    d(WHITE,           kWhite)                                              \
    d(BLACK,           kBlack)                                              \
    d(RED,             kRed)                                                \
    d(GREEN,           kGreen)                                              \
    d(BLUE,            kBlue)                                               \
    d(MAGENTA,         kMagenta)                                            \
    d(CYAN,            kCyan)                                               \
    d(VIOLET,          kViolet)                                             \
    d(PINK,            kPink)                                               \
    d(ORANGE,          kOrange)                                             \
    d(YELLOW,          kYellow)                                             \
    d(SPRING,          kSpring)                                             \
    d(TEAL,            kTeal)                                               \
    d(AZURE,           kAzure)                                              \
    d(GRAY,            kGray)                                               \
    d(DARKRED,         EColor(TColor::GetColorDark(kRed)))                  \
    d(DARKGREEN,       EColor(TColor::GetColorDark(kGreen)))                \
    d(DARKBLUE,        EColor(TColor::GetColorDark(kBlue)))                 \
    d(DARKMAGENTA,     EColor(TColor::GetColorDark(kMagenta)))              \
    d(DARKCYAN,        EColor(TColor::GetColorDark(kCyan)))                 \
    d(DARKVIOLET,      EColor(TColor::GetColorDark(kViolet)))               \
    d(DARKPINK,        EColor(TColor::GetColorDark(kPink)))                 \
    d(DARKORANGE,      EColor(TColor::GetColorDark(kOrange)))               \
    d(DARKYELLOW,      EColor(TColor::GetColorDark(kYellow)))               \
    d(LIGHTGREEN,      EColor(TColor::GetColorBright(kGreen)))              \
    d(LIGHTBLUE,       EColor(TColor::GetColorBright(kBlue)))               \
    d(LIGHTRED,        EColor(TColor::GetColorBright(kRed)))                \
    d(LIGHTMAGENTA,    EColor(TColor::GetColorBright(kMagenta)))            \
    d(LIGHTCYAN,       EColor(TColor::GetColorBright(kCyan)))               \
    d(LIGHTVIOLET,     EColor(TColor::GetColorBright(kViolet)))             \
    d(LIGHTPINK,       EColor(TColor::GetColorBright(kPink)))               \
    d(LIGHTORANGE,     EColor(TColor::GetColorBright(kOrange)))             \
    d(LIGHTYELLOW,     EColor(TColor::GetColorBright(kYellow)))

//------------------------------------------------------------------------------------------------------------------------------------------

// Specify (name, pdg code, color)
#define PARTICLE_DATA_COLOR_TABLE(d)                                        \
    d(PHOTON,               22,             DARKYELLOW)                     \
    d(E_MINUS,              11,             LIGHTBLUE)                      \
    d(E_PLUS,              -11,             LIGHTRED)                       \
    d(MU_MINUS,             13,             BLUE)                           \
    d(MU_PLUS,             -13,             RED)                            \
    d(TAU_MINUS,            15,             DARKBLUE)                       \
    d(TAU_PLUS,            -15,             DARKRED)                        \
    d(NU_E,                 12,             DARKBLUE)                       \
    d(NU_E_BAR,            -12,             DARKRED)                        \
    d(NU_MU,                14,             DARKBLUE)                       \
    d(NU_MU_BAR,           -14,             DARKRED)                        \
    d(NU_TAU,               16,             DARKBLUE)                       \
    d(NU_TAU_BAR,          -16,             DARKRED)                        \
    d(PI_PLUS,             211,             MAGENTA)                        \
    d(PI_MINUS,           -211,             VIOLET)                         \
    d(PI_ZERO,             111,             LIGHTGREEN)                     \
    d(LAMBDA,             3122,             DARKGREEN)                      \
    d(LAMBDA_BAR,        -3122,             DARKGREEN)                      \
    d(K_PLUS,              321,             DARKGREEN)                      \
    d(K_MINUS,            -321,             DARKGREEN)                      \
    d(K_SHORT,             310,             LIGHTGREEN)                     \
    d(K_LONG,              130,             GREEN)                          \
    d(SIGMA_MINUS,        3112,             GREEN)                          \
    d(SIGMA_PLUS,         3222,             GREEN)                          \
    d(PROTON,             2212,             ORANGE)                         \
    d(NEUTRON,            2112,             CYAN)

/**
 *  @brief  The mass switch statement macro
 */
#define GET_ROOT_COLOR(a, b)                                                \
    case a : return b;

/**
 *  @brief  The mass switch statement macro
 */
#define GET_PARTICLE_COLOR_SWITCH(a, b, c)                                  \
    case a : return c;

/**
 *  @brief  The mass switch statement macro
 */
#define GET_PDG_COLOR_SWITCH(a, b, c)                                       \
    case b : return c;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline EColor PandoraMonitoring::GetROOTColor(Color color) const
{
    switch (color)
    {
        COLOR_TABLE(GET_ROOT_COLOR)
        default : return kGray;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Color PandoraMonitoring::GetColorForPdgCode(int particleId) const
{
    switch (particleId)
    {
        PARTICLE_DATA_COLOR_TABLE(GET_PDG_COLOR_SWITCH)
        default : return GRAY;
    }
}

} // namespace pandora_monitoring

#endif // #ifndef PANDORA_MONITORING_H
