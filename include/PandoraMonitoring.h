/**
 *  @file   PandoraMonitoring/include/PandoraMonitoring.h
 *
 *  @brief  Header file for the pandora monitoring class.
 *
 *  $Log: $
 */
#ifndef PANDORA_MONITORING_H
#define PANDORA_MONITORING_H 1

#include "TApplication.h"
#include "TColor.h"

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/PandoraInternal.h"

#include "PandoraMonitoringApi.h"
#include "TTreeWrapper.h"

class TH1;
class TH2F;

class TEveElement;
class TGeoShape;
class TGeoVolume;
class TGeoMedium;

namespace pandora{class CartesianVector;}

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
     *  @brief  Get the pandora monitoring singleton
     */
    static PandoraMonitoring *GetInstance();

    /**
     *  @brief  Create a 1D histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  title the histogram title
     *  @param  nBinsX the number of x bins
     *  @param  xLow the the lower bound for the x axis
     *  @param  xUp the upper bound for the x axis
     */
    void Create1DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp,
        const std::string &xAxisTitle, const std::string &yAxisTitle);

    /**
     *  @brief  Create a 2D histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  title the histogram title
     *  @param  nBinsX the number of x bins
     *  @param  xLow the the lower bound for the x axis
     *  @param  xUp the upper bound for the x axis
     *  @param  nBinsY the number of y bins
     *  @param  yLow the the lower bound for the y axis
     *  @param  yUp the upper bound for the y axis
     */
    void Create2DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp, int nBinsY,
        double yLow, double yUp, const std::string &xAxisTitle, const std::string &yAxisTitle);

    /**
     *  @brief  Add a single entry to a 1D histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  xValue the x value to add
     *  @param  weight the weight to apply to the entry
     */
    void Fill1DHistogram(const std::string &name, float xValue, float weight);

    /**
     *  @brief  Add a single entry to a 2D histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  xValue the x value to add
     *  @param  yValue the y value to add
     *  @param  weight the weight to apply to the entry
     */
    void Fill2DHistogram(const std::string &name, float xValue, float yValue, float weight);

    /**
     *  @brief  Add, multiply or divide two histograms
     * 
     *  @param  nameHisto0 the name associated with the first histogram
     *  @param  nameHisto1 the name associated with the second histogram
     *  @param  coeff0 coefficient for the first histogram
     *  @param  coeff1 coefficient for the second histogram
     *  @param  add add histograms
     *  @param  multiply if "add" is false, multiply the histograms (if multiply is true), else divide
     */
    void AddMultiplyOrDivideHistograms(const std::string &nameHisto0, const std::string &nameHisto1, double coeff0, double coeff1,
        bool add, bool multiply );

    /**
     *  @brief  Set a variable in a tree (create the tree and the branch if not yet existing)
     * 
     *  @param  treeName name of the tree (is created if it does not exist yet)
     *  @param  variableName name of the branch in the tree (the branch is created if it does not exist yet)
     *  @param  variable sets value of the variable (permitted types are float/double/int and std::vector<float>*,std::vector<double>*,std::vector<int>*
     */
    template <typename VariableType>
    void SetTreeVariable(const std::string &treeName, const std::string &variableName, VariableType  variable);

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
     *  @brief  Draw a histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  options the drawing options
     */
    void DrawHistogram(const std::string &name, const std::string &options) const;

    /**
     *  @brief  Draw a pandora histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  options the drawing options
     */
    template <typename T>
    void DrawPandoraHistogram(const T &t, const std::string &options);

    /**
     *  @brief  Save a histogram to a file and tidy up
     * 
     *  @param  name the name associated with the histogram
     *  @param  fileName the file name under which to save the histogram
     *  @param  fileOptions the options associated with opening/recreating a file
     */
    void SaveAndCloseHistogram(const std::string &name, const std::string &fileName, const std::string &fileOptions);

    /**
     *  @brief  Delete a histogram
     * 
     *  @param  name the name associated with the histogram
     */
    void DeleteHistogram(const std::string &name);

    /**
     *  @brief Show the Eve Event-display and pause.
     */
    void ViewEvent();

    /**
     *  @brief  Set TEve display parameters
     * 
     *  @param  blackBackground whether to use a black background color, rather than white
     *  @param  showDetectors turns the visibility of the detector geometry on or off
     *  @param  transparencyThresholdE cell energy for which transparency is saturated (0%, fully opaque)
     *  @param  energyScaleThresholdE cell energy for which color is at top end of continous color palette
     */
    void SetEveDisplayParameters(const bool blackBackground, const bool showDetectors, const float transparencyThresholdE,
        const float energyScaleThresholdE);

    /**
     *  @brief Add MC particles to the Eve event-display
     * 
     *  @param pMCParticleList list of tracks to be added to the event display
     *  @param name of the MC particle list
     *  @param parent pointer to the parent TEveElement. If NULL, the track will be parent element
     *  @param color The color the track elements are drawn with
     *  @param pParticleSuppressionMap map from pdg-codes to energy for suppression of particles types below specific energies
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeMCParticles(const pandora::MCParticleList *const pMCParticleList, std::string name, TEveElement *parent, Color color, 
        const PandoraMonitoringApi::PdgCodeToEnergyMap *pParticleSuppressionMap);

    /**
     *  @brief Add Tracks to the Eve event-display
     * 
     *  @param pTrackList list of tracks to be added to the event display
     *  @param name of the track list
     *  @param parent pointer to the parent TEveElement. If NULL, the track will be parent element
     *  @param color The color the track elements are drawn with
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeTracks(const pandora::TrackList *const pTrackList, std::string name, TEveElement* parent, Color color);

    /**
     *  @brief Add CaloHits to the Eve event-display
     * 
     *  @param pCaloHitList list of calohits to be added to the event display
     *  @param parent name of the calohitlist
     *  @param parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param color The color the cluster elements are drawn with
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeCaloHits(const pandora::CaloHitList *const pCaloHitList, std::string name, TEveElement* parent, Color color, int pfoId = 0);

    /**
     *  @brief Add Particle flow objects to the Eve event-display
     * 
     *  @param pPfoList list of particle flow objects to be added to the event display
     *  @param name of the pfo list
     *  @param parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param color The color the cluster elements are drawn with
     *  @param showAssociatedTracks draw the tracks associated to the cluster
     */
    TEveElement *VisualizeParticleFlowObjects(const pandora::PfoList *const pPfoList, std::string name, TEveElement* parent,
        Color color, bool showAssociatedTracks);

    /**
     *  @brief Add Clusters to the Eve event-display
     * 
     *  @param pClusterList list of clusters to be added to the event display
     *  @param name of the cluster list
     *  @param parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param color The color the cluster elements are drawn with
     *  @param showAssociatedTracks draw the tracks associated to the cluster
     *
     *  @return pointer to created TEveElement
     */
    TEveElement *VisualizeClusters(const pandora::ClusterList *const pClusterList, std::string name, TEveElement* parent, Color color,
        bool showAssociatedTracks, int pfoId = 0);

    /**
     *  @brief  Pause until user enters 'return'
     */
    void Pause() const;

    /**
     *  @brief  Destructor
     */
    ~PandoraMonitoring();

    /**
     *  @brief Delete instance of PandoraMonitoring
     */
    void DeleteInstance();

private:
    typedef std::vector< std::pair< double,double > > DoublePairVector;

    /**
     *  @brief  Default constructor
     */
    PandoraMonitoring();

    /**
     *  @brief Computes the corners of a box in 3D
     * 
     *  @param pCaloHit address of the calo hit
     *  @param corners will be filled with the x,y and z-coordinates of all 8 corners of the box
     */
    void MakeCaloHitCell(const pandora::CaloHit *const pCaloHit, float corners[24]);

    /**
     *  @brief compute the polygon corners for the detector outline
     * 
     *  @param symmetryOrder is the number of polygon corners
     *  @param closestDistanceToIp is the distance to the closest points on the polygon
     *  @param phi0 reference angle where to start the polygon
     *  @param coordinates vector of double,double pairs which is filled with the x and y coordinates of the polygon corners
     */
    void ComputePolygonCorners(int symmetryOrder, double closestDistanceToIp, double phi0, std::vector<std::pair<double,double> > &coordinates);

    /**
     *  @brief Creates a "tube" volume with the given symmetry inside and outside. If a symmetryOrder <= 2 is chosen, a circle is
     *         used instead of a polygon
     * 
     *  @param name of the volume
     *  @param innerSymmetryOrder symmetry order of the inner polygon (circle)
     *  @param outerSymmetryOrder symmetry order of the outer polygon (circle)
     *  @param innerClosestDistanceToIp closest distance between IP and polygon for the inner part of the tube
     *  @param outerClosestDistanceToIp closest distance between IP and polygon for the outer part of the tube
     *  @param innerPhi0 starting angle of the inner polygon
     *  @param outerPhi0 starting angle of the outer polygon
     *  @param halfLength half length (z) of the tube
     *  @param medium TGeoMedium of the volume
     */
    TGeoVolume *MakePolygonTube(std::string name, int innerSymmetryOrder, int outerSymmetryOrder, double innerClosestDistanceToIp,
        double outerClosestDistanceToIp, double innerPhi0, double outerPhi0, double halfLength, TGeoMedium *medium = 0);

    /**
     *  @brief Creates a extruded polygonal (or cylindrical) shape with the given symmetry. If a symmetryOrder <= 2 is chosen, a
     *         circle is used instead of a polygon
     * 
     *  @param symmetryOrder symmetry order of the polygon (circle if <=2)
     *  @param closestDistanceToIp closest distance between IP and polygon (circle radius)
     *  @param phi starting angle of the polygon
     *  @param halfLength half length (z) of the tube
     */
    TGeoShape *MakePolygonTube(int symmetryOrder, double closestDistanceToIp, double phi, double halfLength);

    /**
     *  @brief  Transform a Pandora monitoring API color enum into a ROOT color enum
     * 
     *  @param  color in Pandora monitoring API enum
     */
    EColor GetROOTColor(Color color);

    /**
     *  @brief  Get a color for a PDG code
     * 
     *  @param  particleId pdgCode of the particle
     */
    Color GetColorForPdgCode(int particleId);

    /**
     *  @brief  Initialize eve elements
     * 
     *  @param  transparency the transparency
     */
    void InitializeEve(Char_t transparency = 70);

    /**
     *  @brief  Initialize subdetector eve elements
     * 
     *  @param  pMainDetectorVolume address of the main detector volume
     *  @param  pSubDetectorMedium address of the medium to be used for the subdetectors
     *  @param  transparency the transparency
     */
    void InitializeSubDetectors(TGeoVolume *pMainDetectorVolume, TGeoMedium *pSubDetectorMedium, Char_t transparency);

    /**
     *  @brief  Initialize detector gap eve elements
     * 
     *  @param  pMainDetectorVolume address of the main detector volume
     *  @param  pSubDetectorMedium address of the medium to be used for the gaps
     *  @param  transparency the transparency
     */
    void InitializeGaps(TGeoVolume *pMainDetectorVolume, TGeoMedium *pGapMedium, Char_t transparency);

    typedef std::map<const std::string, TH1 *> HistogramMap;

    static bool                 m_instanceFlag;             ///< The pandora monitoring instance flag
    static PandoraMonitoring   *m_pPandoraMonitoring;       ///< The pandora monitoring instance
    TApplication               *m_pApplication;             ///< The root application

    HistogramMap                m_histogramMap;             ///< The histogram map
    TTreeWrapper                m_treeWrapper;              ///< wrapper around TTree functionality

    static bool                 m_eveInitialized;           ///< is set if ROOT Eve is initialized
    static float                m_scalingFactor;            ///< TEve works with [cm], Pandora works with [mm]
    static bool                 m_openEveEvent;             ///< is set if an Event is open to store objects (hits, clusters,...) in it.
    static int                  m_eventDisplayCounter;      ///< counter for the event displays

    float                       m_transparencyThresholdE;   ///< Cell energy for which transparency is saturated (0%, fully opaque)
    float                       m_energyScaleThresholdE;    ///< Cell energy for which color is at top end of continous color palette
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline PandoraMonitoring::PandoraMonitoring() :
    m_transparencyThresholdE(-1.f),
    m_energyScaleThresholdE(-1.f)
{
    int argc = 0;
    char* argv = (char *)"";

    m_pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
    m_pApplication->SetReturnFromRun(kTRUE);
}

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

inline EColor PandoraMonitoring::GetROOTColor(Color color)
{
    switch (color)
    {
        COLOR_TABLE(GET_ROOT_COLOR)
        default : return kGray;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Color PandoraMonitoring::GetColorForPdgCode(int particleId)
{
    switch (particleId)
    {
        PARTICLE_DATA_COLOR_TABLE(GET_PDG_COLOR_SWITCH)
        default : return GRAY;
    }
}

} // namespace pandora_monitoring

#endif // #ifndef PANDORA_MONITORING_H
