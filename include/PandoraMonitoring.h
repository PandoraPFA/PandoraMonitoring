/**
 *  @file   PandoraMonitoring/include/PandoraMonitoring.h
 *
 *  @brief  Header file for the pandora monitoring class.
 *
 *  $Log: $
 */
#ifndef PANDORA_MONITORING_H
#define PANDORA_MONITORING_H 1

#include "Pandora/PandoraInternal.h"

#include "PandoraMonitoringApi.h"

#include "TApplication.h"
#include "TColor.h"

#include "TTreeWrapper.h"


#include <iostream>
#include <map>

class TH1;
class TH2F;
class TArrow;
class TObject;
class TPolyMarker;
class TTree;
class TBranch;
class TCanvas;

class TEveElement;
class TGeoShape;
class TGeoVolume;
class TGeoMedium;

namespace pandora
{
    class CartesianVector;
} // namespace pandora


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
     *  @brief  Reset pandora monitoring, deleting the singleton
     */
    void Reset();

    /**
     *  @brief  Create a 1D histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  title the histogram title
     *  @param  nBinsX the number of x bins
     *  @param  xLow the the lower bound for the x axis
     *  @param  xUp the upper bound for the x axis
     */
    void Create1DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp);

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
        double yLow, double yUp);

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
     *  @brief  Add ClusterList to the output. The canvas is automatically created if not existing. Pause the output with ViewEvent()
     * 
     *  @param  detectorView the detector view
     *  @param  pClusterList address of the cluster list
     *  @param  color in which the clusters should be drawn
     */ 
    void AddClusterList(DetectorView detectorView, const pandora::ClusterList *const pClusterList, Color color = AUTO );

    /**
     *  @brief  Add TrackList to the output. The canvas is automatically created if not existing. Pause the output with ViewEvent()
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     *  @param  color in which the tracks should be drawn
     */ 
    void AddTrackList(DetectorView detectorView, const pandora::TrackList *const pTrackList, Color color = AUTO );

    /**
     *  @brief  Add OrderedCaloHitList to the output. The canvas is automatically created if not existing. Pause the output with ViewEvent()
     * 
     *  @param  detectorView the detector view
     *  @param  pOrderedCaloHitList address of the OrderedCaloHit list
     *  @param  color in which the OrderedCaloHits should be drawn
     */ 
    void AddCaloHitList(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList, Color color = AUTO );

    /**
     *  @brief  Pauses the monitoring, such that the user can see the output. Clears and deletes the canvases.
     * 
     */ 
    void ViewEvent();

    /**
     *  @brief Show the Eve Event-display and pause.
     * 
     */  
    void View();


    /**
     *  @brief Add Particle flow objects to the Eve event-display
     * 
     *  @param pPfoList list of particle flow objects to be added to the event display
     *  @param name of the pfo list
     *  @param parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param color The color the cluster elements are drawn with
     *  @param showAssociatedTracks draw the tracks associated to the cluster
     *  @param showFit draw an arrow representing the fit through the calorimeterhits (the fit is computed within pandora)
     */  
    TEveElement* VisualizeParticleFlowObjects(const pandora::ParticleFlowObjectList *const pPfoList, std::string name, TEveElement* parent, Color color, 
					      bool showAssociatedTracks, bool showFit  );

    /**
     *  @brief Add Clusters to the Eve event-display
     * 
     *  @param pClusterList list of clusters to be added to the event display
     *  @param name of the cluster list
     *  @param parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param color The color the cluster elements are drawn with
     *  @param showAssociatedTracks draw the tracks associated to the cluster
     *  @param showFit draw an arrow representing the fit through the calorimeterhits (the fit is computed within pandora)
     */  
    TEveElement* VisualizeClusters(const pandora::ClusterList *const pClusterList, std::string name, TEveElement* parent, Color color, 
				   bool showAssociatedTracks, bool showFit  );

    /**
     *  @brief Add Tracks to the Eve event-display
     * 
     *  @param pTrackList list of tracks to be added to the event display
     *  @param name of the track list
     *  @param parent pointer to the parent TEveElement. If NULL, the track will be parent element
     *  @param color The color the track elements are drawn with
     */  
    TEveElement* VisualizeTracks(const pandora::TrackList *const pTrackList, std::string name, TEveElement* parent, Color color );

    /**
     *  @brief Add CaloHits to the Eve event-display
     * 
     *  @param pOrderedCaloHitList list of calohits to be added to the event display
     *  @param parent name of the calohitlist
     *  @param parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param color The color the cluster elements are drawn with
     *  
     *  @return pointer to created TEveElement
     */  
    TEveElement* VisualizeCaloHits(const pandora::OrderedCaloHitList *const pOrderedCaloHitList, std::string name, TEveElement* parent, Color color );

    /**
     *  @brief  Temporary function - just draw the detector outline
     * 
     *  @param  detectorView the detector view
     */ 
    void DetectorOutlineTest(DetectorView detectorView);

    /**
     *  @brief  Temporary function - draw a test canvas and histogram
     */
    void DrawCanvas() const;

    /**
     *  @brief  Temporary function - display the parent addresses of calo hits contained in a cluster list
     * 
     *  @param  pClusterList address of the cluster list
     */
    void LookAtClusters(const pandora::ClusterList *const pClusterList) const;

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
    /**
     *  @brief  The x-y outline parameters class
     */
    class XYOutlineParameters
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  symmetryOrder the order of symmetry of the detector outline
         *  @param  phi0 the angle through which the detector outline polygon has been rotated wrt the cartesian y axis
         *  @param  closestDistanceToIp the closest distance of the detector outline to the interaction point
         */
        XYOutlineParameters(int symmetryOrder, float phi0, float closestDistanceToIp);

        /**
         *  @brief  Sort a list of x-y outline parameters, by decreasing radial distance from the interaction point
         */
        static bool Sort(const XYOutlineParameters &lhs, const XYOutlineParameters &rhs);

        /**
         *  @brief  Get the order of symmetry of the detector outline
         * 
         *  @return The order of symmetry of the detector outline
         */
        int GetSymmetryOrder() const;

        /**
         *  @brief  Get the angle through which the detector outline polygon has been rotated wrt the cartesian y axis
         * 
         *  @return The angle through which the detector outline polygon has been rotated wrt the cartesian y axis
         */
        float GetPhi0() const;

        /**
         *  @brief  Get the closest distance of the detector outline to the interaction point
         * 
         *  @return The closest distance of the detector outline to the interaction point
         */
        float GetClosestDistanceToIp() const;

    private:
        int     m_symmetryOrder;        ///< the order of symmetry of the detector outline
        float   m_phi0;                 ///< The angle through which the detector outline polygon has been rotated wrt the cartesian y axis
        float   m_closestDistanceToIp;  ///< The closest distance of the detector outline to the interaction point
    };

    typedef std::vector<XYOutlineParameters> XYOutlineParametersList;

    /**
     *  @brief  Default constructor
     */
    PandoraMonitoring();

    /**
     *  @brief  Draw the detector outline
     * 
     *  @param  detectorView the detector view
     */ 
    void DrawDetectorOutline(DetectorView detectorView);

    /**
     *  @brief  Draw clusters in an event
     *
     *  @param  detectorView the detector view
     *  @param  pClusterList address of the cluster list
     *  @param  color of the cluster
     */ 
    void DrawClusters(DetectorView detectorView, const pandora::ClusterList *const pClusterList, Color color = AUTO);

    /**
     *  @brief  Draw tracks in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     */ 
    void DrawTracks(DetectorView detectorView, const pandora::TrackList *const pTrackList, Color color = AUTO);

    /**
     *  @brief  Draw calo hits in an event
     *
     *  @param  detectorView the detector view
     *  @param  pOrderedCaloHitList address of the ordered calo hit list
     */ 
    void DrawCaloHits(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList, Color color = AUTO);

    /**
     *  @brief  Get a canvas with a particular detectorView, create and initialize it if not existing
     */
    TCanvas* GetCanvas(DetectorView detectorView);

    /**
     *  @brief  Construct the detector outline
     */
    void MakeDetectorOutline();

    /**
     *  @brief  Construct the xy outline of a detector layer
     * 
     *  @param  xyOutlineParameters x-y outline parameters
     */
    void MakeXYLayerOutline(const XYOutlineParameters &xyOutlineParameters);

    /**
     *  @brief  Construct the xz outline of a detector layer
     * 
     *  @param  innerRCoordinate the inner r coordinate
     *  @param  outerRCoordinate the outer r coordinate
     *  @param  innerZCoordinate the inner z coordinate
     *  @param  outerZCoordinate the outer z coordinate
     */
    void MakeXZLayerOutline(float innerRCoordinate, float outerRCoordinate, float innerZCoordinate, float outerZCoordinate);

    /**
     *  @brief compute the polygon corners for the detector outline
     * 
     *  @param symmetryOrder is the number of polygon corners
     *  @param closestDistanceToIp is the distance to the closest points on the polygon
     *  @param phi0 reference angle where to start the polygon
     *  @param coordinates vector of double,double pairs which is filled with the x and y coordinates of the polygon corners
     */  
    void ComputePolygonCorners( int symmetryOrder, double closestDistanceToIp, double phi0, std::vector<std::pair<Double_t,Double_t> > &coordinates );


    /**
     *  @brief Creates a "tube" volume with the given symmetry inside and outside. If a symmetryOrder <= 2 is chosen, a circle is used instead of a polygon
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
    TGeoVolume* MakePolygoneTube( std::string name, int innerSymmetryOrder, int outerSymmetryOrder, 
				  double innerClosestDistanceToIp, double outerClosestDistanceToIp, double innerPhi0, 
				  double outerPhi0, double halfLength, TGeoMedium* medium = 0 );

    /**
     *  @brief Creates a extruded polygonal (or cylindrical) shape with the given symmetry. If a symmetryOrder <= 2 is chosen, a circle is used instead of a polygon
     * 
     *  @param symmetryOrder symmetry order of the polygon (circle if <=2)
     *  @param closestDistanceToIp closest distance between IP and polygon (circle radius)
     *  @param phi starting angle of the polygon
     *  @param halfLength half length (z) of the tube
     */  
    TGeoShape* MakePolygoneTube( int symmetryOrder, double closestDistanceToIp, double phi, double halfLength );


    /**
     *  @brief Computes the corners of a box in 3D
     * 
     *  @param position center of the box
     *  @param normal normal vector of the box (thickness is drawn along this direction)
     *  @param directionU second direction to define the orientation of the box (corresponds to cellSizeU)
     *  @param cellSizeU size of box in direction U
     *  @param cellSizeV size of box in direction V (direction V is computet automatically from normal and directionU)
     *  @param cellSizeThickness thickness of the box
     *  @param corners will be filled with the x,y and z-coordinates of all 8 corners of the box
     */  
    void MakeCaloHitCell( const pandora::CartesianVector& position, 
			  const pandora::CartesianVector& normal, 
			  const pandora::CartesianVector& directionU, 
			  const float cellSizeU, const float cellSizeV, const float cellSizeThickness,
			  Float_t corners[24] );

    /**
     *  @brief  Transform a Pandora monitoring API color enum into a ROOT color enum
     * 
     *  @param  color in Pandora monitoring API enum
     */
    EColor GetROOTColor( Color color );

    /**
     *  @brief  Get a color for a PDG code
     * 
     *  @param  pdgCode of the particle
     */
    Color GetColorForPdgCode(int pdgCode);


    void InitializeEve(Char_t transparency = 70);
    

    static bool                 m_instanceFlag;         ///< The pandora monitoring instance flag
    static PandoraMonitoring    *m_pPandoraMonitoring;  ///< The pandora monitoring instance
    TApplication                *m_pApplication;        ///< The root application

    typedef std::map<const std::string, TH1 *> HistogramMap;

    HistogramMap                m_histogramMap;         ///< The histogram map

    TTreeWrapper                m_treeWrapper;          ///< wrapper around TTree functionality


    typedef std::vector<TObject *> TObjectVector;
    typedef std::vector<TArrow *> TArrowVector;
    typedef std::vector<TPolyMarker *> TPolyMarkerVector;

    typedef std::map<DetectorView,TCanvas*> CanvasMap;

    bool                        m_isOutlineConstructed; ///< Whether the detector outline has been constructed
    TH2F                        *m_pXYAxes;             ///< The xy axes
    TH2F                        *m_pXZAxes;             ///< The xz axes

    TObjectVector               m_2DObjectsXY;          ///< The 2d xy graphs vector
    TObjectVector               m_2DObjectsXZ;          ///< The 2d xz graphs vector
    TArrowVector                m_eventArrows;          ///< The event arrows vector
    TPolyMarkerVector           m_eventMarkers;         ///< The event markers vector

    CanvasMap                   m_canvasMap;            ///< The canvases for each of the detector-views the user has requested
    static bool                 m_eveInitialized;       ///< is set if ROOT Eve is initialized

    static float                m_scalingFactor;        ///< TEve works with [cm], Pandora works with [mm]
    static bool                 m_openEveEvent;         ///< is set if an Event is open to store objects (hits, clusters,...) in it.
    static int                  m_eventDisplayCounter;  ///< counter for the event displays
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline PandoraMonitoring::PandoraMonitoring()
{
    int argc = 0;
    char* argv = (char *)"";

    m_pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
    m_pApplication->SetReturnFromRun(kTRUE);
    m_isOutlineConstructed = false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline PandoraMonitoring::XYOutlineParameters::XYOutlineParameters(int symmetryOrder, float phi0, float closestDistanceToIp) :
    m_symmetryOrder(symmetryOrder),
    m_phi0(phi0),
    m_closestDistanceToIp(closestDistanceToIp)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int PandoraMonitoring::XYOutlineParameters::GetSymmetryOrder() const
{
    return m_symmetryOrder;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PandoraMonitoring::XYOutlineParameters::GetPhi0() const
{
    return m_phi0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float PandoraMonitoring::XYOutlineParameters::GetClosestDistanceToIp() const
{
    return m_closestDistanceToIp;
}



//------------------------------------------------------------------------------------------------------------------------------------------

// Specify (color, ROOT-color)
#define COLOR_TABLE(d)                         \
    d(WHITE,           kWhite)                 \
    d(BLACK,           kBlack)                 \
    d(RED,             kRed)                   \
    d(GREEN,           kGreen)                 \
    d(BLUE,            kBlue)                  \
    d(MAGENTA,         kMagenta)               \
    d(CYAN,            kCyan)                  \
    d(VIOLET,          kViolet)                \
    d(PINK,            kPink)                  \
    d(ORANGE,          kOrange)                \
    d(YELLOW,          kYellow)                \
    d(SPRING,          kSpring)                \
    d(TEAL,            kTeal)                  \
    d(AZURE,           kAzure)                 \
    d(GRAY,            kGray)                  \
    d(DARKRED,         EColor(TColor::GetColorDark(kRed)))       \
    d(DARKGREEN,       EColor(TColor::GetColorDark(kGreen)))     \
    d(DARKBLUE,        EColor(TColor::GetColorDark(kBlue)))      \
    d(DARKMAGENTA,     EColor(TColor::GetColorDark(kMagenta)))   \
    d(DARKCYAN,        EColor(TColor::GetColorDark(kCyan)))      \
    d(DARKVIOLET,      EColor(TColor::GetColorDark(kViolet)))    \
    d(DARKPINK,        EColor(TColor::GetColorDark(kPink)))      \
    d(DARKORANGE,      EColor(TColor::GetColorDark(kOrange)))    \
    d(DARKYELLOW,      EColor(TColor::GetColorDark(kYellow)))    \
    d(LIGHTGREEN,      EColor(TColor::GetColorBright(kGreen)))   \
    d(LIGHTBLUE,       EColor(TColor::GetColorBright(kBlue)))    \
    d(LIGHTMAGENTA,    EColor(TColor::GetColorBright(kMagenta))) \
    d(LIGHTCYAN,       EColor(TColor::GetColorBright(kCyan)))    \
    d(LIGHTVIOLET,     EColor(TColor::GetColorBright(kViolet)))  \
    d(LIGHTPINK,       EColor(TColor::GetColorBright(kPink)))    \
    d(LIGHTORANGE,     EColor(TColor::GetColorBright(kOrange)))  \
    d(LIGHTYELLOW,     EColor(TColor::GetColorBright(kYellow))) 



// Specify (name, pdg code, color)
#define PARTICLE_DATA_COLOR_TABLE(d)                         \
    d(PHOTON,               22,             DARKYELLOW)      \
    d(E_MINUS,              11,             LIGHTBLUE)       \
    d(E_PLUS,              -11,             LIGHTRED)      \
    d(MU_MINUS,             13,             BLUE)      \
    d(MU_PLUS,             -13,             RED)      \
    d(TAU_MINUS,            15,             DARKBLUE)      \
    d(TAU_PLUS,            -15,             DARKRED)      \
    d(NU_E,                 12,             DARKBLUE)      \
    d(NU_E_BAR,            -12,             DARKRED)      \
    d(NU_MU,                14,             DARKBLUE)      \
    d(NU_MU_BAR,           -14,             DARKRED)      \
    d(NU_TAU,               16,             DARKBLUE)      \
    d(NU_TAU_BAR,          -16,             DARKRED)      \
    d(PI_PLUS,             211,             MAGENTA)      \
    d(PI_MINUS,           -211,             VIOLET)      \
    d(PI_ZERO,             111,             LIGHTGREEN)      \
    d(PI_ZERO_BAR,        -111,             LIGHTGREEN)      \
    d(LAMBDA,             3122,             DARKGREEN)      \
    d(LAMBDA_BAR,        -3122,             DARKGREEN)      \
    d(K_PLUS,              321,             DARKGREEN)      \
    d(K_MINUS,            -321,             DARKGREEN)      \
    d(K_SHORT,             310,             LIGHTGREEN)      \
    d(K_SHORT_BAR,        -310,             LIGHTGREEN)      \
    d(K_LONG,              130,             GREEN)      \
    d(K_LONG_BAR,         -130,             GREEN)      \
    d(SIGMA_MINUS,        3112,             GREEN)      \
    d(SIGMA_PLUS,         3222,             GREEN)      \
    d(PROTON,             2212,             ORANGE)      \
    d(NEUTRON,            2112,             CYAN)



/**
 *  @brief  The mass switch statement macro
 */
#define GET_ROOT_COLOR(a, b)                                         \
    case a : return b;

/**
 *  @brief  The mass switch statement macro
 */
#define GET_PARTICLE_COLOR_SWITCH(a, b, c)                                         \
    case a : return c;

/**
 *  @brief  The mass switch statement macro
 */
#define GET_PDG_COLOR_SWITCH(a, b, c)                                         \
    case b : return c;


inline EColor PandoraMonitoring::GetROOTColor(Color color)
{
    switch (color)
    {
        COLOR_TABLE(GET_ROOT_COLOR)
        default : return kGray;
    }
}

inline Color PandoraMonitoring::GetColorForPdgCode(int pdgCode)
{
    switch (pdgCode)
    {
        PARTICLE_DATA_COLOR_TABLE(GET_PDG_COLOR_SWITCH)
        default : return GRAY;
    }
}





} // namespace pandora_monitoring

#endif // #ifndef PANDORA_MONITORING_H



