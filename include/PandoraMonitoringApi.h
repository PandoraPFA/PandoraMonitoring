/**
 *  @file   PandoraMonitoring/include/PandoraMonitoringApi.h
 *
 *  @brief  Header file for the pandora monitoring api class.
 *
 *  $Log: $
 */
#ifndef PANDORA_MONITORING_API_H
#define PANDORA_MONITORING_API_H 1

#include "Pandora/PandoraInternal.h"

/**
 *  @brief  The detector view enum
 */
enum DetectorView
{
    DETECTOR_VIEW_XY,
    DETECTOR_VIEW_XZ
};

/**
 *  @brief  The color enum
 */
enum Color
{
    WHITE,
    BLACK,
    RED,
    GREEN,
    BLUE,
    MAGENTA,
    CYAN,
    VIOLET,
    PINK,
    ORANGE,
    YELLOW,
    SPRING,
    TEAL,
    AZURE,
    GRAY,
    DARKRED,
    DARKGREEN,
    DARKBLUE,
    DARKMAGENTA,
    DARKCYAN,
    DARKVIOLET,
    DARKPINK,
    DARKORANGE,
    DARKYELLOW,
    LIGHTRED,
    LIGHTGREEN,
    LIGHTBLUE,
    LIGHTMAGENTA,
    LIGHTCYAN,
    LIGHTVIOLET,
    LIGHTPINK,
    LIGHTORANGE,
    LIGHTYELLOW,
    AUTO // automatic choice of colors
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  PandoraMonitoringApi class
 */
class PandoraMonitoringApi
{
public:
    /**
     *  @brief  Create a 1D histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  title the histogram title
     *  @param  nBinsX the number of x bins
     *  @param  xLow the the lower bound for the x axis
     *  @param  xUp the upper bound for the x axis
     */
    static void Create1DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp);

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
    static void Create2DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp, int nBinsY,
        double yLow, double yUp);

    /**
     *  @brief  Add a single entry to a 1D histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  xValue the x value to add
     *  @param  weight the weight to apply to the entry
     */
    static void Fill1DHistogram(const std::string &name, float xValue, float weight = 1);

    /**
     *  @brief  Add a single entry to a 2D histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  xValue the x value to add
     *  @param  yValue the y value to add
     *  @param  weight the weight to apply to the entry
     */
    static void Fill2DHistogram(const std::string &name, float xValue, float yValue, float weight = 1);

    /**
     *  @brief  Set a variable in a tree (create the tree and the branch if not yet existing)
     * 
     *  @param  treeName name of the tree (is created if it does not exist yet)
     *  @param  variableName name of the branch in the tree (the branch is created if it does not exist yet)
     *  @param  variable sets value of the variable (permitted types are float/double/int and std::vector<float>*,std::vector<double>*,std::vector<int>*
     */
    template <typename T>
    static void SetTreeVariable(const std::string &treeName, const std::string &variableName, T t);

    /**
     *  @brief  Fill the tree with the variables which have been set before with SetTreeVariable
     * 
     *  @param  treeName name of the tree to be filled
     */
    static void FillTree(const std::string &treeName);

    /**
     *  @brief  Print the tree
     * 
     *  @param  treeName name of the tree to be printed
     */
    static void PrintTree(const std::string &treeName);

    /**
     *  @brief  Scan the tree (print a list of all values of all branches)
     * 
     *  @param  treeName name of the tree to be scanned
     */
    static void ScanTree(const std::string &treeName);

    /**
     *  @brief  Save the tree to a file
     *  @param  fileName the file name under which to save the histogram
     *  @param  fileOptions the options associated with opening/recreating a file
     * 
     *  @param  treeName name of the tree to be written to a file
     */
    static void SaveTree(const std::string &treeName, const std::string &fileName, const std::string &fileOptions);

    /**
     *  @brief  Draw a histogram
     * 
     *  @param  name the name associated with the histogram
     */
    static void DrawHistogram(const std::string &name);

    /**
     *  @brief  Draw a histogram
     * 
     *  @param  name the name associated with the histogram
     *  @param  options the drawing options
     */
    static void DrawHistogram(const std::string &name, const std::string &options);

    /**
     *  @brief  Save a histogram to a file and tidy up
     * 
     *  @param  name the name associated with the histogram
     *  @param  fileName the file name under which to save the histogram
     *  @param  fileOptions the options associated with opening/recreating a file
     */
    static void SaveAndCloseHistogram(const std::string &name, const std::string &fileName, const std::string &fileOptions);

    /**
     *  @brief  Delete a histogram
     * 
     *  @param  name the name associated with the histogram
     */
    static void DeleteHistogram(const std::string &name);

    /**
     *  @brief Show the Eve Event-display and pause.
     * 
     */  
    static void View();


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
    static void VisualizeParticleFlowObjects(const pandora::ParticleFlowObjectList *const pPfoList, std::string name, Color color, 
					     bool showAssociatedTracks = true, bool showFit = true );

    /**
     *  @brief Add Clusters to the Eve event-display
     * 
     *  @param pClusterList list of clusters to be added to the event display
     *  @param name of the cluster list
     *  @param color The color the cluster elements are drawn with
     *  @param showAssociatedTracks draw the tracks associated to the cluster
     *  @param showFit draw an arrow representing the fit through the calorimeterhits (the fit is computed within pandora)
     */  
    static void VisualizeClusters(const pandora::ClusterList *const pClusterList, std::string name, Color color,
				  bool showAssociatedTracks = true, bool showFit = true );

    /**
     *  @brief Add Tracks to the Eve event-display
     * 
     *  @param pTrackList list of tracks to be added to the event display
     *  @param name of the track list
     *  @param color The color the track elements are drawn with
     */  
    static void VisualizeTracks(const pandora::TrackList *const pTrackList, std::string name, Color color );

    /**
     *  @brief Add CaloHits to the Eve event-display
     * 
     *  @param pOrderedCaloHitList list of calohits to be added to the event display
     *  @param name of the calohit list
     *  @param color The color the cluster elements are drawn with
     */  
    static void VisualizeCaloHits(const pandora::OrderedCaloHitList *const pOrderedCaloHitList, std::string name, Color color );

    /**
     *  @brief  Pauses the monitoring, such that the user can see the output. Clear the canvases and other data. Waits for key-press
     * 
     */ 
    static void ViewEvent();

    /**
     *  @brief  Add ClusterList to the output. The canvas is automatically created if not existing. Pause the output with ViewEvent()
     * 
     *  @param  detectorView the detector view
     *  @param  pClusterList address of the cluster list
     *  @param  color in which the clusters should be drawn
     */ 
    static void AddClusterList(DetectorView detectorView, const pandora::ClusterList *const pClusterList, Color color = AUTO);

    /**
     *  @brief  Add TrackList to the output. The canvas is automatically created if not existing. Pause the output with ViewEvent()
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     *  @param  color in which the tracks should be drawn
     */ 
    static void AddTrackList(DetectorView detectorView, const pandora::TrackList *const pTrackList, Color color = AUTO);

    /**
     *  @brief  Add OrderedCaloHitList to the output. The canvas is automatically created if not existing. Pause the output with ViewEvent()
     * 
     *  @param  detectorView the detector view
     *  @param  pOrderedCaloHitList address of the orderedCaloHit list
     *  @param  color in which the orderedCaloHits should be drawn
     */ 
    static void AddCaloHitList(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList, Color color = AUTO);

    /**
     *  @brief  Draw the detector outline
     * 
     *  @param  detectorView the detector view
     */ 
    static void DrawDetectorOutline(DetectorView detectorView);

    /**
     *  @brief  Temporary function - draw a test canvas and histogram
     */
    static void Test();

    /**
     *  @brief  Temporary function - display the parent addresses of calo hits contained in a cluster list
     * 
     *  @param  pClusterList address of the cluster list
     */
    static void Test2(const pandora::ClusterList *const pClusterList);

    /**
     *  @brief  Delete monitoring instance
     * 
     */
    static void Delete();
};

#endif // #ifndef PANDORA_MONITORING_API_H
