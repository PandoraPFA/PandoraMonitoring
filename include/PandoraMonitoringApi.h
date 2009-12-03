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
    template <typename VariableType>
    static void SetTreeVariable(const std::string &treeName, const std::string &variableName, VariableType  variable);

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
     *  @brief  Save the tree to a file
     *  @param  fileName the file name under which to save the histogram
     *  @param  fileOptions the options associated with opening/recreating a file
     * 
     *  @param  treeName name of the tree to be written to a file
     */
    static void SaveTree(const std::string &treeName, const std::string &fileName, const std::string &fileOptions );

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
     *  @brief  Draw tracks in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     */ 
    static void DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList);

    /**
     *  @brief  Draw calo hits in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pOrderedCaloHitList address of the ordered calo hit list
     */ 
    static void DrawEvent(DetectorView detectorView, const pandora::OrderedCaloHitList *const pOrderedCaloHitList);

    /**
     *  @brief  Draw clusters in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pClusterList address of the cluster list
     */ 
    static void DrawEvent(DetectorView detectorView, const pandora::ClusterList *const pClusterList);

    /**
     *  @brief  Draw tracks and calo hits in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     *  @param  pOrderedCaloHitList address of the ordered calo hit list
     */ 
    static void DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList, const pandora::OrderedCaloHitList *const pOrderedCaloHitList);

    /**
     *  @brief  Draw tracks and clusters in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     *  @param  pClusterList address of the cluster list
     */ 
    static void DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList, const pandora::ClusterList *const pClusterList);

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
};

#endif // #ifndef PANDORA_MONITORING_API_H
