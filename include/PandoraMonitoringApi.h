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
 *  @brief  The color enum
 */
enum Color
{
    WHITE=0,
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
    AUTO,       // automatic choice of colors
    AUTOID,     // automatic choice of colors depending on the particle ID
    AUTOTYPE,   // automatic choice of colors depending on the particle type 
    AUTOITER,   // automatic choice of colors iterating through colors
    AUTOENERGY  // continuous color palette indicating hit energies
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
    static void Create1DHistogram(const std::string &name, const std::string &title, int nBinsX, float xLow, float xUp,
        const std::string xAxisTitle = "", const std::string yAxisTitle = "");

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
        double yLow, double yUp, const std::string xAxisTitle = "", const std::string yAxisTitle = "");

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
     *  @brief  Add
     * 
     *  @param  nameHisto0 the name associated with the first histogram
     *  @param  nameHisto1 the name associated with the second histogram
     *  @param  coeff0 coefficient for the first histogram
     *  @param  coeff1 coefficient for the second histogram
     */
    static void AddHistograms(const std::string &nameHisto0, const std::string &nameHisto1, double coeff0, double coeff1);

    /**
     *  @brief  Multiply
     * 
     *  @param  nameHisto0 the name associated with the first histogram
     *  @param  nameHisto1 the name associated with the second histogram
     *  @param  coeff0 coefficient for the first histogram
     *  @param  coeff1 coefficient for the second histogram
     */
    static void MultiplyHistograms(const std::string &nameHisto0, const std::string &nameHisto1, double coeff0, double coeff1);

    /**
     *  @brief  Divide two histograms
     * 
     *  @param  nameHisto0 the name associated with the first histogram
     *  @param  nameHisto1 the name associated with the second histogram
     *  @param  coeff0 coefficient for the first histogram
     *  @param  coeff1 coefficient for the second histogram
     */
    static void DivideHistograms(const std::string &nameHisto0, const std::string &nameHisto1, double coeff0, double coeff1);

    /**
     *  @brief  Set a variable in a tree (create the tree and the branch if not yet existing)
     * 
     *  @param  treeName name of the tree (is created if it does not exist yet)
     *  @param  variableName name of the branch in the tree (the branch is created if it does not exist yet)
     *  @param  t sets value of the variable (permitted types are float/double/int and std::vector<float>*,std::vector<double>*,std::vector<int>*
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
     *  @brief  Draw a pandora histogram
     * 
     *  @param  t the pandora histogram
     */
    template <typename T>
    static void DrawPandoraHistogram(const T &t);

    /**
     *  @brief  Draw a pandora histogram
     * 
     *  @param  t the pandora histogram
     *  @param  options the drawing options
     */
    template <typename T>
    static void DrawPandoraHistogram(const T &t, const std::string &options);

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
     *  @brief  Show the Eve Event-display and pause.
     */
    static void ViewEvent();

    typedef std::map<int, float> PdgCodeToEnergyMap;

    /**
     *  @brief  Set TEve display parameters
     * 
     *  @param  blackBackground whether to use a black background color, rather than white
     *  @param  showDetectors turns the visibility of the detector geometry on or off
     *  @param  transparencyThresholdE cell energy for which transparency is saturated (0%, fully opaque)
     *  @param  energyScaleThresholdE cell energy for which color is at top end of continous color palette
     */
    static void SetEveDisplayParameters(const bool blackBackground, const bool showDetectors, const float transparencyThresholdE = -1.f,
        const float energyScaleThresholdE = -1.f);

    /**
     *  @brief Add MCParticles to the Eve event-display
     * 
     *  @param pMCParticleList list of MC particles to be added to the event display
     *  @param name of the MC particle list
     *  @param color The color the track elements are drawn with
     *  @param pParticleSuppressionMap map from pdg-codes to energy for suppression of particles types below specific energies
     */  
    static void VisualizeMCParticles(const pandora::MCParticleList *const pMCParticleList, std::string name, Color color,
        const PdgCodeToEnergyMap *pParticleSuppressionMap = NULL);

    /**
     *  @brief Add Tracks to the Eve event-display
     * 
     *  @param pTrackList list of tracks to be added to the event display
     *  @param name of the track list
     *  @param color The color the track elements are drawn with
     */  
    static void VisualizeTracks(const pandora::TrackList *const pTrackList, std::string name, Color color);

    /**
     *  @brief Add CaloHits to the Eve event-display
     * 
     *  @param pCaloHitList list of calohits to be added to the event display
     *  @param name of the calohit list
     *  @param color The color the cluster elements are drawn with
     */  
    static void VisualizeCaloHits(const pandora::CaloHitList *const pCaloHitList, std::string name, Color color);

    /**
     *  @brief Add Particle flow objects to the Eve event-display
     * 
     *  @param pPfoList list of particle flow objects to be added to the event display
     *  @param name of the pfo list
     *  @param parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param color The color the cluster elements are drawn with
     *  @param showAssociatedTracks draw the tracks associated to the cluster
     */
    static void VisualizeParticleFlowObjects(const pandora::PfoList *const pPfoList, std::string name, Color color,
        bool showAssociatedTracks = true);

    /**
     *  @brief Add Clusters to the Eve event-display
     * 
     *  @param pClusterList list of clusters to be added to the event display
     *  @param name of the cluster list
     *  @param color The color the cluster elements are drawn with
     *  @param showAssociatedTracks draw the tracks associated to the cluster
     */  
    static void VisualizeClusters(const pandora::ClusterList *const pClusterList, std::string name, Color color,
        bool showAssociatedTracks = true);

    /**
     *  @brief  Delete monitoring instance
     */
    static void Delete();
};

#endif // #ifndef PANDORA_MONITORING_API_H
