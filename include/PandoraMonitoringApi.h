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

namespace pandora { class Pandora; }

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  The detector view enum
 */
enum DetectorView
{
    DETECTOR_VIEW_DEFAULT,
    DETECTOR_VIEW_XY,
    DETECTOR_VIEW_XZ
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  The color enum
 */
enum Color
{
    WHITE = 0,
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
     *  @brief  Set a variable in a tree (create the tree and the branch if not yet existing)
     * 
     *  @param  pandora the calling pandora instance
     *  @param  treeName name of the tree (is created if it does not exist yet)
     *  @param  variableName name of the branch in the tree (the branch is created if it does not exist yet)
     *  @param  t sets value of the variable (permitted types are float/double/int and std::vector<float>*,std::vector<double>*,std::vector<int>*
     */
    template <typename T>
    static void SetTreeVariable(const pandora::Pandora &pandora, const std::string &treeName, const std::string &variableName, T t);

    /**
     *  @brief  Fill the tree with the variables which have been set before with SetTreeVariable
     * 
     *  @param  pandora the calling pandora instance
     *  @param  treeName name of the tree to be filled
     */
    static void FillTree(const pandora::Pandora &pandora, const std::string &treeName);

    /**
     *  @brief  Print the tree
     * 
     *  @param  pandora the calling pandora instance
     *  @param  treeName name of the tree to be printed
     */
    static void PrintTree(const pandora::Pandora &pandora, const std::string &treeName);

    /**
     *  @brief  Scan the tree (print a list of all values of all branches)
     * 
     *  @param  pandora the calling pandora instance
     *  @param  treeName name of the tree to be scanned
     */
    static void ScanTree(const pandora::Pandora &pandora, const std::string &treeName);

    /**
     *  @brief  Save the tree to a file
     * 
     *  @param  pandora the calling pandora instance
     *  @param  fileName the file name under which to save the histogram
     *  @param  fileOptions the options associated with opening/recreating a file
     * 
     *  @param  treeName name of the tree to be written to a file
     */
    static void SaveTree(const pandora::Pandora &pandora, const std::string &treeName, const std::string &fileName,
        const std::string &fileOptions);

    /**
     *  @brief  Draw a pandora histogram
     * 
     *  @param  pandora the calling pandora instance
     *  @param  t the pandora histogram
     */
    template <typename T>
    static void DrawPandoraHistogram(const pandora::Pandora &pandora, const T &t);

    /**
     *  @brief  Draw a pandora histogram
     * 
     *  @param  pandora the calling pandora instance
     *  @param  t the pandora histogram
     *  @param  options the drawing options
     */
    template <typename T>
    static void DrawPandoraHistogram(const pandora::Pandora &pandora, const T &t, const std::string &options);

    /**
     *  @brief  Set TEve display parameters
     * 
     *  @param  pandora the calling pandora instance
     *  @param  showDetectors turns the visibility of the detector geometry on or off
     *  @param  detectorView the detector view
     *  @param  transparencyThresholdE cell energy for which transparency is saturated (0%, fully opaque)
     *  @param  energyScaleThresholdE cell energy for which color is at top end of continous color palette
     */
    static void SetEveDisplayParameters(const pandora::Pandora &pandora, const bool showDetectors = true,
        const DetectorView detectorView = DETECTOR_VIEW_DEFAULT, const float transparencyThresholdE = -1.f, const float energyScaleThresholdE = -1.f);

    typedef std::map<int, float> PdgCodeToEnergyMap;

    /**
     *  @brief  Add MCParticles to the Eve event-display
     * 
     *  @param  pandora the calling pandora instance
     *  @param  pMCParticleList list of MC particles to be added to the event display
     *  @param  name of the MC particle list
     *  @param  color The color the track elements are drawn with
     *  @param  pParticleSuppressionMap map from pdg-codes to energy for suppression of particles types below specific energies
     */
    static void VisualizeMCParticles(const pandora::Pandora &pandora, const pandora::MCParticleList *const pMCParticleList,
        const std::string &name, const Color color, const PdgCodeToEnergyMap *pParticleSuppressionMap = NULL);

    /**
     *  @brief  Add Tracks to the Eve event-display
     * 
     *  @param  pandora the calling pandora instance
     *  @param  pTrackList list of tracks to be added to the event display
     *  @param  name of the track list
     *  @param  color The color the track elements are drawn with
     */
    static void VisualizeTracks(const pandora::Pandora &pandora, const pandora::TrackList *const pTrackList,
        const std::string &name, const Color color);

    /**
     *  @brief Add CaloHits to the Eve event-display
     * 
     *  @param  pandora the calling pandora instance
     *  @param pCaloHitList list of calohits to be added to the event display
     *  @param name of the calohit list
     *  @param color The color the cluster elements are drawn with
     */  
    static void VisualizeCaloHits(const pandora::Pandora &pandora, const pandora::CaloHitList *const pCaloHitList,
        const std::string &name, const Color color);

    /**
     *  @brief  Add Clusters to the Eve event-display
     * 
     *  @param  pandora the calling pandora instance
     *  @param  pClusterList list of clusters to be added to the event display
     *  @param  name of the cluster list
     *  @param  color The color the cluster elements are drawn with
     *  @param  showAssociatedTracks draw the tracks associated to the cluster
     */
    static void VisualizeClusters(const pandora::Pandora &pandora, const pandora::ClusterList *const pClusterList,
        const std::string &name, const Color color, bool showAssociatedTracks = false);

    /**
     *  @brief  Add Particle flow objects to the Eve event-display
     * 
     *  @param  pandora the calling pandora instance
     *  @param  pPfoList list of particle flow objects to be added to the event display
     *  @param  name of the pfo list
     *  @param  parent pointer to the parent TEveElement. If NULL, the cluster will be parent element
     *  @param  color The color the cluster elements are drawn with
     *  @param  showAssociatedTracks draw the tracks associated to the cluster
     *  @param  showVertices whether to draw markers to represent the pfo vertices
     *  @param  displayPfoHierarchy whether to draw daughter particles (only) within parent elements
     */
    static void VisualizeParticleFlowObjects(const pandora::Pandora &pandora, const pandora::PfoList *const pPfoList,
        const std::string &name, const Color color, bool showVertices = true, bool displayPfoHierarchy = true);

    /**
     *  @brief  Add Vertices to the Eve event-display
     * 
     *  @param  pandora the calling pandora instance
     *  @param  pVertexList list of vertices to be added to the event display
     *  @param  name of the vertex list
     *  @param  color The color the cluster elements are drawn with
     */
    static void VisualizeVertices(const pandora::Pandora &pandora, const pandora::VertexList *const pVertexList,
        const std::string &name, const Color color);

    /**
     *  @brief Add marker to visualization
     * 
     *  @param  pandora the calling pandora instance
     *  @param pMarkerPoint address of the marker point
     *  @param name of the marker
     *  @param color the marker color
     *  @param showAssociatedTracks draw the tracks associated to the cluster
     */
    static void AddMarkerToVisualization(const pandora::Pandora &pandora, const pandora::CartesianVector *const pMarkerPoint,
        const std::string &name, const Color color, const unsigned int markerSize);

    /**
     *  @brief  Show the Eve Event-display and pause.
     * 
     *  @param  pandora the calling pandora instance
     */
    static void ViewEvent(const pandora::Pandora &pandora);

    /**
     *  @brief  Pause thread until user enters 'return'
     * 
     *  @param  pandora the calling pandora instance
     */
    static void Pause(const pandora::Pandora &pandora);

    /**
     *  @brief  Create monitoring instance, with its associated ROOT TApplication
     * 
     *  @param  pandora the calling pandora instance
     */
    static void Create(const pandora::Pandora &pandora);

    /**
     *  @brief  Delete monitoring instance, terminating its associated ROOT TApplication (unless it is in use by another instance)
     * 
     *  @param  pandora the calling pandora instance
     */
    static void Delete(const pandora::Pandora &pandora);
};

#endif // #ifndef PANDORA_MONITORING_API_H
