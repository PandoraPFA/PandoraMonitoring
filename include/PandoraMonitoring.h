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

#include <iostream>
#include <map>

class TH1;
class TH2F;
class TObject;
class TPolyMarker;

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
     *  @brief  Draw tracks in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     */ 
    void DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList);

    /**
     *  @brief  Draw clusters in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pClusterList address of the cluster list
     */ 
    void DrawEvent(DetectorView detectorView, const pandora::ClusterList *const pClusterList);

    /**
     *  @brief  Draw tracks and clusters in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     *  @param  pClusterList address of the cluster list
     */ 
    void DrawEvent(DetectorView detectorView, const pandora::TrackList *const pTrackList, const pandora::ClusterList *const pClusterList);

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
        XYOutlineParameters(int symmetryOrder, float phi0, float closestDistanceToIp) :
            m_symmetryOrder(symmetryOrder),
            m_phi0(phi0),
            m_closestDistanceToIp(closestDistanceToIp)
        {
        }

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
     *  @brief  Draw tracks in an event
     * 
     *  @param  detectorView the detector view
     *  @param  pTrackList address of the track list
     */ 
    void DrawTracks(DetectorView detectorView, const pandora::TrackList *const pTrackList);

    /**
     *  @brief  Draw clusters in an event
     *
     *  @param  detectorView the detector view
     *  @param  pClusterList address of the cluster list
     */ 
    void DrawClusters(DetectorView detectorView, const pandora::ClusterList *const pClusterList);

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

    static bool                 m_instanceFlag;         ///< The pandora monitoring instance flag
    static PandoraMonitoring    *m_pPandoraMonitoring;  ///< The pandora monitoring instance
    TApplication                *m_pApplication;        ///< The root application

    typedef std::map<const std::string, TH1 *> HistogramMap;

    HistogramMap                m_histogramMap;         ///< The histogram map

    typedef std::vector<TObject *> TObjectVector;
    typedef std::vector<TPolyMarker *> TPolyMarkerVector;

    bool                        m_isOutlineConstructed; ///< Whether the detector outline has been constructed
    TH2F                        *m_pXYAxes;             ///< The xy axes
    TH2F                        *m_pXZAxes;             ///< The xz axes

    TObjectVector               m_2DObjectsXY;          ///< The 2d xy graphs vector
    TObjectVector               m_2DObjectsXZ;          ///< The 2d xz graphs vector
    TPolyMarkerVector           m_eventMarkers;         ///< The event markers vector
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

} // namespace pandora_monitoring

#endif // #ifndef PANDORA_MONITORING_H
