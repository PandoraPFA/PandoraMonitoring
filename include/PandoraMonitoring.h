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

#include "TApplication.h"

#include <iostream>
#include <map>

class TH1;

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
     *  @brief  Save a histogram to a file
     * 
     *  @param  name the name associated with the histogram
     *  @param  fileName the file name under which to save the histogram
     */
    void SaveHistogram(const std::string &name, const std::string &fileName);

    /**
     *  @brief  Delete a histogram
     * 
     *  @param  name the name associated with the histogram
     */
    void DeleteHistogram(const std::string &name);

    /**
     *  @brief  Temporary function - draw the detector outline
     */ 
    void DrawDetectorOutline() const;

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
     *  @brief  Default constructor
     */
    PandoraMonitoring();

    static bool                 m_instanceFlag;         ///< The pandora monitoring instance flag
    static PandoraMonitoring    *m_pPandoraMonitoring;  ///< The pandora monitoring instance
    TApplication                *m_pApplication;        ///< The root application

    typedef std::map<const std::string, TH1 *> HistogramMap;

    HistogramMap                m_histogramMap;         ///< The histogram map
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline PandoraMonitoring::PandoraMonitoring()
{
    int argc = 0;
    char* argv = (char *)"";

    m_pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
    m_pApplication->SetReturnFromRun(kTRUE);
}

} // namespace pandora_monitoring

#endif // #ifndef PANDORA_MONITORING_H
