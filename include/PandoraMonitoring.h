/**
 *  @file   PandoraMonitoring/include/PandoraMonitoring.h
 *
 *  @brief  Header file for the pandora monitoring class.
 *
 *  $Log: $
 */
#ifndef PANDORA_MONITORING_H
#define PANDORA_MONITORING_H 1

#include "PandoraInternal.h"

#include "TApplication.h"

#include <iostream>

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
     *  @brief  Destructor
     */
    ~PandoraMonitoring();

    /**
     *  @brief  Draw a test canvas
     */
    void DrawCanvas();

    /**
     *  @brief  Loop through clusters, displaying parent addresses of constituent calo hits
     */
    void LookAtClusters(const pandora::ClusterList *const pClusterList);

    /**
     *  @brief  Pause until user enters 'return'
     */
    void Pause();

private:
    /**
     *  @brief  Default constructor
     */
    PandoraMonitoring();

    static bool                 m_instanceFlag;         ///< The pandora monitoring instance flag
    static PandoraMonitoring    *m_pPandoraMonitoring;  ///< The pandora monitoring instance
    TApplication                *m_pApplication;        ///< The root application
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline PandoraMonitoring::~PandoraMonitoring()
{
    m_instanceFlag = false;
    delete m_pPandoraMonitoring;
    m_pPandoraMonitoring = NULL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline PandoraMonitoring::PandoraMonitoring()
{
    int argc = 0;
    char* argv = "";

    m_pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
    m_pApplication->SetReturnFromRun(kTRUE);
}

} // namespace pandora_monitoring

#endif // #ifndef PANDORA_MONITORING_H
