// ----------------------------------------------------------------------------
// 'JBarrelHCalTreeMaker.h'
// Derek Anderson
// 04.13.2023
//
// A JANA plugin to construct a tree of BHCal tiles
// for training a ML clusterizer.
// ----------------------------------------------------------------------------

// root includes
#include <TTree.h>
#include <TFile.h>
// JANA includes
#include <JANA/JEventProcessor.h>
#include <JANA/JEventProcessorSequentialRoot.h>
// EDM includes
#include <edm4eic/ReconstructedParticle.h>
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/ProtoCluster.h>
#include <edm4eic/Cluster.h>
// ACTS includes
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"
#include <Acts/Surfaces/DiscSurface.hpp>
#include "Acts/Geometry/TrackingGeometry.hpp"
#include <Acts/Material/IMaterialDecorator.hpp>
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
// misc includes
#include <spdlog/spdlog.h>
#include <services/log/Log_service.h>

// global constants
static const size_t NTiles    = 15360;
static const size_t NMaxPars  = 1000;
static const size_t NMaxClust = 20000;



class JBarrelHCalTreeMakerProcessor : public JEventProcessorSequentialRoot {

  private:

    // Data objects we will need from JANA e.g.
    PrefetchT<edm4eic::ReconstructedParticle> genParticles  = {this, "GeneratedParticles"};
    PrefetchT<edm4eic::CalorimeterHit>        bhcalRecHits  = {this, "HcalBarrelRecHits"};
    PrefetchT<edm4eic::Cluster>               bhcalClusters = {this, "HcalBarrelClusters"};
    PrefetchT<edm4eic::Cluster>               becalClusters = {this, "EcalBarrelSciGlassClusters"};

    // i/o members
    TTree *m_tEventTree;
    TTree *m_tClusterTree;

    // tile members (for event tree)
    uint64_t m_numTiles;
    float    m_tileEne[NTiles];
    float    m_tileTime[NTiles];
    float    m_tileTilt[NTiles];
    float    m_tileBarycenter[NTiles];
    uint64_t m_tileTrueID[NTiles];
    short    m_tileIndex[NTiles];
    short    m_tileTower[NTiles];
    short    m_tileSector[NTiles];
    short    m_tileClustIDA[NTiles];
    short    m_tileClustIDB[NTiles];

    // particle members (for cluster tree)
    uint64_t m_numParticles;
    float    m_parEne[NMaxPars];
    float    m_parPhi[NMaxPars];
    float    m_parEta[NMaxPars];

    // bhcal cluster members (for cluster tree)
    uint64_t m_numClustBHCal;
    uint64_t m_bhcalClustNumCells[NMaxClust];
    float    m_bhcalClustEne[NMaxClust];
    float    m_bhcalClustEta[NMaxClust];
    float    m_bhcalClustPhi[NMaxClust];

    // becal cluster members (for cluster tree)
    uint64_t m_numClustBECal;
    uint64_t m_becalClustNumCells[NMaxClust];
    float    m_becalClustEne[NMaxClust];
    float    m_becalClustEta[NMaxClust];
    float    m_becalClustPhi[NMaxClust];

    // utility members
    double                 m_tileTilts[NTiles];
    double                 m_tileBarycenters[NTiles];
    dd4hep::BitFieldCoder* m_decoder;

    // private methods
    void InitializeDecoder();
    void InitializeTrees();
    void InitializeMaps();
    void ResetVariables();

  public:

    // ctor
    JBarrelHCalTreeMakerProcessor() { SetTypeName(NAME_OF_THIS); }

    // inhereited methods
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;

};

// end ------------------------------------------------------------------------
