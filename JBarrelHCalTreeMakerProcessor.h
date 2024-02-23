// ----------------------------------------------------------------------------
// 'JBarrelHCalTreeMaker.h'
// Derek Anderson
// 04.13.2023
//
// A JANA plugin to construct a tree of BHCal tiles
// for training a ML clusterizer.
// ----------------------------------------------------------------------------

// c includes
#include <cmath>
#include <vector>
#include <cstdlib>
// root includes
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
// JANA includes
#include <JANA/JEventProcessor.h>
#include <JANA/JEventProcessorSequentialRoot.h>
// EDM4EIC includes
#include <edm4eic/Cluster.h>
#include <edm4eic/TrackerHit.h>
#include <edm4eic/ProtoCluster.h>
#include <edm4eic/RawTrackerHit.h>
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/ReconstructedParticle.h>
// misc includes
#include <spdlog/spdlog.h>
#include <services/log/Log_service.h>
#include <services/rootfile/RootFile_service.h>
#include <services/geometry/dd4hep/DD4hep_service.h>
// DD4HEP includes
#include "DD4hep/Objects.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/IDDescriptor.h"
// DDRec includes
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/CellIDPositionConverter.h"
// DDG4 includes
#include "DDG4/Geant4Data.h"

// global constants
static const size_t NTiles    = 15360;
static const size_t NMaxPars  = 1000;
static const size_t NMaxClust = 20000;

// calculation parameters
static const bool AddBECalClusters = false;



class JBarrelHCalTreeMakerProcessor : public JEventProcessorSequentialRoot {

  private:

    // data objects we will definitely need from JANA
    PrefetchT<edm4eic::ReconstructedParticle> genParticles  = {this, "GeneratedParticles"};
    PrefetchT<edm4eic::CalorimeterHit>        bhcalRecHits  = {this, "HcalBarrelRecHits"};
    PrefetchT<edm4eic::Cluster>               bhcalClusters = {this, "HcalBarrelClusters"};

    // name of BECal cluster collection (only used if needed)
    std::string m_becalClustName = "EcalBarrelClusters";

    // i/o members
    TTree      *m_tEventTree;
    TTree      *m_tClusterTree;
    TDirectory *m_dPluginDir;

    // vector of hit tiles (for matching to clusters)
    std::vector<uint64_t> m_vecTileID;
    std::vector<uint64_t> m_vecTileIsMatched;

    // tile members (for event tree)
    uint64_t              m_numTiles;
    std::vector<float>    m_tileEne;
    std::vector<float>    m_tileTime;
    std::vector<float>    m_tileTilt;
    std::vector<float>    m_tileBarycenter;
    std::vector<uint64_t> m_tileCellID;
    std::vector<uint64_t> m_tileTrueID;
    std::vector<short>    m_tileIndex;
    std::vector<short>    m_tileTower;
    std::vector<short>    m_tileSector;
    std::vector<short>    m_tileClustIDA;
    std::vector<short>    m_tileClustIDB;
    /* TODO add
     *   - x, y, z 
     *   - hit branches for clusters
     */

    // particle members (for cluster tree)
    uint64_t           m_numParticles;
    std::vector<float> m_parEne;
    std::vector<float> m_parPhi;
    std::vector<float> m_parEta;
    /* TODO add
     *   - PDG
     *   - mass
     *   - MCParticles w/
     *     + PDG
     *     + generator, simulator
     *     + energy, eta, phi, mass
     *   - particle branches for clusters
     *     + those sitting in associations
     *     + those connected via contributions
     */

    // bhcal cluster members (for cluster tree)
    uint64_t              m_numClustBHCal;
    std::vector<uint64_t> m_bhcalClustNumCells;
    std::vector<float>    m_bhcalClustEne;
    std::vector<float>    m_bhcalClustEta;
    std::vector<float>    m_bhcalClustPhi;
    /* TODO add
     *   - above branches
     *   - time
     */

    // becal cluster members (for cluster tree)
    uint64_t              m_numClustBECal;
    std::vector<uint64_t> m_becalClustNumCells;
    std::vector<float>    m_becalClustEne;
    std::vector<float>    m_becalClustEta;
    std::vector<float>    m_becalClustPhi;

    // utility members
    std::vector<double>    m_tileTilts;
    std::vector<double>    m_tileBarycenters;
    dd4hep::BitFieldCoder* m_decoder;

    // private methods
    void InitializeDecoder();
    void InitializeTrees();
    void InitializeMaps();
    void ResetVariables();
    void ResizeEvtTiles(const size_t nTiles);
    void ResizeEvtPars(const size_t nPars);
    void ResizeBHCalClusts(const size_t nBHCalClust);
    void ResizeBECalClusts(const size_t nBECalClust);

  public:

    // ctor
    JBarrelHCalTreeMakerProcessor() { SetTypeName(NAME_OF_THIS); }

    // inhereited methods
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;

};

// end ------------------------------------------------------------------------
