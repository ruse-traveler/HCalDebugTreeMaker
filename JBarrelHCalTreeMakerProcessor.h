// ----------------------------------------------------------------------------
// 'JBarrelHCalTreeMaker.h'
// Derek Anderson
// 04.13.2023
//
// A JANA plugin to construct a tree of BHCal tiles
// for training a ML clusterizer.
// ----------------------------------------------------------------------------

// c++ utilities
#include <map>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <utility>
// root libraries
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>
// JANA libraries
#include <JANA/JEventProcessor.h>
#include <JANA/JEventProcessorSequentialRoot.h>
// data model types
#include <edm4eic/Cluster.h>
#include <edm4eic/TrackerHit.h>
#include <edm4eic/ProtoCluster.h>
#include <edm4eic/RawTrackerHit.h>
#include <edm4eic/CalorimeterHit.h>
#include <edm4eic/ReconstructedParticle.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimCalorimeterHit.h>
// eicrecon services
#include <spdlog/spdlog.h>
#include <services/log/Log_service.h>
#include <services/rootfile/RootFile_service.h>
#include <services/geometry/dd4hep/DD4hep_service.h>
// dd4hep utilities
#include "DD4hep/Objects.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/IDDescriptor.h"
#include "DDG4/Geant4Data.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/CellIDPositionConverter.h"

// calculation parameters
static const bool AddBECalClusters = false;



class JBarrelHCalTreeMakerProcessor : public JEventProcessorSequentialRoot {

  public:

    // ctor
    JBarrelHCalTreeMakerProcessor() { SetTypeName(NAME_OF_THIS); }

    // inhereited methods
    void InitWithGlobalRootLock() override;
    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    void FinishWithGlobalRootLock() override;

  private:

    // data objects we will definitely need from JANA
    PrefetchT<edm4hep::MCParticle>            mcParticles   = {this, "MCParticles"};
    PrefetchT<edm4eic::ReconstructedParticle> genParticles  = {this, "GeneratedParticles"};
    PrefetchT<edm4hep::SimCalorimeterHit>     bhcalSimHits  = {this, "HcalBarrelHits"};
    PrefetchT<edm4eic::CalorimeterHit>        bhcalRecHits  = {this, "HcalBarrelRecHits"};
    PrefetchT<edm4eic::Cluster>               bhcalClusters = {this, "HcalBarrelClusters"};

    // name of BECal cluster collection (only used if needed)
    std::string m_becalClustName = "EcalBarrelClusters";

    // i/o members
    TTree      *m_tEventTree;
    TTree      *m_tClusterTree;
    TDirectory *m_dPluginDir;

    // maps for matching to clusters & particles
    std::map<int64_t, bool>    m_mapTileIDToIsMatched;
    std::map<int64_t, int64_t> m_mapTileIndexToID;
    std::map<int64_t, int64_t> m_mapTileIDToClustIndex;
    std::map<int64_t, int64_t> m_mapContribIDToTileID;
    std::map<int64_t, int64_t> m_mapContribIDToClustIndex;

    // DEPRECATED
    std::vector<uint64_t> m_vecTileID;
    std::vector<uint64_t> m_vecTileIsMatched;

    // tile members (for event tree)
    uint64_t              m_evtIndex;
    uint64_t              m_numTiles;
    std::vector<float>    m_tileEne;
    std::vector<float>    m_tilePosX;
    std::vector<float>    m_tilePosY;
    std::vector<float>    m_tilePosZ;
    std::vector<float>    m_tileTime;
    std::vector<uint64_t> m_tileCellID;
    std::vector<short>    m_tileIndex;
    std::vector<short>    m_tileTower;
    std::vector<short>    m_tileSector;
    std::vector<short>    m_tileClustIDA;
    std::vector<short>    m_tileClustIDB;
    std::vector<uint64_t> m_tileTrueID;

    // generated particle members (for event tree)
    uint64_t             m_numGenParticles;
    std::vector<int32_t> m_genParType;
    std::vector<int32_t> m_genParPDG;
    std::vector<float>   m_genParEne;
    std::vector<float>   m_genParPhi;
    std::vector<float>   m_genParEta;
    std::vector<float>   m_genParMass;

    // mc particle members (for event tree)
    uint64_t             m_numMCParticles;
    std::vector<int32_t> m_mcParGenStat;
    std::vector<int32_t> m_mcParSimStat;
    std::vector<int32_t> m_mcParPDG;
    std::vector<float>   m_mcParEne;
    std::vector<float>   m_mcParPhi;
    std::vector<float>   m_mcParEta;
    std::vector<float>   m_mcParMass;
    std::vector<float>   m_mcParStartVx;
    std::vector<float>   m_mcParStartVy;
    std::vector<float>   m_mcParStartVz;
    std::vector<float>   m_mcParStopVx;
    std::vector<float>   m_mcParStopVy;
    std::vector<float>   m_mcParStopVz;
    std::vector<float>   m_mcParTime;

    // bhcal cluster members (for cluster tree)
    uint64_t              m_numClustBHCal;
    std::vector<uint64_t> m_bhcalClustIndex;
    std::vector<uint64_t> m_bhcalClustNumCells;
    std::vector<float>    m_bhcalClustEne;
    std::vector<float>    m_bhcalClustEta;
    std::vector<float>    m_bhcalClustPhi;
    std::vector<float>    m_bhcalClustPosX;
    std::vector<float>    m_bhcalClustPosY;
    std::vector<float>    m_bhcalClustPosZ;
    std::vector<float>    m_bhcalClustTime;
    std::vector<uint64_t> m_bhcalClustNumTileInClust;
    std::vector<uint64_t> m_bhcalClustNumMCParAssocToClust;
    std::vector<uint64_t> m_bhcalClustNumMCParContribToClust;

    // bhcal tiles in cluster members (for cluster tree)
    std::vector<uint64_t> m_tileInClustID;
    std::vector<uint64_t> m_tileInClustClustIndex;
    std::vector<float>    m_tileInClustEne;
    std::vector<float>    m_tileInClustPosX;
    std::vector<float>    m_tileInClustPosY;
    std::vector<float>    m_tileInClustPosZ;
    std::vector<float>    m_tileInClustTime;
    std::vector<uint64_t> m_tileInClustCellID;
    std::vector<uint64_t> m_tileInClustTrueID;
    std::vector<short>    m_tileInClustIndex;
    std::vector<short>    m_tileInClustTower;
    std::vector<short>    m_tileInClustSector;
    std::vector<short>    m_tileInClustClustIDA;
    std::vector<short>    m_tileInClustClustIDB;

    // mc particle associated to cluster members (for cluster tree)
    std::vector<uint64_t> m_mcParAssocToClustClustIndex;
    std::vector<int32_t>  m_mcParAssocToClustGenStat;
    std::vector<int32_t>  m_mcParAssocToClustSimStat;
    std::vector<int32_t>  m_mcParAssocToClustPDG;
    std::vector<float>    m_mcParAssocToClustEne;
    std::vector<float>    m_mcParAssocToClustPhi;
    std::vector<float>    m_mcParAssocToClustEta;
    std::vector<float>    m_mcParAssocToClustMass;
    std::vector<float>    m_mcParAssocToClustStartVx;
    std::vector<float>    m_mcParAssocToClustStartVy;
    std::vector<float>    m_mcParAssocToClustStartVz;
    std::vector<float>    m_mcParAssocToClustStopVx;
    std::vector<float>    m_mcParAssocToClustStopVy;
    std::vector<float>    m_mcParAssocToClustStopVz;
    std::vector<float>    m_mcParAssocToClustTime;

    // mc particle contributing to cluster members (for cluster tree)
    std::vector<uint64_t> m_mcParContribToClustClustIndex;
    std::vector<uint64_t> m_mcParContribToClustCellID;
    std::vector<int32_t>  m_mcParContribToClustGenStat;
    std::vector<int32_t>  m_mcParContribToClustSimStat;
    std::vector<int32_t>  m_mcParContribToClustPDG;
    std::vector<float>    m_mcParContribToClustEne;
    std::vector<float>    m_mcParContribToClustPhi;
    std::vector<float>    m_mcParContribToClustEta;
    std::vector<float>    m_mcParContribToClustMass;
    std::vector<float>    m_mcParContribToClustStartVx;
    std::vector<float>    m_mcParContribToClustStartVy;
    std::vector<float>    m_mcParContribToClustStartVz;
    std::vector<float>    m_mcParContribToClustStopVx;
    std::vector<float>    m_mcParContribToClustStopVy;
    std::vector<float>    m_mcParContribToClustStopVz;
    std::vector<float>    m_mcParContribToClustTime;

    // becal cluster members (for cluster tree)
    uint64_t              m_numClustBECal;
    std::vector<uint64_t> m_becalClustNumCells;
    std::vector<float>    m_becalClustEne;
    std::vector<float>    m_becalClustEta;
    std::vector<float>    m_becalClustPhi;
    std::vector<float>    m_becalClustPosX;
    std::vector<float>    m_becalClustPosY;
    std::vector<float>    m_becalClustPosZ;
    std::vector<float>    m_becalClustTime;

    // utility members
    std::vector<double>    m_tileTilts;
    std::vector<double>    m_tileBarycenters;
    dd4hep::BitFieldCoder* m_decoder;

    // private methods
    void   InitializeDecoder();
    void   InitializeTrees();
    void   ResetVariables();
    double GetTime(const double tIn);

};

// end ------------------------------------------------------------------------
