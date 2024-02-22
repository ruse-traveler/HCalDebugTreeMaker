// ----------------------------------------------------------------------------
// 'JBarrelHCalTreeMakerProcessor.cc'
// Derek Anderson
// 04.13.2023
//
//
// A JANA plugin to construct a tree of BHCal tiles
// for training a ML clusterizer.
//
// Derived from code by Frederike Bock (thanks!!)
// ----------------------------------------------------------------------------

// user includes
#include "JBarrelHCalTreeMakerProcessor.h"

// the following just makes this a JANA plugin
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app -> Add(new JBarrelHCalTreeMakerProcessor);
  }
}



// inherited methods ----------------------------------------------------------

void JBarrelHCalTreeMakerProcessor::InitWithGlobalRootLock(){

  //  grab output file
  auto rootfile_svc = GetApplication() -> GetService<RootFile_service>();
  auto rootfile     = rootfile_svc     -> GetHistFile();

  // create directory
  m_dPluginDir = rootfile -> mkdir("JBarrelHCalTreeMaker");
  m_dPluginDir -> cd();

  // reset tree variables
  ResetVariables();

  // initialize output trees and tile maps
  InitializeDecoder();
  InitializeTrees();
  InitializeMaps();
  return;

}  // end 'InitWithGlobalRootLock()'



void JBarrelHCalTreeMakerProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  // reset tree variables
  ResetVariables();

  // loop over generated particles
  size_t nPar = 0;
  for (auto par : genParticles()) {

    // grab particle information
    const int   parType  = par -> getType();
    const float parEne   = par -> getEnergy();
    const float parPx    = par -> getMomentum().x;
    const float parPy    = par -> getMomentum().y;
    const float parPz    = par -> getMomentum().z;
    const float parPt    = std::sqrt((parPx * parPx) + (parPy * parPy));
    const float parEta   = -1. * std::log(std::tan(std::atan2(parPt, parPz) / 2.));
    const float parPhi   = std::atan2(parPy, parPx);

    // only accept truth particles
    const bool isTruth = (parType == 1);
    if (!isTruth) continue;

    // set output variables
    if (nPar < NMaxPars) {
      m_parEne[nPar] = parEne;
      m_parEta[nPar] = parEta;
      m_parPhi[nPar] = parPhi;
      ++nPar;
    }
  }  // end particle loop

  // loop over bhcal hits
  size_t nHit = 0;
  for (auto hit : bhcalRecHits()) {

    // grab hit info
    const uint64_t hitID = hit -> getCellID();
    const float    eHit  = hit -> getEnergy();

    // get hit indices
    const auto hitTile   = m_decoder -> get(hitID, 3);
    const auto hitTower  = m_decoder -> get(hitID, 2);
    const auto hitSector = m_decoder -> get(hitID, 1);

    // grab hit time
    const double maxTime = std::numeric_limits<double>::max();

    double hitTime = hit -> getTime();
    if (hitTime > maxTime) {
      hitTime = maxTime;
    }

    // add to vectors
    m_vecTileID.push_back(hitID);
    m_vecTileIsMatched.push_back(false);

    // set output variables
    m_tileEne[nHit]        = eHit;
    m_tileTime[nHit]       = hitTime;
    m_tileTilt[nHit]       = -1.;
    m_tileBarycenter[nHit] = -1.;
    m_tileCellID[nHit]     = hitID;
    m_tileTrueID[nHit]     = 0;  // FIXME this should be set to associated truth particle
    m_tileIndex[nHit]      = hitTile;
    m_tileTower[nHit]      = hitTower;
    m_tileSector[nHit]     = hitSector;
    ++nHit;
  }  // end bhcal hit loop

  // loop over bhcal clusters
  size_t nClustHCal = 0;
  for (auto bhClust : bhcalClusters()) {

    // grab cluster info
    const double nHitsClustHCal = bhClust -> getNhits();
    const double eClustHCal     = bhClust -> getEnergy();
    const double fClustHCal     = bhClust -> getIntrinsicPhi();
    const double tClustHCal     = bhClust -> getIntrinsicTheta();
    const double hClustHCal     = -1. * std::log(std::tan(tClustHCal / 2.));

    // associate each hit with corresponding cluster
    const auto bhcalClustHits = bhClust -> getHits();
    for (auto clustHit : bhcalClustHits) {

      // get hit ID
      const uint64_t clustHitID = clustHit.getCellID();

      // check if tile was hit
      size_t iAssocTile = -1;
      for (size_t iHit = 0; iHit < nHit; iHit++) {
        const bool isSameCell = (clustHitID == m_vecTileID.at(iHit));
        const bool wasMatched = m_vecTileIsMatched.at(iHit);
        if (isSameCell && !wasMatched) {
          m_tileClustIDA[iHit]     = nClustHCal;
          m_tileClustIDB[iHit]     = 0;
          m_vecTileIsMatched[iHit] = true;
        }
      }  // end hit tile loop
    }  // end cluster hit loop

    // set output variables
    m_bhcalClustNumCells[nClustHCal] = nHitsClustHCal;
    m_bhcalClustEne[nClustHCal]      = eClustHCal;
    m_bhcalClustEta[nClustHCal]      = hClustHCal;
    m_bhcalClustPhi[nClustHCal]      = fClustHCal;
    ++nClustHCal;
  }  // end bhcal cluster loop



  // fill becal cluster branches if needed
  size_t nClustECal = 0;
  if (AddBECalClusters) {

    // grab becal clusters and loop over them
    auto becalClusters = event -> Get<edm4eic::Cluster>(m_becalClustName.data());
    for (auto beClust : becalClusters) {

      // grab cluster info
      const double nHitsClustECal = beClust -> getNhits();
      const double eClustECal     = beClust -> getEnergy();
      const double fClustECal     = beClust -> getIntrinsicPhi();
      const double tClustECal     = beClust -> getIntrinsicTheta();
      const double hClustECal     = -1. * std::log(std::tan(tClustECal / 2.));

      // set output variables
      m_becalClustNumCells[nClustECal] = nHitsClustECal;
      m_becalClustEne[nClustECal]      = eClustECal;
      m_becalClustEta[nClustECal]      = hClustECal;
      m_becalClustPhi[nClustECal]      = fClustECal;
      ++nClustECal;
    }  // end becal cluster loop
  }  // end if (AddBECalClusters)

  // set output event variables
  m_numTiles      = nHit;
  m_numParticles  = nPar;
  m_numClustBHCal = nClustHCal;
  m_numClustBECal = nClustECal;

  // fill output trees
  m_tEventTree   -> Fill();
  m_tClusterTree -> Fill();
  return;

}  // end 'ProcessSequential(std::shared_ptr<JEvent>&)'



void JBarrelHCalTreeMakerProcessor::FinishWithGlobalRootLock() {

  // clean up variables
  ResetVariables();
  return;

}  // end 'FinishWithGlobalRootLock()'



// private methods ------------------------------------------------------------

void JBarrelHCalTreeMakerProcessor::InitializeDecoder() {

  // grab geometry service
  auto geom_svc = GetApplication() -> GetService<DD4hep_service>();

  // make sure readout is available
  dd4hep::IDDescriptor idDescriptor;
  try {
    idDescriptor  = geom_svc -> detector() -> readout("HcalBarrelHits").idSpec();
  } catch (...) {
    throw std::runtime_error("PANIC: readout class is not in output!");
  }

  // grab decoder
  short indexTile   = 0;
  short indexTower  = 0;
  short indexSector = 0;
  try {
    m_decoder   = geom_svc  -> detector() -> readout("HcalBarrelHits").idSpec().decoder();
    indexTile   = m_decoder -> index("tile");
    indexTower  = m_decoder -> index("tower");
    indexSector = m_decoder -> index("sector");
  } catch (...) {
    throw std::runtime_error("PANIC: something went wrong grabbing the decoder!");
  }
  return;

}  // end 'InitializeDecoder()'



void JBarrelHCalTreeMakerProcessor::InitializeTrees() {

  // switch to output directory
  m_dPluginDir -> cd();

  // initialize event tree
  m_tEventTree = new TTree("event_tree", "event_tree");
  m_tEventTree -> SetDirectory(m_dPluginDir);
  m_tEventTree -> Branch("cell_BHCAL_N",       &m_numTiles,       "cell_BHCAL_N/I");
  m_tEventTree -> Branch("cell_BHCAL_E",        m_tileEne,        "cell_BHCAL_E[cell_BHCAL_N]/F");
  m_tEventTree -> Branch("cell_BHCAL_T",        m_tileTime,       "cell_BHCAL_T[cell_BHCAL_N]/F");
  m_tEventTree -> Branch("cell_BHCAL_tilt",     m_tileTilt,       "cell_BHCAL_tilt[cell_BHCAL_N]/F");
  m_tEventTree -> Branch("cell_BHCAL_gravCent", m_tileBarycenter, "cell_BHCAL_gravCent[cell_BHCAL_N]/F");
  m_tEventTree -> Branch("cell_BHCAL_ID",       m_tileCellID,     "cell_BHCAL_ID[cell_BHCAL_N]/I");
  m_tEventTree -> Branch("cell_BHCAL_tile",     m_tileIndex,      "cell_BHCAL_tile[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_tower",    m_tileTower,      "cell_BHCAL_tower[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_sector",   m_tileSector,     "cell_BHCAL_sector[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_clusIDA",  m_tileClustIDA,   "cell_BHCAL_clusIDA[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_clusIDB",  m_tileClustIDB,   "cell_BHCAL_clusIDB[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_trueID",   m_tileTrueID,     "cell_BHCAL_trueID[cell_BHCAL_N]/I");

  // initialize cluster tree
  m_tClusterTree = new TTree("cluster_tree", "cluster_tree");
  m_tClusterTree -> SetDirectory(m_dPluginDir);
  m_tClusterTree -> Branch("mc_N",                 &m_numParticles,       "mc_N/I");
  m_tClusterTree -> Branch("mc_E",                  m_parEne,             "mc_E[mc_N]/F");
  m_tClusterTree -> Branch("mc_Phi",                m_parPhi,             "mc_Phi[mc_N]/F");
  m_tClusterTree -> Branch("mc_Eta",                m_parEta,             "mc_Eta[mc_N]/F");
  m_tClusterTree -> Branch("cluster_BHCAL_N",      &m_numClustBHCal,      "cluster_BHCAL_N/I");
  m_tClusterTree -> Branch("cluster_BHCAL_E",       m_bhcalClustEne,      "cluster_BHCAL_E[cluster_BHCAL_N]/F");
  m_tClusterTree -> Branch("cluster_BHCAL_Ncells",  m_bhcalClustNumCells, "cluster_BHCAL_Ncells[cluster_BHCAL_N]/I");
  m_tClusterTree -> Branch("cluster_BHCAL_Eta",     m_bhcalClustEta,      "cluster_BHCAL_Eta[cluster_BHCAL_N]/F");
  m_tClusterTree -> Branch("cluster_BHCAL_Phi",     m_bhcalClustPhi,      "cluster_BHCAL_Phi[cluster_BHCAL_N]/F");

  // add BECal branches if needed
  if (AddBECalClusters) {
    m_tClusterTree -> Branch("cluster_BECAL_N",      &m_numClustBECal,      "cluster_BECAL_N/I");
    m_tClusterTree -> Branch("cluster_BECAL_E",       m_becalClustEne,      "cluster_BECAL_E[cluster_BECAL_N]/F");
    m_tClusterTree -> Branch("cluster_BECAL_Ncells",  m_becalClustNumCells, "cluster_BECAL_Ncells[cluster_BECAL_N]/I");
    m_tClusterTree -> Branch("cluster_BECAL_Eta",     m_becalClustEta,      "cluster_BECAL_Eta[cluster_BECAL_N]/F");
    m_tClusterTree -> Branch("cluster_BECAL_Phi",     m_becalClustPhi,      "cluster_BECAL_Phi[cluster_BECAL_N]/F");
  }
  return;

} // end 'InitializeTrees()'



void JBarrelHCalTreeMakerProcessor::InitializeMaps() {

  /* initialize maps here */
  return;

}  // end 'InitializeMaps()'



void JBarrelHCalTreeMakerProcessor::ResetVariables() {

  // reset vectors of hit tiles
  m_vecTileID.clear();
  m_vecTileIsMatched.clear();

  // reset tile variables
  m_numTiles = 0;
  for (size_t iTile = 0; iTile < NTiles; iTile++) {
    m_tileEne[iTile]        = 0.;
    m_tileTime[iTile]       = 0.;
    m_tileTilt[iTile]       = 0.;
    m_tileBarycenter[iTile] = 0.;
    m_tileCellID[iTile]     = 0;
    m_tileTrueID[iTile]     = 0;
    m_tileIndex[iTile]      = 0;
    m_tileTower[iTile]      = 0;
    m_tileSector[iTile]     = 0;
    m_tileClustIDA[iTile]   = 0;
    m_tileClustIDB[iTile]   = 0;
  }

  // reset particle variables
  m_numParticles = 0;
  for (size_t iPar = 0; iPar < NMaxPars; iPar++) {
    m_parEne[iPar] = 0.;
    m_parPhi[iPar] = 0.;
    m_parEta[iPar] = 0.;
  }

  // reset bhcal/becal cluster variables
  m_numClustBHCal = 0;
  m_numClustBECal = 0;
  for (size_t iClust = 0; iClust < NMaxClust; iClust++) {
    m_bhcalClustNumCells[iClust] = 0;
    m_bhcalClustEne[iClust]      = 0.;
    m_bhcalClustEta[iClust]      = 0.;
    m_bhcalClustPhi[iClust]      = 0.;
    m_becalClustNumCells[iClust] = 0;
    m_becalClustEne[iClust]      = 0.;
    m_becalClustEta[iClust]      = 0.;
    m_becalClustPhi[iClust]      = 0.;
  }
  return;

}  // end 'ResetVariables()'

// end ------------------------------------------------------------------------
