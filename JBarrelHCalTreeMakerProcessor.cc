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

  // initialize event index
  m_evtIndex = 0;

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
  size_t nGenPar = 0;
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
    m_genParEne.push_back( parEne );
    m_genParEta.push_back( parEta );
    m_genParPhi.push_back( parPhi );
    ++nGenPar;
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
    m_vecTileID.push_back( hitID );
    m_vecTileIsMatched.push_back( false );

    // set output variables
    m_tileEne.push_back( eHit );
    m_tileTime.push_back( hitTime );
    m_tileTilt.push_back( -1. );
    m_tileBarycenter.push_back( -1. );
    m_tileCellID.push_back( hitID );
    m_tileTrueID.push_back( 0 );  // FIXME this should be set to associated truth particle
    m_tileIndex.push_back( hitTile );
    m_tileTower.push_back( hitTower );
    m_tileSector.push_back( hitSector );
    ++nHit;
  }  // end bhcal hit loop

  // resize remaining tile vectors
  m_tileClustIDA.resize(nHit);
  m_tileClustIDB.resize(nHit);

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
    m_bhcalClustNumCells.push_back( nHitsClustHCal );
    m_bhcalClustEne.push_back( eClustHCal );
    m_bhcalClustEta.push_back( hClustHCal );
    m_bhcalClustPhi.push_back( fClustHCal );
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
      m_becalClustNumCells.push_back( nHitsClustECal );
      m_becalClustEne.push_back( eClustECal );
      m_becalClustEta.push_back( hClustECal );
      m_becalClustPhi.push_back( fClustECal );
      ++nClustECal;
    }  // end becal cluster loop
  }  // end if (AddBECalClusters)

  // set output event variables
  m_numTiles        = nHit;
  m_numGenParticles = nGenPar;
  m_numClustBHCal   = nClustHCal;
  m_numClustBECal   = nClustECal;

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
  m_tEventTree -> Branch("evt_index",           &m_evtIndex,        "evt_index/I");
  m_tEventTree -> Branch("cell_BHCAL_N",        &m_numTiles,        "cell_BHCAL_N/I");
  m_tEventTree -> Branch("gen_N",               &m_numGenParticles, "gen_N/I");
  m_tEventTree -> Branch("mc_N",                &m_numMCParticles,  "mc_N/I");
  m_tEventTree -> Branch("cell_BHCAL_E",        &m_tileEne);
  m_tEventTree -> Branch("cell_BHCAL_PosX",     &m_tilePosX);
  m_tEventTree -> Branch("cell_BHCAL_PosY",     &m_tilePosY);
  m_tEventTree -> Branch("cell_BHCAL_PosZ",     &m_tilePosZ);
  m_tEventTree -> Branch("cell_BHCAL_Time",     &m_tileTime);
  m_tEventTree -> Branch("cell_BHCAL_Tilt",     &m_tileTilt);
  m_tEventTree -> Branch("cell_BHCAL_GravCent", &m_tileBarycenter);
  m_tEventTree -> Branch("cell_BHCAL_ID",       &m_tileCellID);
  m_tEventTree -> Branch("cell_BHCAL_Tile",     &m_tileIndex);
  m_tEventTree -> Branch("cell_BHCAL_Tower",    &m_tileTower);
  m_tEventTree -> Branch("cell_BHCAL_Sector",   &m_tileSector);
  m_tEventTree -> Branch("cell_BHCAL_ClusIDA",  &m_tileClustIDA);
  m_tEventTree -> Branch("cell_BHCAL_ClusIDB",  &m_tileClustIDB);
  m_tEventTree -> Branch("cell_BHCAL_TrueID",   &m_tileTrueID);
  m_tEventTree -> Branch("gen_Type",            &m_genParType);
  m_tEventTree -> Branch("gen_PDG",             &m_genParPDG);
  m_tEventTree -> Branch("gen_E",               &m_genParEne);
  m_tEventTree -> Branch("gen_Phi",             &m_genParPhi);
  m_tEventTree -> Branch("gen_Eta",             &m_genParEta);
  m_tEventTree -> Branch("gen_Mass",            &m_genParMass);
  m_tEventTree -> Branch("mc_GenStatus",        &m_mcParGenStatus);
  m_tEventTree -> Branch("mc_SimStatus",        &m_mcParSimStatus);
  m_tEventTree -> Branch("mc_PDG",              &m_mcParPDG);
  m_tEventTree -> Branch("mc_E",                &m_mcParEne);
  m_tEventTree -> Branch("mc_Phi",              &m_mcParPhi);
  m_tEventTree -> Branch("mc_Eta",              &m_mcParEta);
  m_tEventTree -> Branch("mc_Mass",             &m_mcParMass);
  m_tEventTree -> Branch("mc_StartVx",          &m_mcParStartVx);
  m_tEventTree -> Branch("mc_StartVy",          &m_mcParStartVy);
  m_tEventTree -> Branch("mc_StartVz",          &m_mcParStartVz);
  m_tEventTree -> Branch("mc_StopVx",           &m_mcParStopVx);
  m_tEventTree -> Branch("mc_StopVy",           &m_mcParStopVy);
  m_tEventTree -> Branch("mc_StopVz",           &m_mcParStopVz);
  m_tEventTree -> Branch("mc_Time",             &m_mcParTime);

  // initialize cluster tree
  m_tClusterTree = new TTree("cluster_tree", "cluster_tree");
  m_tClusterTree -> SetDirectory(m_dPluginDir);
  m_tClusterTree -> Branch("evt_index",                    &m_evtIndex,      "evt_index/I");
  m_tClusterTree -> Branch("cluster_BHCAL_N",              &m_numClustBHCal, "cluster_BHCAL_N/I");
  m_tClusterTree -> Branch("cluster_BHCAL_index",          &m_bhcalClustIndex);
  m_tClusterTree -> Branch("cluster_BHCAL_E",              &m_bhcalClustEne);
  m_tClusterTree -> Branch("cluster_BHCAL_NCells",         &m_bhcalClustNumCells);
  m_tClusterTree -> Branch("cluster_BHCAL_Eta",            &m_bhcalClustEta);
  m_tClusterTree -> Branch("cluster_BHCAL_Phi",            &m_bhcalClustPhi);
  m_tClusterTree -> Branch("cluster_BHCAL_NTile",          &m_numTileInClust);
  m_tClusterTree -> Branch("cluster_BHCAL_NAssoc",         &m_numMCParAssocToClust);
  m_tClusterTree -> Branch("cluster_BHCAL_NContrib",       &m_numMCParContribToClust);
  m_tClusterTree -> Branch("cluster_BHCAL_TileClustIndex", &m_tileInClustClustIndex);
  m_tClusterTree -> Branch("cluster_BHCAL_TileE",          &m_tileInClustEne);
  m_tClusterTree -> Branch("cluster_BHCAL_TilePosX",       &m_tileInClustPosX);
  m_tClusterTree -> Branch("cluster_BHCAL_TilePosY",       &m_tileInClustPosY);
  m_tClusterTree -> Branch("cluster_BHCAL_TilePosZ",       &m_tileInClustPosZ);
  m_tClusterTree -> Branch("cluster_BHCAL_TileTime",       &m_tileInClustTime);
  m_tClusterTree -> Branch("cluster_BHCAL_TileTilt",       &m_tileInClustTilt);
  m_tClusterTree -> Branch("cluster_BHCAL_TileGravCent",   &m_tileInClustBarycenter);
  m_tClusterTree -> Branch("cluster_BHCAL_TileCellID",     &m_tileInClustCellID);
  m_tClusterTree -> Branch("cluster_BHCAL_TileTrueID",     &m_tileInClustTrueID);
  m_tClusterTree -> Branch("cluster_BHCAL_TileTilID",      &m_tileInClustIndex);
  m_tClusterTree -> Branch("cluster_BHCAL_TileTwrID",      &m_tileInClustTower);
  m_tClusterTree -> Branch("cluster_BHCAL_TileSecID",      &m_tileInClustSector);
  m_tClusterTree -> Branch("cluster_BHCAL_TileClustIDA",   &m_tileInClustClustIDA);
  m_tClusterTree -> Branch("cluster_BHCAL_TileClustIDB",   &m_tileInClustClustIDB);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocClustIndex", &m_mcParAssocToClustClustIndex);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocGenStat",   &m_mcParAssocToClustGenStatus);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocSimStat",   &m_mcParAssocToClustSimStatus);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocPDG",       &m_mcParAssocToClustPDG);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocE",         &m_mcParAssocToClustEne);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocPhi",       &m_mcParAssocToClustPhi);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocEta",       &m_mcParAssocToClustEta);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocMass",      &m_mcParAssocToClustMass);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocStartVx",   &m_mcParAssocToClustStartVx);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocStartVy",   &m_mcParAssocToClustStartVy);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocStartVz",   &m_mcParAssocToClustStartVz);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocStopVx",    &m_mcParAssocToClustStopVx);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocStopVy",    &m_mcParAssocToClustStopVy);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocStopVz",    &m_mcParAssocToClustStopVz);
  m_tClusterTree -> Branch("cluster_BHCAL_AssocTime",      &m_mcParAssocToClustTime);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribClustIndex", &m_mcParContribToClustClustIndex);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribGenStat", &m_mcParContribToClustGenStatus);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribSimStat", &m_mcParContribToClustSimStatus);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribPDG",     &m_mcParContribToClustPDG);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribE",       &m_mcParContribToClustEne);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribPhi",     &m_mcParContribToClustPhi);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribEta",     &m_mcParContribToClustEta);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribMass",    &m_mcParContribToClustMass);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribStartVx", &m_mcParContribToClustStartVx);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribStartVy", &m_mcParContribToClustStartVy);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribStartVz", &m_mcParContribToClustStartVz);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribStopVx",  &m_mcParContribToClustStopVx);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribStopVy",  &m_mcParContribToClustStopVy);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribStopVz",  &m_mcParContribToClustStopVz);
  m_tClusterTree -> Branch("cluster_BHCAL_ContribTime",    &m_mcParContribToClustTime);


  // add BECal branches if needed
  if (AddBECalClusters) {
    m_tClusterTree -> Branch("cluster_BECAL_N",      &m_numClustBECal, "cluster_BECAL_N/I");
    m_tClusterTree -> Branch("cluster_BECAL_E",      &m_becalClustEne);
    m_tClusterTree -> Branch("cluster_BECAL_Ncells", &m_becalClustNumCells);
    m_tClusterTree -> Branch("cluster_BECAL_Eta",    &m_becalClustEta);
    m_tClusterTree -> Branch("cluster_BECAL_Phi",    &m_becalClustPhi);
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

  // reset event tile variables
  m_numTiles = 0;
  m_tileEne.clear();
  m_tilePosX.clear();
  m_tilePosY.clear();
  m_tilePosZ.clear();
  m_tileTime.clear();
  m_tileTilt.clear();
  m_tileBarycenter.clear();
  m_tileCellID.clear();
  m_tileTrueID.clear();
  m_tileIndex.clear();
  m_tileTower.clear();
  m_tileSector.clear();
  m_tileClustIDA.clear();
  m_tileClustIDB.clear();

  // reset event gen particle variables
  m_numGenParticles = 0;
  m_genParType.clear();
  m_genParPDG.clear();
  m_genParEne.clear();
  m_genParPhi.clear();
  m_genParEta.clear();
  m_genParMass.clear();

  // reset event mc particle variables
  m_numMCParticles = 0;
  m_mcParGenStatus.clear();
  m_mcParSimStatus.clear();
  m_mcParPDG.clear();
  m_mcParEne.clear();
  m_mcParPhi.clear();
  m_mcParEta.clear();
  m_mcParMass.clear();
  m_mcParStartVx.clear();
  m_mcParStartVy.clear();
  m_mcParStartVz.clear();
  m_mcParStopVx.clear();
  m_mcParStopVy.clear();
  m_mcParStopVz.clear();
  m_mcParTime.clear();

  // reset bhcal cluster variables
  m_numClustBHCal = 0;
  m_numClustBECal = 0;
  m_bhcalClustNumCells.clear();
  m_bhcalClustEne.clear();
  m_bhcalClustEta.clear();
  m_bhcalClustPhi.clear();
  m_numTileInClust.clear();
  m_numMCParContribToClust.clear();
  m_numMCParAssocToClust.clear();

  // reset bhcal tiles in cluster variables
  m_tileInClustClustIndex.clear();
  m_tileInClustEne.clear();
  m_tileInClustPosX.clear();
  m_tileInClustPosY.clear();
  m_tileInClustPosZ.clear();
  m_tileInClustTime.clear();
  m_tileInClustTilt.clear();
  m_tileInClustBarycenter.clear();
  m_tileInClustCellID.clear();
  m_tileInClustTrueID.clear();
  m_tileInClustIndex.clear();
  m_tileInClustTower.clear();
  m_tileInClustSector.clear();
  m_tileInClustClustIDA.clear();
  m_tileInClustClustIDB.clear();

  // reset mc particle associated to cluster variables
  m_mcParAssocToClustClustIndex.clear();
  m_mcParAssocToClustGenStatus.clear();
  m_mcParAssocToClustSimStatus.clear();
  m_mcParAssocToClustPDG.clear();
  m_mcParAssocToClustEne.clear();
  m_mcParAssocToClustPhi.clear();
  m_mcParAssocToClustEta.clear();
  m_mcParAssocToClustMass.clear();
  m_mcParAssocToClustStartVx.clear();
  m_mcParAssocToClustStartVy.clear();
  m_mcParAssocToClustStartVz.clear();
  m_mcParAssocToClustStopVx.clear();
  m_mcParAssocToClustStopVy.clear();
  m_mcParAssocToClustStopVz.clear();
  m_mcParAssocToClustTime.clear();

  // reset mc particle contributing to cluster variables
  m_mcParContribToClustClustIndex.clear();
  m_mcParContribToClustGenStatus.clear();
  m_mcParContribToClustSimStatus.clear();
  m_mcParContribToClustPDG.clear();
  m_mcParContribToClustEne.clear();
  m_mcParContribToClustPhi.clear();
  m_mcParContribToClustEta.clear();
  m_mcParContribToClustMass.clear();
  m_mcParContribToClustStartVx.clear();
  m_mcParContribToClustStartVy.clear();
  m_mcParContribToClustStartVz.clear();
  m_mcParContribToClustStopVx.clear();
  m_mcParContribToClustStopVy.clear();
  m_mcParContribToClustStopVz.clear();
  m_mcParContribToClustTime.clear();

  // reset becal cluster variables
  m_becalClustNumCells.clear();
  m_becalClustEne.clear();
  m_becalClustEta.clear();
  m_becalClustPhi.clear();
  return;

}  // end 'ResetVariables()'

// end ------------------------------------------------------------------------
