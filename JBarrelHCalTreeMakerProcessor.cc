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
    m_parEne.push_back( parEne );
    m_parEta.push_back( parEta );
    m_parPhi.push_back( parPhi );
    ++nPar;
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
  m_tEventTree -> Branch("cell_BHCAL_N",        &m_numTiles, "cell_BHCAL_N/I");
  m_tEventTree -> Branch("cell_BHCAL_E",        &m_tileEne);
  m_tEventTree -> Branch("cell_BHCAL_T",        &m_tileTime);
  m_tEventTree -> Branch("cell_BHCAL_tilt",     &m_tileTilt);
  m_tEventTree -> Branch("cell_BHCAL_gravCent", &m_tileBarycenter);
  m_tEventTree -> Branch("cell_BHCAL_ID",       &m_tileCellID);
  m_tEventTree -> Branch("cell_BHCAL_tile",     &m_tileIndex);
  m_tEventTree -> Branch("cell_BHCAL_tower",    &m_tileTower);
  m_tEventTree -> Branch("cell_BHCAL_sector",   &m_tileSector);
  m_tEventTree -> Branch("cell_BHCAL_clusIDA",  &m_tileClustIDA);
  m_tEventTree -> Branch("cell_BHCAL_clusIDB",  &m_tileClustIDB);
  m_tEventTree -> Branch("cell_BHCAL_trueID",   &m_tileTrueID);

  // initialize cluster tree
  m_tClusterTree = new TTree("cluster_tree", "cluster_tree");
  m_tClusterTree -> SetDirectory(m_dPluginDir);
  m_tClusterTree -> Branch("mc_N",                 &m_numParticles, "mc_N/I");
  m_tClusterTree -> Branch("mc_E",                 &m_parEne);
  m_tClusterTree -> Branch("mc_Phi",               &m_parPhi);
  m_tClusterTree -> Branch("mc_Eta",               &m_parEta);
  m_tClusterTree -> Branch("cluster_BHCAL_N",      &m_numClustBHCal);
  m_tClusterTree -> Branch("cluster_BHCAL_E",      &m_bhcalClustEne);
  m_tClusterTree -> Branch("cluster_BHCAL_Ncells", &m_bhcalClustNumCells);
  m_tClusterTree -> Branch("cluster_BHCAL_Eta",    &m_bhcalClustEta);
  m_tClusterTree -> Branch("cluster_BHCAL_Phi",    &m_bhcalClustPhi);

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

  // reset tile variables
  m_numTiles = 0;
  m_tileEne.clear();
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

  // reset particle variables
  m_numParticles = 0;
  m_parEne.clear();
  m_parPhi.clear();
  m_parEta.clear();

  // reset bhcal/becal cluster variables
  m_numClustBHCal = 0;
  m_numClustBECal = 0;
  m_bhcalClustNumCells.clear();
  m_bhcalClustEne.clear();
  m_bhcalClustEta.clear();
  m_bhcalClustPhi.clear();
  m_becalClustNumCells.clear();
  m_becalClustEne.clear();
  m_becalClustEta.clear();
  m_becalClustPhi.clear();
  return;

}  // end 'ResetVariables()'



void JBarrelHCalTreeMakerProcessor::ResizeEvtTiles(const size_t nTiles) {

  m_tileEne.resize(nTiles);
  m_tileTime.resize(nTiles);
  m_tileTilt.resize(nTiles);
  m_tileBarycenter.resize(nTiles);
  m_tileCellID.resize(nTiles);
  m_tileTrueID.resize(nTiles);
  m_tileIndex.resize(nTiles);
  m_tileTower.resize(nTiles);
  m_tileSector.resize(nTiles);
  m_tileClustIDA.resize(nTiles);
  m_tileClustIDB.resize(nTiles);
  return;

}  // end 'ResizeEvtTiles(size_t)'



void JBarrelHCalTreeMakerProcessor::ResizeEvtPars(const size_t nPars) {

  m_parEne.resize(nPars);
  m_parPhi.resize(nPars);
  m_parEta.resize(nPars);
  return;

}  // end 'ResizeEvtPars(size_t)'



void JBarrelHCalTreeMakerProcessor::ResizeBHCalClusts(const size_t nBHCalClust) {

  m_bhcalClustNumCells.resize(nBHCalClust);
  m_bhcalClustEne.resize(nBHCalClust);
  m_bhcalClustEta.resize(nBHCalClust);
  m_bhcalClustPhi.resize(nBHCalClust);
  return;

}  // end 'ResizeBHCalClusts(size_t)'



void JBarrelHCalTreeMakerProcessor::ResizeBECalClusts(const size_t nBECalClust) {

  m_becalClustNumCells.resize(nBECalClust);
  m_becalClustEne.resize(nBECalClust);
  m_becalClustEta.resize(nBECalClust);
  m_becalClustPhi.resize(nBECalClust);
  return;

}  // end 'ResizeBECalClusts(size_t)'

// end ------------------------------------------------------------------------
