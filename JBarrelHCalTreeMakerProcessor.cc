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

// root libraries
#include <Math/Vector4D.h>
// processor definition
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

  // reset tree/matching variables
  ResetVariables();

  // initialize event index
  m_evtIndex = 0;

  // initialize output trees
  InitializeDecoder();
  InitializeTrees();
  return;

}  // end 'InitWithGlobalRootLock()'



void JBarrelHCalTreeMakerProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

  // reset tree variables
  ResetVariables();

  // loop over generated particles
  size_t nGenPar = 0;
  for (auto gen : genParticles()) {

    // get particle kinematics
    ROOT::Math::PxPyPzE4D<float> pGen(
      gen -> getMomentum().x,
      gen -> getMomentum().y,
      gen -> getMomentum().z,
      gen -> getEnergy()
    );

    // only accept truth particles
    const int  typeGen = gen -> getType();
    const bool isTruth = (typeGen == 1);
    if (!isTruth) continue;

    // set output variables
    m_genParType.push_back( typeGen );
    m_genParPDG.push_back( gen -> getPDG() );
    m_genParEne.push_back( pGen.E() );
    m_genParEta.push_back( pGen.Eta() );
    m_genParPhi.push_back( pGen.Phi() );
    m_genParMass.push_back( gen -> getMass() );
    ++nGenPar;
  }  // end generated particle loop

  // loop over mc particles
  size_t nMCPar = 0;
  for (auto mc : mcParticles()) {

    // grab particle kinematics
    ROOT::Math::PxPyPzM4D<float> pMC(
      mc -> getMomentum().x,
      mc -> getMomentum().y,
      mc -> getMomentum().z,
      mc -> getMass()
    );

    // set output variables
    m_mcParGenStat.push_back( mc -> getGeneratorStatus() );
    m_mcParSimStat.push_back( mc -> getSimulatorStatus() );
    m_mcParPDG.push_back( mc -> getPDG() );
    m_mcParEne.push_back( pMC.E() );
    m_mcParPhi.push_back( pMC.Phi() );
    m_mcParEta.push_back( pMC.Eta() );
    m_mcParMass.push_back( mc -> getMass() );
    m_mcParStartVx.push_back( mc -> getVertex().x );
    m_mcParStartVx.push_back( mc -> getVertex().y );
    m_mcParStartVy.push_back( mc -> getVertex().z );
    m_mcParStopVx.push_back( mc -> getEndpoint().x );
    m_mcParStopVx.push_back( mc -> getEndpoint().y );
    m_mcParStopVy.push_back( mc -> getEndpoint().z );
    m_mcParTime.push_back( mc -> getTime() );
    ++nMCPar;
  }  // end mc particle loop

  // loop over bhcal hits
  size_t nHit = 0;
  for (auto hit : bhcalRecHits()) {

    // get hit indices
    const int64_t hitID     = hit -> getCellID();
    const auto    hitTile   = m_decoder -> get(hitID, 3);
    const auto    hitTower  = m_decoder -> get(hitID, 2);
    const auto    hitSector = m_decoder -> get(hitID, 1);

    // grab hit time
    const double maxTime = std::numeric_limits<double>::max();

    double hitTime = hit -> getTime();
    if (hitTime > maxTime) {
      hitTime = maxTime;
    }

    // add to maps
    m_mapTileIndexToID.insert( std::pair<int64_t, int64_t>(nHit, hitID) );
    m_mapTileIDToIsMatched.insert( std::pair<int64_t, bool>(hitID, false) );
    m_mapTileIDToClustIndex.insert( std::pair<int64_t, int64_t>(hitID, -1) );

    // set output variables
    m_tileEne.push_back( hit -> getEnergy() );
    m_tilePosX.push_back( hit -> getPosition().x );
    m_tilePosY.push_back( hit -> getPosition().y );
    m_tilePosZ.push_back( hit -> getPosition().z );
    m_tileTime.push_back( hitTime );
    m_tileCellID.push_back( hitID );
    m_tileTrueID.push_back( -1 );
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

    // calculate eta
    const double thClustHCal = bhClust -> getIntrinsicTheta();
    const double hClustHCal  = -1. * std::log(std::tan(thClustHCal / 2.));

    // grab consituent hits
    const auto bhcalClustHits = bhClust -> getHits();

    // associate each hit, contributing particle with corresponding cluster
    int64_t nHitsInClust = 0;
    for (auto clustHit : bhcalClustHits) {

      // get hit ID
      const uint64_t clustHitID = clustHit.getCellID();

      // set max possible time
      const double maxTime = std::numeric_limits<double>::max();

      // grab hit time
      double hitTime = clustHit.getTime();
      if (hitTime > maxTime) {
        hitTime = maxTime;
      }

      // check if tile was hit
      bool    isMatched  = false;
      int64_t iClustTile = -1;
      for (size_t iHit = 0; iHit < nHit; iHit++) {
        const bool isSameCell = (clustHitID == m_mapTileIndexToID[iHit]);
        if (isSameCell && !m_mapTileIDToIsMatched[clustHitID]) {
          isMatched  = true;
          iClustTile = iHit;
          break;
        }
      }  // end hit tile loop

      // if matched, set values appropriately
      if (isMatched) {
        m_tileClustIDA[iClustTile]          = nClustHCal;
        m_tileClustIDB[iClustTile]          = 0;
        m_mapTileIDToIsMatched[clustHitID]  = isMatched;
        m_mapTileIDToClustIndex[clustHitID] = nClustHCal;
      }

      // get hit indices
      const auto hitTile   = m_decoder -> get(clustHitID, 3);
      const auto hitTower  = m_decoder -> get(clustHitID, 2);
      const auto hitSector = m_decoder -> get(clustHitID, 1);

      /* loop over sim hits here */

      // set output variables
      m_tileInClustEne.push_back( clustHit.getEnergy() );
      m_tileInClustPosX.push_back( clustHit.getPosition().x );
      m_tileInClustPosY.push_back( clustHit.getPosition().y );
      m_tileInClustPosZ.push_back( clustHit.getPosition().z );
      m_tileInClustTime.push_back( hitTime );
      m_tileInClustCellID.push_back( clustHitID );
      m_tileInClustTrueID.push_back( -1 );
      m_tileInClustIndex.push_back( hitTile );
      m_tileInClustTower.push_back( hitTower );
      m_tileInClustSector.push_back( hitSector );
      m_tileInClustClustIDA.push_back( nClustHCal );
      m_tileInClustClustIDB.push_back( 0 );
      ++nHitsInClust;
    }  // end cluster hit loop

    // set output variables
    m_bhcalClustIndex.push_back( nClustHCal );
    m_bhcalClustNumCells.push_back( bhClust -> getNhits() );
    m_bhcalClustNumTileInClust.push_back( nHitsInClust );
    m_bhcalClustEne.push_back( bhClust -> getEnergy() );
    m_bhcalClustEta.push_back( hClustHCal );
    m_bhcalClustPhi.push_back( bhClust -> getIntrinsicPhi() );
    m_bhcalClustPosX.push_back( bhClust -> getPosition().x );
    m_bhcalClustPosY.push_back( bhClust -> getPosition().y );
    m_bhcalClustPosZ.push_back( bhClust -> getPosition().z );
    m_bhcalClustTime.push_back( bhClust -> getTime() );
    ++nClustHCal;
  }  // end bhcal cluster loop

  // fill becal cluster branches if needed
  size_t nClustECal = 0;
  if (AddBECalClusters) {

    // grab becal clusters and loop over them
    auto becalClusters = event -> Get<edm4eic::Cluster>(m_becalClustName.data());
    for (auto beClust : becalClusters) {

      // calculate eta
      const double thClustECal = beClust -> getIntrinsicTheta();
      const double hClustECal  = -1. * std::log(std::tan(thClustECal / 2.));

      // set output variables
      m_becalClustNumCells.push_back( beClust -> getNhits() );
      m_becalClustEne.push_back( beClust -> getEnergy() );
      m_becalClustEta.push_back( hClustECal );
      m_becalClustPhi.push_back( beClust -> getIntrinsicPhi() );
      m_becalClustPosX.push_back( beClust -> getPosition().x );
      m_becalClustPosY.push_back( beClust -> getPosition().y );
      m_becalClustPosZ.push_back( beClust -> getPosition().z );
      ++nClustECal;
    }  // end becal cluster loop
  }  // end if (AddBECalClusters)

  // set output event variables
  m_numTiles        = nHit;
  m_numGenParticles = nGenPar;
  m_numMCParticles  = nMCPar;
  m_numClustBHCal   = nClustHCal;
  m_numClustBECal   = nClustECal;

  // fill output trees
  m_tEventTree   -> Fill();
  m_tClusterTree -> Fill();

  // increment event index and exit routine
  ++m_evtIndex;
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
  m_tEventTree = new TTree("FlatTree", "Flat Tree with Clusters, Cells, and Particles");
  m_tEventTree -> SetDirectory(m_dPluginDir);
  m_tEventTree -> Branch("EvtIndex",      &m_evtIndex,        "EvtIndex/I");
  m_tEventTree -> Branch("EvtNHCell",     &m_numTiles,        "EvtNHCell/I");
  m_tEventTree -> Branch("EvtNGenPar",    &m_numGenParticles, "EvtNGenPar/I");
  m_tEventTree -> Branch("EvtNMCPar",     &m_numMCParticles,  "EvtNMCPar/I");
  m_tEventTree -> Branch("HCellEne",      &m_tileEne);
  m_tEventTree -> Branch("HCellPosX",     &m_tilePosX);
  m_tEventTree -> Branch("HCellPosY",     &m_tilePosY);
  m_tEventTree -> Branch("HCellPosZ",     &m_tilePosZ);
  m_tEventTree -> Branch("HCellTime",     &m_tileTime);
  m_tEventTree -> Branch("HCellID",       &m_tileCellID);
  m_tEventTree -> Branch("HCellTile",     &m_tileIndex);
  m_tEventTree -> Branch("HCellTower",    &m_tileTower);
  m_tEventTree -> Branch("HCellSector",   &m_tileSector);
  m_tEventTree -> Branch("HCellClusIDA",  &m_tileClustIDA);
  m_tEventTree -> Branch("HCellClusIDB",  &m_tileClustIDB);
  m_tEventTree -> Branch("HCellTrueID",   &m_tileTrueID);
  m_tEventTree -> Branch("HClustIndex",   &m_bhcalClustIndex);
  m_tEventTree -> Branch("HClustEne",     &m_bhcalClustEne);
  m_tEventTree -> Branch("HClustNCell",   &m_bhcalClustNumCells);
  m_tEventTree -> Branch("HClustEta",     &m_bhcalClustEta);
  m_tEventTree -> Branch("HClustPhi",     &m_bhcalClustPhi);
  m_tEventTree -> Branch("HClustPosX",    &m_bhcalClustPosX);
  m_tEventTree -> Branch("HClustPosY",    &m_bhcalClustPosY);
  m_tEventTree -> Branch("HClustPosZ",    &m_bhcalClustPosZ);
  m_tEventTree -> Branch("HClustTime",    &m_bhcalClustTime);
  m_tEventTree -> Branch("GenParType",    &m_genParType);
  m_tEventTree -> Branch("GenParPDG",     &m_genParPDG);
  m_tEventTree -> Branch("GenParE",       &m_genParEne);
  m_tEventTree -> Branch("GenParPhi",     &m_genParPhi);
  m_tEventTree -> Branch("GenParEta",     &m_genParEta);
  m_tEventTree -> Branch("GenParMass",    &m_genParMass);
  m_tEventTree -> Branch("MCParGenStat",  &m_mcParGenStat);
  m_tEventTree -> Branch("MCParSimStat",  &m_mcParSimStat);
  m_tEventTree -> Branch("MCParPDG",      &m_mcParPDG);
  m_tEventTree -> Branch("MCParE",        &m_mcParEne);
  m_tEventTree -> Branch("MCParPhi",      &m_mcParPhi);
  m_tEventTree -> Branch("MCParEta",      &m_mcParEta);
  m_tEventTree -> Branch("MCParMass",     &m_mcParMass);
  m_tEventTree -> Branch("MCParStartVx",  &m_mcParStartVx);
  m_tEventTree -> Branch("MCParStartVy",  &m_mcParStartVy);
  m_tEventTree -> Branch("MCParStartVz",  &m_mcParStartVz);
  m_tEventTree -> Branch("MCParStopVx",   &m_mcParStopVx);
  m_tEventTree -> Branch("MCParStopVy",   &m_mcParStopVy);
  m_tEventTree -> Branch("MCParStopVz",   &m_mcParStopVz);
  m_tEventTree -> Branch("MCParTime",     &m_mcParTime);

  // initialize cluster tree
  m_tClusterTree = new TTree("RelationalTree", "Tree with Cluster-Cell/Particle Relations");
  m_tClusterTree -> SetDirectory(m_dPluginDir);
  m_tClusterTree -> Branch("EvtIndex",                &m_evtIndex,      "EvtIndex/I");
  m_tClusterTree -> Branch("EvtNHClust",              &m_numClustBHCal, "EvtNHClust/I");
  m_tClusterTree -> Branch("HClustIndex",             &m_bhcalClustIndex);
  m_tClusterTree -> Branch("HClustEne",               &m_bhcalClustEne);
  m_tClusterTree -> Branch("HClustNCell",             &m_bhcalClustNumCells);
  m_tClusterTree -> Branch("HClustEta",               &m_bhcalClustEta);
  m_tClusterTree -> Branch("HClustPhi",               &m_bhcalClustPhi);
  m_tClusterTree -> Branch("HClustPosX",              &m_bhcalClustPosX);
  m_tClusterTree -> Branch("HClustPosY",              &m_bhcalClustPosY);
  m_tClusterTree -> Branch("HClustPosZ",              &m_bhcalClustPosZ);
  m_tClusterTree -> Branch("HClustTime",              &m_bhcalClustTime);
  m_tClusterTree -> Branch("HClustNTile",             &m_bhcalClustNumTileInClust);
  m_tClusterTree -> Branch("HClustNAssocPar",         &m_bhcalClustNumMCParAssocToClust);
  m_tClusterTree -> Branch("HClustNContribPar",       &m_bhcalClustNumMCParContribToClust);
  m_tClusterTree -> Branch("HClustCellClustIndex",    &m_tileInClustClustIndex);
  m_tClusterTree -> Branch("HClustCellEne",           &m_tileInClustEne);
  m_tClusterTree -> Branch("HClustCellPosX",          &m_tileInClustPosX);
  m_tClusterTree -> Branch("HClustCellPosY",          &m_tileInClustPosY);
  m_tClusterTree -> Branch("HClustCellPosZ",          &m_tileInClustPosZ);
  m_tClusterTree -> Branch("HClustCellTime",          &m_tileInClustTime);
  m_tClusterTree -> Branch("HClustCellTileCellID",    &m_tileInClustCellID);
  m_tClusterTree -> Branch("HClustCellTileTrueID",    &m_tileInClustTrueID);
  m_tClusterTree -> Branch("HClustCellTileTilID",     &m_tileInClustIndex);
  m_tClusterTree -> Branch("HClustCellTileTwrID",     &m_tileInClustTower);
  m_tClusterTree -> Branch("HClustCellTileSecID",     &m_tileInClustSector);
  m_tClusterTree -> Branch("HClustCellClustIDA",      &m_tileInClustClustIDA);
  m_tClusterTree -> Branch("HClustCellClustIDB",      &m_tileInClustClustIDB);
  m_tClusterTree -> Branch("HClustAssocClustIndex",   &m_mcParAssocToClustClustIndex);
  m_tClusterTree -> Branch("HClustAssocGenStat",      &m_mcParAssocToClustGenStat);
  m_tClusterTree -> Branch("HClustAssocSimStat",      &m_mcParAssocToClustSimStat);
  m_tClusterTree -> Branch("HClustAssocPDG",          &m_mcParAssocToClustPDG);
  m_tClusterTree -> Branch("HClustAssocEne",          &m_mcParAssocToClustEne);
  m_tClusterTree -> Branch("HClustAssocPhi",          &m_mcParAssocToClustPhi);
  m_tClusterTree -> Branch("HClustAssocEta",          &m_mcParAssocToClustEta);
  m_tClusterTree -> Branch("HClustAssocMass",         &m_mcParAssocToClustMass);
  m_tClusterTree -> Branch("HClustAssocStartVx",      &m_mcParAssocToClustStartVx);
  m_tClusterTree -> Branch("HClustAssocStartVy",      &m_mcParAssocToClustStartVy);
  m_tClusterTree -> Branch("HClustAssocStartVz",      &m_mcParAssocToClustStartVz);
  m_tClusterTree -> Branch("HClustAssocStopVx",       &m_mcParAssocToClustStopVx);
  m_tClusterTree -> Branch("HClustAssocStopVy",       &m_mcParAssocToClustStopVy);
  m_tClusterTree -> Branch("HClustAssocStopVz",       &m_mcParAssocToClustStopVz);
  m_tClusterTree -> Branch("HClustAssocTime",         &m_mcParAssocToClustTime);
  m_tClusterTree -> Branch("HClustContribClustIndex", &m_mcParContribToClustClustIndex);
  m_tClusterTree -> Branch("HClustContribGenStat",    &m_mcParContribToClustGenStat);
  m_tClusterTree -> Branch("HClustContribSimStat",    &m_mcParContribToClustSimStat);
  m_tClusterTree -> Branch("HClustContribPDG",        &m_mcParContribToClustPDG);
  m_tClusterTree -> Branch("HClustContribE",          &m_mcParContribToClustEne);
  m_tClusterTree -> Branch("HClustContribPhi",        &m_mcParContribToClustPhi);
  m_tClusterTree -> Branch("HClustContribEta",        &m_mcParContribToClustEta);
  m_tClusterTree -> Branch("HClustContribMass",       &m_mcParContribToClustMass);
  m_tClusterTree -> Branch("HClustContribStartVx",    &m_mcParContribToClustStartVx);
  m_tClusterTree -> Branch("HClustContribStartVy",    &m_mcParContribToClustStartVy);
  m_tClusterTree -> Branch("HClustContribStartVz",    &m_mcParContribToClustStartVz);
  m_tClusterTree -> Branch("HClustContribStopVx",     &m_mcParContribToClustStopVx);
  m_tClusterTree -> Branch("HClustContribStopVy",     &m_mcParContribToClustStopVy);
  m_tClusterTree -> Branch("HClustContribStopVz",     &m_mcParContribToClustStopVz);
  m_tClusterTree -> Branch("HClustContribTime",       &m_mcParContribToClustTime);

  // add BECal branches if needed
  if (AddBECalClusters) {

    // add to flat tree
    m_tEventTree -> Branch("EvtNEClust",  &m_numClustBECal, "EvtNEClust/I");
    m_tEventTree -> Branch("EClustEne",   &m_becalClustEne);
    m_tEventTree -> Branch("EClustNCell", &m_becalClustNumCells);
    m_tEventTree -> Branch("EClustEta",   &m_becalClustEta);
    m_tEventTree -> Branch("EClustPhi",   &m_becalClustPhi);
    m_tEventTree -> Branch("EClustPosX",  &m_becalClustPosX);
    m_tEventTree -> Branch("EClustPosY",  &m_becalClustPosY);
    m_tEventTree -> Branch("EClustPosZ",  &m_becalClustPosZ);
    m_tEventTree -> Branch("EClustTime",  &m_becalClustTime);

    // add to relational tree
    m_tClusterTree -> Branch("EvtNEClust",  &m_numClustBECal, "EvtNEClust/I");
    m_tClusterTree -> Branch("EClustEne",   &m_becalClustEne);
    m_tClusterTree -> Branch("EClustNCell", &m_becalClustNumCells);
    m_tClusterTree -> Branch("EClustEta",   &m_becalClustEta);
    m_tClusterTree -> Branch("EClustPhi",   &m_becalClustPhi);
    m_tClusterTree -> Branch("EClustPosX",  &m_becalClustPosX);
    m_tClusterTree -> Branch("EClustPosY",  &m_becalClustPosY);
    m_tClusterTree -> Branch("EClustPosZ",  &m_becalClustPosZ);
    m_tClusterTree -> Branch("EClustTime",  &m_becalClustTime);
  }
  return;

} // end 'InitializeTrees()'



void JBarrelHCalTreeMakerProcessor::ResetVariables() {

  // reset vectors of hit tiles [DEPRECATED]
  m_vecTileID.clear();
  m_vecTileIsMatched.clear();

  // reset maps
  m_mapTileIDToIsMatched.clear();
  m_mapTileIndexToID.clear();
  m_mapTileIDToClustIndex.clear();
  m_mapContribIDToTileID.clear();
  m_mapContribIDToClustIndex.clear();

  // reset event tile variables
  m_numTiles = 0;
  m_tileEne.clear();
  m_tilePosX.clear();
  m_tilePosY.clear();
  m_tilePosZ.clear();
  m_tileTime.clear();
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
  m_mcParGenStat.clear();
  m_mcParSimStat.clear();
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
  m_bhcalClustPosX.clear();
  m_bhcalClustPosY.clear();
  m_bhcalClustPosZ.clear();
  m_bhcalClustTime.clear();
  m_bhcalClustNumTileInClust.clear();
  m_bhcalClustNumMCParContribToClust.clear();
  m_bhcalClustNumMCParAssocToClust.clear();

  // reset bhcal tiles in cluster variables
  m_tileInClustClustIndex.clear();
  m_tileInClustEne.clear();
  m_tileInClustPosX.clear();
  m_tileInClustPosY.clear();
  m_tileInClustPosZ.clear();
  m_tileInClustTime.clear();
  m_tileInClustCellID.clear();
  m_tileInClustTrueID.clear();
  m_tileInClustIndex.clear();
  m_tileInClustTower.clear();
  m_tileInClustSector.clear();
  m_tileInClustClustIDA.clear();
  m_tileInClustClustIDB.clear();

  // reset mc particle associated to cluster variables
  m_mcParAssocToClustClustIndex.clear();
  m_mcParAssocToClustGenStat.clear();
  m_mcParAssocToClustSimStat.clear();
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
  m_mcParContribToClustGenStat.clear();
  m_mcParContribToClustSimStat.clear();
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
  m_becalClustPosX.clear();
  m_becalClustPosY.clear();
  m_becalClustPosZ.clear();
  m_becalClustTime.clear();
  return;

}  // end 'ResetVariables()'

// end ------------------------------------------------------------------------
