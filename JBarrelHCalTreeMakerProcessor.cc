// ----------------------------------------------------------------------------
// 'JBarreHCalTreeMakerProcessor.cc'
// Derek Anderson
// 04.13.2023
//
//
// Template for this file generated with eicmkplugin.py
// ----------------------------------------------------------------------------

// user includes
#include "JBarreHCalTreeMakerProcessor.h"
// JANA includes
#include <services/rootfile/RootFile_service.h>

using namespace std;

// the following just makes this a JANA plugin
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app -> Add(new JBarreHCalTreeMakerProcessor);
  }
}



// inherited methods ----------------------------------------------------------

void JBarreHCalTreeMakerProcessor::InitWithGlobalRootLock(){

  //  create directory in output file
  auto rootfile_svc = GetApplication() -> GetService<RootFile_service>();
  auto rootfile     = rootfile_svc     -> GetHistFile();
  rootfile -> mkdir("JBarreHCalTreeMaker") -> cd();

  // reset tree variables
  ResetVariables();

  // initialize output trees and tile maps
  InitializeDecoder();
  InitializeTrees();
  InitializeMaps();
  return;

}  // end 'InitWithGlobalRootLock()'



void JBarreHCalTreeMakerProcessor::ProcessSequential(const shared_ptr<const JEvent>& event) {

  // Fill histograms here. e.g.
  // for( auto hit : rawhits()  ) hEraw->Fill(  hit->getEnergy());
  // for( auto hit : digihits() ) hEdigi->Fill( hit->getAmplitude(), hit->getEnergy());

  return;

}  // end 'ProcessSequential(shared_ptr<JEvent>&)'



void JBarreHCalTreeMakerProcessor::FinishWithGlobalRootLock() {

  /* finalizing done here */
  return;

}  // end 'FinishWithGlobalRootLock()'



// private methods ------------------------------------------------------------

void JBarreHCalTreeMakerProcessor::InitializeDecoder() {

  dd4hep::Detector                     &detector  = dd4hep::Detector::getInstance();
  dd4hep::rec::CellIDPositionConverter  converter = new dd4hep::rec::CellIDPositionConverter(detector);
  try {

    // grab bhcal
    m_decoder = detector.readout("HcalBarrelHits").idSpec().decoder();

    // grab tiles
    auto tile   = m_decoder -> index("tile");
    auto tower  = m_decoder -> index("tower");
    auto sector = m_decoder -> index("sector");
    cout << "full list: " << " " << m_decoder -> fieldDescription() << endl;
  } catch (...) {
      cout <<"2nd: "  << m_decoder << endl;
      m_log- > error("PANIC: readout class not in the output");
      throw runtime_error("PANIC: readout class is not in the output!");
  }
  return;

}  // end 'InitializeDecoder()'



void JBarreHCalTreeMakerProcessor::InitializeTrees() {

  // initialize event tree
  m_tEventTree = new TTree("event_tree", "event_tree");
  m_tEventTree -> Branch("cell_BHCAL_N",       &m_numTiles,       "cell_BHCAL_N/I");
  m_tEventTree -> Branch("cell_BHCAL_E",        m_tileEne,        "cell_BHCAL_E[cell_BHCAL_N]/F");
  m_tEventTree -> Branch("cell_BHCAL_T",        m_tileTime,       "cell_BHCAL_T[cell_BHCAL_N]/F");
  m_tEventTree -> Branch("cell_BHCAL_tilt",     m_tileTilt,       "cell_BHCAL_tilt[cell_BHCAL_N]/F");
  m_tEventTree -> Branch("cell_BHCAL_gravCent", m_tileBarycenter, "cell_BHCAL_gravCent[cell_BHCAL_N]/F");
  m_tEventTree -> Branch("cell_BHCAL_tile",     m_tileIndex,      "cell_BHCAL_tile[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_tower",    m_tileTower,      "cell_BHCAL_tower[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_sector",   m_tileSector,     "cell_BHCAL_sector[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_clusIDA",  m_tileClustIDA,   "cell_BHCAL_clusIDA[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_clusIDB",  m_tileClustIDB,   "cell_BHCAL_clusIDB[cell_BHCAL_N]/S");
  m_tEventTree -> Branch("cell_BHCAL_trueID",   m_tileTrueID,     "cell_BHCAL_trueID[cell_BHCAL_N]/I");

  // initialize cluster tree
  m_tClusterTree -> Branch("mc_N",                 &m_numParticles,       "mc_N/I");
  m_tClusterTree -> Branch("mc_E",                  m_parEne,             "mc_E[mc_N]/F");
  m_tClusterTree -> Branch("mc_Phi",                m_parPhi,             "mc_Phi[mc_N]/F");
  m_tClusterTree -> Branch("mc_Eta",                m_parEta,             "mc_Eta[mc_N]/F");
  m_tClusterTree -> Branch("cluster_BHCAL_N",      &m_numClustBHCal,      "cluster_BHCAL_N/I");
  m_tClusterTree -> Branch("cluster_BHCAL_E",       m_bhcalClustEne,      "cluster_BHCAL_E[cluster_BHCAL_N]/F");
  m_tClusterTree -> Branch("cluster_BHCAL_Ncells",  m_bhcalClustNumCells, "cluster_BHCAL_Ncells[cluster_BHCAL_N]/I");
  m_tClusterTree -> Branch("cluster_BHCAL_Eta",     m_bhcalClustEta,      "cluster_BHCAL_Eta[cluster_BHCAL_N]/F");
  m_tClusterTree -> Branch("cluster_BHCAL_Phi",     m_bhcalCLustPhi,      "cluster_BHCAL_Phi[cluster_BHCAL_N]/F");
  m_tClusterTree -> Branch("cluster_BECAL_N",      &m_numClustBECal,      "cluster_BECAL_N/I");
  m_tClusterTree -> Branch("cluster_BECAL_E",       m_becalClustEne,      "cluster_BECAL_E[cluster_BECAL_N]/F");
  m_tClusterTree -> Branch("cluster_BECAL_Ncells",  m_becalClustNumCells, "cluster_BECAL_Ncells[cluster_BECAL_N]/I");
  m_tClusterTree -> Branch("cluster_BECAL_Eta",     m_becalClustEta,      "cluster_BECAL_Eta[cluster_BECAL_N]/F");
  m_tClusterTree -> Branch("cluster_BECAL_Phi",     m_becalClustPhi,      "cluster_BECAL_Phi[cluster_BECAL_N]/F");
  return;

} // end 'InitializeTrees()'



void JBarreHCalTreeMakerProcessor::InitializeMaps() {

  return;

}  // end 'InitializeMaps()'



void JBarreHCalTreeMakerProcessor::ResetVariables() {

  // reset tile variables
  m_numTiles = 0;
  for (size_t iTile = 0; iTile < NTiles; iTile++) {
    m_tileEne[iTile]        = 0.;
    m_tileTime[iTile]       = 0.;
    m_tileTilt[iTile]       = 0.;
    m_tileBarycenter[iTile] = 0.;
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
