
//
// Template for this file generated with eicmkplugin.py
//

#include "JBHCalMLClusterizerProcessor.h"
#include <services/rootfile/RootFile_service.h>


// The following just makes this a JANA plugin
extern "C" {
    void InitPlugin(JApplication *app) {
        InitJANAPlugin(app);
        app->Add(new JBHCalMLClusterizerProcessor);
    }
}

//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
void JBHCalMLClusterizerProcessor::InitWithGlobalRootLock(){
    auto rootfile_svc = GetApplication()->GetService<RootFile_service>();
    auto rootfile = rootfile_svc->GetHistFile();
    rootfile->mkdir("JBHCalMLClusterizer")->cd();

    // Create histograms here. e.g.
    // hEraw  = new TH1D("Eraw",  "BEMC hit energy (raw)",  100, 0, 0.075);
    // hEdigi = new TH2D("Edigi", "BEMC hit energy (digi) vs. raw", 200, 0, 2000.0, 100, 0, 0.075);
}

//-------------------------------------------
// ProcessSequential
//-------------------------------------------
void JBHCalMLClusterizerProcessor::ProcessSequential(const std::shared_ptr<const JEvent>& event) {

    // Fill histograms here. e.g.
    // for( auto hit : rawhits()  ) hEraw->Fill(  hit->getEnergy());
    // for( auto hit : digihits() ) hEdigi->Fill( hit->getAmplitude(), hit->getEnergy());
}

//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
void JBHCalMLClusterizerProcessor::FinishWithGlobalRootLock() {

    // Do any final calculations here.
}

