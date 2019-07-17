/* Analyze
 *
 *  * A. Mastbaum <mastbaum@uchicago.edu>, 2019/04/22
 *  * Edited: D. Belayneh <dawit@uchicago.edu>, 2019/12/07
 *
 */

#include <iostream>
#include <map>
#include <string>
#include <vector>

// ROOT classes
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TMath.h"
#include "TStyle.h"
// gallery classes 
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
// LarSoft classes 
#include "lardataobj/Simulation/SimPhotons.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"


int main(int argc, char* argv[]) {
  // Parse command-line arguments
   if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " "
              << "OUTPUT.root INPUT.root [INPUT2.root ...]" << std::endl;
    return 0;
  }

  // graph formatting lines 
  gStyle->SetOptStat(0);
  gStyle->SetHistLineColor(kBlue);
  gStyle->SetFillColor(0);   
  gStyle->SetOptStat("nemr"); 

  // read and print input file name  
  std::string outfile = argv[1];
  std::vector<std::string> filename;
  for (int i=2; i<argc; i++) { 
    std::cout << "FILE " << argv[i] << std::endl; 
    filename.push_back(argv[i]);
  }
 
  size_t nevents = 0;

  //Histograms about  parent particles
  //Note: these are initial values (can extend to final state values)
  TH1F* nparts = new TH1F("nparts", "NParticles", 1000, 0, 10);
  nparts->GetXaxis()->SetTitle("Num_Particles"); 
  TH1F* ppdg = new TH1F("ppdg", "Parent_PdgIDs", 1000, 0, 20);
  ppdg->GetXaxis()->SetTitle("PdgCodes");
  TH1F* px0 = new TH1F("px", "Parent_P_X_0", 1000, -2, 2);
  px0->GetXaxis()->SetTitle("P_X_0 [Gev]");
  TH1F* py0 = new TH1F("py", "Parent_P_Y_0", 1000, -2, 2);
  py0->GetXaxis()->SetTitle("P_Y_0 [Gev]");
  TH1F* pz0 = new TH1F("pz", "Parent_P_Z_0", 1000, 0, 3);
  pz0->GetXaxis()->SetTitle("P_Z_0 [Gev]");
  TH1F* p0 = new TH1F("p", "Parent_P_0", 1000, 0, 3);
  p0->GetXaxis()->SetTitle("P_0 [Gev]");
  TH1F* E = new TH1F("E", "Parent_E", 1000, 0, 3);
  E->GetXaxis()->SetTitle("E [Gev]");
  //TH1F* dist = new TH1F("dist", "Parent_travel_dis", 100, 0, 20);
  //dist->GetXaxis()->SetTitle("distance");
  TH1F* costheta_xy = new TH1F("costheta_xy", "Parent_CosTheta", 10000, -1, 1);
  costheta_xy->GetXaxis()->SetTitle("Costheta");

  // Correlation Plots  
  TH2F* pcostheta = new TH2F("pcostheta", "P vs CosTheta", 20, 0.8, 1, 20, 0, 2);
  pcostheta->GetXaxis()->SetTitle("CosTheta");
  pcostheta->GetYaxis()->SetTitle("P_0 [Gev]");
  pcostheta->SetOption("COLZ"); // draws scatter
  // Clone() method makes a copy of any TObject 
  // TH2F* hB = (TH2F*) hA->Clone("hB");


  // Daughter Lepton Histograms 
  TH1F* lpdg = new TH1F("lpdg", "Lepton_PdgID", 2000, 0, 20);
  lpdg->GetXaxis()->SetTitle("PdgCodes");
  TH1F* lpx0 = new TH1F("lpx", "Lepton_P_X_0", 1000, -2, 2);
  lpx0->GetXaxis()->SetTitle("P_X_0 [Gev]");
  TH1F* lpy0 = new TH1F("lpy", "Lepton_P_Y_0", 1000, -2, 2);
  lpy0->GetXaxis()->SetTitle("P_Y_0 [Gev]");
  TH1F* lpz0 = new TH1F("lpz", "Lepton_P_Z_0", 1000, -2, 2);
  lpz0->GetXaxis()->SetTitle("P_Z_0 [Gev]");
  TH1F* lp0 = new TH1F("lp", "Lepton_P_0", 1000, 0, 3);
  lp0->GetXaxis()->SetTitle("P_0 [Gev]");
  TH1F* lE = new TH1F("lE", "Lepton_E", 1000, 0, 3);
  lE->GetXaxis()->SetTitle("E [Gev]");
  TH1F* lcostheta_xy = new TH1F("lcostheta_xy", "Lepton_CosTheta", 1000, -1, 1);
  lcostheta_xy->GetXaxis()->SetTitle("Costheta");

  TH2F* lpcostheta = new TH2F("lpcostheta", "Lepton P vs CosTheta", 75, -1, 1, 75, 0, 2);
  lpcostheta->GetXaxis()->SetTitle("CosTheta");
  lpcostheta->GetYaxis()->SetTitle("P_0 [Gev]");
  lpcostheta->SetOption("COLZ");  

// Event loop
// get events from filename untill there are none
for (gallery::Event ev(filename) ; !ev.atEnd(); ev.next()) {
    if (nevents % 100 == 0) {
      std::cout << "EVENT " << nevents << std::endl;
    }
    // increment event counter 
    nevents++;

    // events that go into hA 
    gallery::Handle<std::vector<simb::MCTruth> > mctruth;
    ev.getByLabel({"generator"}, mctruth);

    // events that go into hB (don't need second Handle)
    // gallery::Handle<std::vector<sim::SimPhotons> > photonsB;
    // ev.getByLabel({"photonpropagation","","PhotonLibraryPropagation"}, photonsB);

    //assert(photonsA->size() == photonsB->size());

    /*
    for (size_t i=0; i<photonsA->size(); i++) {
      const sim::SimPhotons& pA = photonsA->at(i);
      const sim::SimPhotons& pB = photonsB->at(i);
      std::cout << i << ": " << pA.size() << " " << pB.size() << std::endl;
      hA->Fill(pA.OpChannel(), pA.size());
      hB->Fill(pB.OpChannel(), pB.size());
    }
    */ 

     // Loop through particles in MCTruth 
     // for (auto const& particle : (*mctruth)) (this would be nice)
     for (int i=0; i<1; i++) // was i < num_par (loop not used for now)
	{

         // Get parent info 
         const simb::MCNeutrino& neutrino = mctruth->at(0).GetNeutrino();
         const simb::MCParticle& parent = neutrino.Nu(); 
         int pdg = parent.PdgCode(); 
         // contains momentum as a vector
	 const TLorentzVector momentum = parent.Momentum();
         float p = parent.P();
         float momentum_x = momentum.Px();
         float momentum_y = momentum.Py(); 
         float momentum_z = momentum.Pz(); 
	 float energy = momentum.E();
         float costheta = momentum.CosTheta();  
	 //size_t distance = parent.EndZ(); 
        
         // Fill corresponding parent histograms 
         ppdg->Fill(pdg);
         p0->Fill(p);
         px0->Fill(momentum_x);
         py0->Fill(momentum_y);
         pz0->Fill(momentum_z); 
         E->Fill(energy);
         //dist->Fill(distance); 
         costheta_xy->Fill(costheta);
         pcostheta->Fill(costheta, p);     
	
         // Get daughter Lepton (Electron or Muon)
         const simb::MCParticle& lep = neutrino.Lepton(); 
         int lep_pdg = lep.PdgCode(); 
	 const TLorentzVector lmomentum = lep.Momentum();
         float lp = lep.P(); 
         float lmomentum_x = lmomentum.Px();
         float lmomentum_y = lmomentum.Py();
         float lmomentum_z = lmomentum.Pz();
         float lenergy = lmomentum.E();
         float lcostheta = lmomentum.CosTheta();  
         // Fill daughter  Lepton hists
         lpdg->Fill(lep_pdg);
         lp0->Fill(lp);
         lpx0->Fill(lmomentum_x);
         lpy0->Fill(lmomentum_y);
         lpz0->Fill(lmomentum_z);
         lE->Fill(lenergy);
         lcostheta_xy->Fill(lcostheta); 
         lpcostheta->Fill(lcostheta, lp);
        }
}

         // Draw Hists
         ppdg->Draw();
         gPad->SaveAs("ParentID.png");
         p0->Draw();
         gPad->SaveAs("P0.png");
         px0->Draw();
         gPad->SaveAs("Px0.png");
         py0->Draw();
         gPad->SaveAs("Py0.png");
         pz0->Draw(); 
         gPad->SaveAs("Pz0.png");
         E->Draw();
         gPad->SaveAs("ParentE.png");
         //dist->Fill(distance); 
         costheta_xy->Draw();
         gPad->SaveAs("Parent_CosTheta.png");
         pcostheta->Draw();
         gPad->SaveAs("Parent_PvsCosTheta.png");     
	
         lpdg->Draw();
         gPad->SaveAs("LeptonID.png");
         lp0->Draw();
         gPad->SaveAs("LeptonP0.png");
         lpx0->Draw();
         gPad->SaveAs("LeptonPx.png");
         lpy0->Draw();
         gPad->SaveAs("LeptonPy.png");
         lpz0->Draw();
         gPad->SaveAs("LeptonPz.png");
         lE->Draw();
         gPad->SaveAs("LeptonE.png");
         lcostheta_xy->Draw();
         gPad->SaveAs("LeptonCosTheta.png"); 
         lpcostheta->Draw();
         gPad->SaveAs("LeptonPvsCosTheta.png");
        
// Save histogram (to file and PDF)
  TFile* fout = new TFile(outfile.c_str(), "recreate");

  // Write the hists and plots 
  ppdg->Write();
  p0->Write();
  px0->Write();
  py0->Write();
  pz0->Write();
  E->Write(); 
  costheta_xy->Write(); 
  pcostheta->Write();
 
  lpdg->Write(); 
  lp0->Write();
  lpx0->Write();
  lpy0->Write();
  lpz0->Write();
  lE->Write(); 
  lcostheta_xy->Write(); 
  lpcostheta->Write();


  fout->Close();

  return 0;
}
