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
  gStyle->SetHistLineColor(kBlack);


  // read and print input file name  
  std::string outfile = argv[1];
  std::vector<std::string> filename;
  for (int i=2; i<argc; i++) { 
    std::cout << "FILE " << argv[i] << std::endl; 
    filename.push_back(argv[i]);
  }
 
  size_t nevents = 0;

  // floating point histograms 
  // TH2F::TH2F(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup);
  
  //Histograms about  parent particles
  //Note: these are initial values (can extend to final state values)
  TH1F* nparents = new TH1F("nparents", "NParents", 100, 0, 10);
  nparents->GetXaxis()->SetTitle("Num_Parents");
  nparents->GetYaxis()->SetTitle("Num_Events"); 
  TH1F* ppdg = new TH1F("ppdg", "Parent_PdgIDs", 200, 0, 20);
  ppdg->GetXaxis()->SetTitle("PdgCodes");
  ppdg->GetYaxis()->SetTitle("Num_Events");
  TH1F* px0 = new TH1F("px", "Parent_P_X_0", 100, -20, 20);
  px0->GetXaxis()->SetTitle("P_X_0");
  px0->GetYaxis()->SetTitle("Num_Events");
  TH1F* py0 = new TH1F("py", "Parent_P_Y_0", 100, -20, 20);
  py0->GetXaxis()->SetTitle("P_Y_0");
  py0->GetYaxis()->SetTitle("Num_Events");
  TH1F* p0 = new TH1F("p", "Parent_P_0", 100, -20, 20);
  p0->GetXaxis()->SetTitle("P_0");
  p0->GetYaxis()->SetTitle("Num_Events");
  TH1F* E = new TH1F("E", "Parent_E", 100, -20, 20);
  E->GetXaxis()->SetTitle("E");
  E->GetYaxis()->SetTitle("Num_Events");
  TH1F* dist = new TH1F("dist", "Parent_travel_dis", 100, 0, 20);
  dist->GetXaxis()->SetTitle("distance");
  dist->GetYaxis()->SetTitle("Num_Events"); 
  TH1F* theta_xy = new TH1F("theta_xy", "Parent_Theta_XY", 1000, -5, 5);
  theta_xy->GetXaxis()->SetTitle("theta_xy");
  theta_xy->GetYaxis()->SetTitle("Num_Events");

  /* 
  // Histograms about daughter particles 
  TH1F* ndaughter = new TH1F("ndaughter", "Ndaughters", 20, 0, 20);
  ndaughter->GetXaxis()->SetTitle("Num_daughters");
  ndaughter->GetYaxis()->SetTitle("Num_Events");
  TH1F* dpdg = new TH1F("dpdg", "Daughter_PDgIds", 200, 0, 20); 
  dpdg->GetXaxis()->SetTitle("daughter_PDgIds");
  dpdg->GetYaxis()->SetTitle("Num_Events");
  */ 

  // Correlation Plots  
  TH2F* ptheta = new TH2F("ptheta", "P vs Theta", 20, -1, 4, 20, 0, 20);
  ptheta->GetXaxis()->SetTitle("theta_xy");
  ptheta->GetYaxis()->SetTitle("P_0");
   
  // Clone() method makes a copy of any TObject 
  // TH2F* hB = (TH2F*) hA->Clone("hB");

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
     int num_par = mctruth->at(0).NParticles();
     nparents->Fill(num_par); 
     // for (auto const& particle : (*mctruth)) (this would be nice)
     for (int i=0; i<num_par; i++) 
	{
         //  extend this to the case where mctruth vector might have more 
         //  elements.
        
         // Get parent info 
         const simb::MCParticle& parent = mctruth->at(0).GetParticle(i); 
         int pdg = parent.PdgCode(); 
         // contains momentum as a vector
	 const TLorentzVector momentum = parent.Momentum();
         size_t p = parent.P();
         size_t momentum_x = momentum.Px();
         size_t momentum_y = momentum.Py(); 
	 size_t energy = momentum.E();
         size_t theta = TMath::ACos(momentum.CosTheta());  
	 size_t distance = parent.EndZ(); 
        
         // Fill corresponding parent histograms
         ppdg->Fill(pdg);
         p0->Fill(p);
         px0->Fill(momentum_x);
         py0->Fill(momentum_y);
         E->Fill(energy);
         dist->Fill(distance); 
         theta_xy->Fill(theta);
         int test = ptheta->Fill(theta, p);
         std::cout << "filling? " <<  test << std::endl;        
	
         /*
         // Get of daughters parents
         // doesn't quite work like that 
         auto const ldaughter = parent.fdaughters(); 
         int num_d = parent.NumberDaughters();
         ndaughter->Fill(num_d); 
         for (int j=0; j<num_d; j++)
	    { 
              const simb::MCParticle& daughter = ldaughter->at(j);
              int daughter_pdg = daughter.PdgCode();
              dpdg->Fill(daughter_pdg);
             }  
         */ 
	}  
     
  }
  // has to be DrawClone to survive after execution of macro1
  //ptheta->DrawClone(); 

// Save histogram (to file and PDF)
  TFile* fout = new TFile(outfile.c_str(), "recreate");

  // Write the hists and plots
  nparents->Write(); 
  ppdg->Write();
  p0->Write();
  px0->Write();
  py0->Write();
  E->Write(); 
  //ndaughter->Write();
  //dpdg->Write();   
  theta_xy->Write(); 
  ptheta->Write();
 

  fout->Close();

  return 0;
}
