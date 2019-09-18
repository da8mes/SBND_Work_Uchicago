/* Analyze_g4
 *
 *  * Written by: D. Belayneh <dawit@uchicago.edu>, 2019/12/10
 *
 *
 */

#include <iostream>
#include <map>
#include <string>
#include <vector>

// ROOT classes
#include "TCanvas.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TMath.h"
#include "TStyle.h"
// gallery classes 
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
//#include "gallery/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
// LarSoft classes 
#include "lardataobj/Simulation/SimPhotons.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"

//===============GLOBAL VARIABLES==================================// 

// Detector Dimensions: the outer most co-ordinate limits for the active 
// volumes for TPC0 and TPC1
// width (x)
const double x_low0 = -199; 
const double x_high0 = -2.65;
const double x_low1 = 2.65; 
const double x_high1 = 199.15; 
// height (y)
const double y_low0 = -200; 
const double y_high0 = 200;
const double y_low1 = -200;
const double y_high1 = 200; 
// length (z)
const double z_low0 = 0; 
const double z_high0 = 500;
const double z_low1 = 0;
const double z_high1 = 500;
// ==============BIN DEFS =========================================//
//===============END OF GLOBAL VARIABLES===========================// 
 
// checks if the given position is in the detector (SBND)
// TPCs. returns 1 if true, 0 if false. 
int inTPC(TLorentzVector pos)
{
	double x = pos[0];
	double y = pos[1];
	double z = pos[2];
 
	if (x < x_low0 || x > x_high1) return 0; 
	else if (y < y_low0 || y > y_high1) return 0;
	else if (z < z_low0 || z > z_high1) return 0; 
	else return 1;
}

// populate muEoActiveLen given
// MCTrajectory of the Particle and the corresponding 
// two hists for track length and initial Energy 
void fillEnergyLength(const simb::MCTrajectory& track, TH2F* hist, TH1F* test, TH1F* test1, TH1F* test2)
{

	int muonsInTPC = 0; 

	simb::MCTrajectory::const_iterator pt; 
	TLorentzVector tPos;
	TLorentzVector tMom;
	TLorentzVector tPosInit; // first active point 
	TLorentzVector tPosFinal; // last active point 

	double Eo; 
	int step = 1; // move along track via step size  
	int eInitUndetermined = 1; // flag to find Eo 
	for (pt = track.begin(); pt != track.end(); pt+=step)
	{ 
		tPos = (*pt).first;
		int intpc = inTPC(tPos); 
		if (intpc)
		{
			if (eInitUndetermined)	
			{
				// first tpc point
				tPosInit = tPos;
				test->Fill(tPos[1]); // tPos-y
				tMom = (*pt).second;
				Eo = tMom.E(); // init energy
				// test1->Fill(Eo);
				muonsInTPC++;
				eInitUndetermined = 0; // init point is determined
				continue; 	
			}
	
			else // never init point 
			{
				TLorentzVector npt = (*(pt+step)).first; 
				if (inTPC(npt)) continue; 
				// last point in tpc
				tPosFinal = tPos;
				double activeDist = (tPosFinal.Vect() - tPosInit.Vect()).Mag();
				// test2->Fill(activeDist); 
				hist->Fill(Eo, activeDist); 
				return ; 				
			}	
		}
	}
	return ; 
}

// write hist to TNtuple 
void makeTNtupleFromHist(TNtuple* table, TH2F* hist) 
{

	int binx; 
	int biny; 
	int binz; 
        int nBinx = hist->GetNbinsX();
        int nBiny = hist->GetNbinsY();
	int nBin = hist->GetBin(nBinx, nBiny);
        for (int i = 1; i <= nBin; i++)
        {
                hist->GetBinXYZ(i, binx, biny, binz);
                double N = hist->GetBinContent(i);
		double Eo = hist->GetXaxis()->GetBinCenter(binx);  
		double dist = hist->GetYaxis()->GetBinCenter(biny);  
 	        	
		std::cout << "i: " << i << std::endl; 		
		std::cout << "binx: " << binx << std::endl; 
		std::cout << "biny: " << biny << std::endl; 		
		std::cout << "binz: " << binz << std::endl; 		
		std::cout << "N: " << N << std::endl; 		
		std::cout << "Eo: " << Eo << std::endl; 		
		std::cout << "dist: " << dist << std::endl; 		
		std::cout << "\n ------------ \n " << dist << std::endl; 				 
		 
		table->Fill(N, Eo, dist);  
	} 
	return ;
}


// Calculate the Point of Closest Approach for pos from parent 
// track. Returns the shortest distance between pos and the entire 
// trajectory track. 
double poca(const simb::MCTrajectory& track, TLorentzVector pos)
{
	simb::MCTrajectory::const_iterator pt;
	pt = track.begin(); 
	TLorentzVector tPos = (*pt).first; 
	// advance on the track via step
	int step = 1; // controls our sensitivity to deltas 
	double distMin = (pos.Vect() - tPos.Vect()).Mag(); 
	for ( ; pt != track.end(); pt+=step)
	{ 
		tPos = (*pt).first;
		if (not inTPC(tPos)) continue; 
		double dist = (pos.Vect() - tPos.Vect()).Mag();       
	        //std::cout << "\n----\ndist: " << dist << std::endl;
		if (dist < distMin) distMin = dist; 
	}        
        //std::cout << "\n----\npoca: " << distMin << std::endl;
	return distMin; 
}

int main(int argc, char* argv[]) {
  // Parse command-line arguments
   if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " "
              << "OUTPUT.root INPUT.root [INPUT2.root ...]" << std::endl;
    return 0;
  }

  // graph formatting lines 
  gStyle->SetHistLineColor(kBlue);
  //gStyle->SetFillColor(0);   
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(1111);
  gROOT->ForceStyle(); 

  // read and print input file name  
  std::string outfile = argv[1];
  std::vector<std::string> filename;
  for (int i=2; i<argc; i++) { 
    std::cout << "FILE " << argv[i] << std::endl; 
    filename.push_back(argv[i]);
  }


  
//=============Constants==========================================//
  // Ar pdgCodes 
  //const int Ar_40 = 1000180400;
  //const int Ar_39 = 1000180390; 
  size_t nevents = 0;
  //int nElectrons = 0;
  int nMuonsMinus = 0;
  int nMuonsPlus = 0;  
  int nMuonElectrons = 0; 
  int nMichelElectrons = 0;
  int nDeltaElectrons = 0;
  int nBremElectrons = 0; 


//================ Test Hists  ========================================//
 
  TH1F* TEST_HIST = new TH1F("testy", "y-tPosInit (expect line at 200)", 1000, -300, 300); 
  TH1F* TEST_HIST1 = new TH1F("testEo", "Eo-tPosInit", 100, 0, 100); 
  TH1F* TEST_HIST2 = new TH1F("testActDist", "activeDist (expect line at 400)", 1000, 0, 500); 

//================= Muon Hists=======================================//

  // variable hist for muEoActiveDist energy [0, 1]: 250 Mev ; [1, 30]: 1Gev ; [30,100]: 10 Gev
  const Int_t EMAX = 100; // Gev
  const Int_t XMAX = 500; // Cm 
  
  const Int_t ENBINS = 40 + (30 - 1) + (EMAX - 30)/10;    
  Double_t step = 0.025; 
  Double_t Eedges[ENBINS+1]; 
  int i = 0; 
  Double_t e = 0;  
  for (int i = 0; i < ENBINS+1; i++)
     {
	Eedges[i] = e;
	if (e >= 1 && e < 30) step = 1; 
        if (e >= 30) step = 10;
        e+= step;  
     }
 
  const Int_t XNBINS = 100; 
  Double_t Xedges[XNBINS+1]; 
  for (int i = 0; i < XNBINS+1; i++)
     {
         Xedges[i] = i * (XMAX/XNBINS);  
     }
 
  TH2F* muEoActiveDist = new TH2F("muEoAD", "mu (Eo, ActiveDistance)", ENBINS, Eedges, XNBINS, Xedges); 
  muEoActiveDist->SetXTitle("Eo [Gev]");
  muEoActiveDist->SetYTitle("dist [cm]");
  muEoActiveDist->SetOption("COLZ"); 

  TH2F* muEoActiveDistStop = new TH2F("muEoADstop", "muStop (Eo, ActiveDistance)", 100, 0, 10, 100, 0, 500); 
  muEoActiveDistStop->SetXTitle("Eo [Gev]");
  muEoActiveDistStop->SetYTitle("dist [cm]");
  muEoActiveDistStop->SetOption("COLZ");

  TH2F* muEoActiveDistPen = new TH2F("muEoADpen", "muPen (Eo, ActiveDistance)", 100, 0, 10, 100, 0, 500); 
  muEoActiveDistPen->SetXTitle("Eo [Gev]");
  muEoActiveDistPen->SetYTitle("dist [cm]");
  muEoActiveDistPen->SetOption("COLZ");

  // All Muons; make this TH2F* of (Eo, activeDist)
  TF1* fMuEInit = new TF1("fMuEInit", "[0]*exp(-x/[1])", 0, 100);
  TH1F* muEInit = new TH1F("mu_energy_init", "mu_Energy_Initial", 100, 0, 100); 
  TF1* fMuLen = new TF1("fMuLen", "[0]*exp(-x/[1])", 0, 3000); // might not be exponential
  // max Active TPC Length is 4^2 + 4^2 + 5^2 == 57 [m] == 5700 [cm] 
  TH1F* muActiveLen = new TH1F("mu_len", "mu_ActiveTPC_TLength", 60, 0, 6000); 

  // Penetrating Muons
  TF1* fMuEInitPen = new TF1("fMuEInitPen", "[0]*exp(-x/[1])", 0, 100);
  TH1F* muEInitPen = new TH1F("mu_energy_pen", "mu_Energy_Penetrating", 50, 0, 100); 
  TF1* fMuLenPen = new TF1("fMuLenPen", "[0]*exp(-x/[1])", 0, 3000);
  TH1F* muActiveLenPen = new TH1F("mu_len_pen", "mu_ActiveTPC_TLength_Penetrating", 50, 0, 100);

  // Stopping Muons
  TF1* fMuEInitStop = new TF1("fMuEInitStop", "[0]*exp(-x/[1])", 0, 100);
  TH1F* muEInitStop = new TH1F("mu_energy_stop", "mu_Energy_Stopping", 50, 0, 2); 
  TF1* fMuLenStop = new TF1("fMuLenStop", "[0]*exp(-x/[1])", 0, 3000);
  TH1F* muActiveLenStop = new TH1F("mu_len_stop", "mu_ActiveTPC_TLength_Stopping", 50, 0, 100);

//==============================================//

  TH1F* muEnergy = new TH1F("mu_energy", "mu_Energy", 100, 0, 100); 
  TH1F* muBeta = new TH1F("mu_beta", "mu_Beta", 100, 0.8, 1.1); 
  TH1F* muGamma = new TH1F("mu_gamma", "mu_Gamma", 100, 1, 200);
  muGamma->SetXTitle("gamma"); 
  // define an exponential for lifetime fit f(x)=A*exp(-x/tau) 
  TF1* fMinusLife = new TF1("fMinusLife", "[0]*exp(-x/[1])", 0, 3000);
  TF1* fPlusLife = new TF1("fPlusLife", "[0]*exp(-x/[1])", 0, 3000);
  fMinusLife->SetParameters(10, 2000);
  fPlusLife->SetParameters(10, 2000); 
  TH1F* muMinusLife = new TH1F("mu_minus_life", "Muon- Lifetime [ns]", 500, 0, 3000);
  muMinusLife->SetFillColor(kBlue); 
  muMinusLife->SetXTitle("t [ns]");  
  TH1F* muPlusLife = new TH1F("mu_plus_life", "Muon+ Lifetime [ns]", 500, 0, 3000);
  muPlusLife->SetFillColor(kBlue); 
  muPlusLife->SetXTitle("t [ns]"); 

//=================Electron Hists=====================================//
 
  // All Electron E
  TH1F* eEnergy = new TH1F("eE", "e_Energy", 1000, 0, 1);
  
  eEnergy->SetFillColor(kBlue);
  // All Muon E
  TH1F* eEnergyMuon = new TH1F("eE_mu", "e_Energy_mu", 500, 0.001, 0.5); 
  eEnergyMuon->SetFillColor(kRed);
  eEnergyMuon->SetXTitle("e Energy [Gev]"); 

  // Stack to compare MC delta vs Theory delta
  THStack* deltaStack = new THStack("deltaStack", "Delta Probability Dist");  
  // Log Binning
  // bins = {10^-4, ..., 10^-1}, 100 log spaced bins  
  const Int_t NBINS = 100;
  Double_t edges[NBINS + 1];  
  for (int i = 0; i < NBINS+1; i++) 
     {
        edges[i] = pow(10, -2.0 + 2.0*i/NBINS); 	
     }

  TH1F* eEnergyDelta = new TH1F("eE_mu_delta", "e_Energy_delta", NBINS, edges);
  eEnergyDelta->SetFillColor(kGreen);
  eEnergyDelta->SetXTitle("e Energy [Gev]"); 

  TH1F* eEnergyBrem = new TH1F("eE_mu_brem", "e_Energy_Brem", 1000, 0, 2);
  eEnergyBrem->SetFillColor(kMagenta);
  eEnergyBrem->SetXTitle("e Energy [Gev]"); 

  TH1F* eEnergyPP = new TH1F("eE_mu_pp", "e_Energy_Brem_PairProd", 1000, 0, 2);
  eEnergyBrem->SetFillColor(kBlue);
  eEnergyBrem->SetXTitle("e Energy [Gev]"); 

  //Theoretical Expectation of Delta Energy Spectrum
  //Identical binning with eEnergyDelta    
  //TH1F* eEnergyDeltaTheory = new TH1F("eE_theory_delta", "e_Energy_Delta_Theory", NBINS, edges);
  //eEnergyDelta->SetFillColor(kDeepSea);

  // Michelle Electron energies 
  TH1F* eEnergyMich = new TH1F("e_michE", "e_Michel_Energy", 100, 0, 0.06);
  eEnergyMich->SetFillColor(kBlue); 
  eEnergyMich->SetXTitle("e Energy [Gev]");  
  
  TH1F* eCosDelta = new TH1F("e_delta_costheta", "e_CosTheta_Delta", 20, -1, 1); 
  TH1F* eDistMu = new TH1F("e_dist_mu", "Electron Distance from Muon End Pos", 1000, -1, 10); 
  eDistMu->SetFillColor(kDeepSea);
  eDistMu->SetXTitle("dist [cm]");  

  // electron Point of Closest Approcach from Muon track 
  TH1F* ePocaMu = new TH1F("e_poca_mu", "Electron POCA Muon", 1000, -0.01, 0.5);
  ePocaMu->SetFillColor(kBlue); 
  ePocaMu->SetXTitle("dist [cm]");  
  TH1F* eTrackDist = new TH1F("e_Track_Dist", "e_Total_TrackLength [no Michels]", 1000, 0, 14); 
  eTrackDist->SetFillColor(kBlue);
  eTrackDist->SetXTitle("dist [cm]");  

//=============Event loop=================================================//
// get events from filename untill there are none
for (gallery::Event ev(filename) ; !ev.atEnd(); ev.next()) {
    if (nevents % 100 == 0) {
      //std::cout << "EVENT " << nevents << std::endl;
    }
    // increment event counter 
    nevents++; 
     
    gallery::Handle<std::vector<simb::MCTruth> > mctruth;                       
    ev.getByLabel({"generator"}, mctruth);

    // get simulated particle tracks from g4
    gallery::Handle<std::vector<sim::MCTrack> > mctrackHandle;
    ev.getByLabel({"largeant"}, mctrackHandle); 

    // get particles in MCParticles from g4 (also contains all genie info) 
    gallery::Handle<std::vector<simb::MCParticle>> mcparticleHandle;
    ev.getByLabel({"largeant"}, mcparticleHandle);
    
    // This is a map between a particle's trackId and
    // its MCParticle pointer i.e has to be used as ->  
    std::map<int, const simb::MCParticle* > particleMap; 
    
    // Loop through MCParticle and create the map 
    for (auto const& particle : (*mcparticleHandle))
      {

        int particleTrackId = particle.TrackId();
        // Add the address of the MCParticle to the map, with the 
        // trackId as the key.
        particleMap[particleTrackId] = &particle;

              
        /*
	std::cout << "Track ID: " << particle.TrackId() << std::endl;
	std::cout << "Start Process: " << particle.Process() << std::endl; 
	std::cout << "End Process: " << particle.EndProcess() << std::endl; 
	std::cout << "Parent track ID: " << particle.Mother() << std::endl; 
	std::cout << "Start energy: " << particle.E() << " GeV"<< std::endl;
        */  
      }
    
    for (auto const& particle: (*mcparticleHandle)) 
      { 

        // select for u- (pdgId = 13)
        int pdgId = particle.PdgCode(); 
        if (abs(pdgId) != 13)
            continue; 
  
	// Hist mu Gamma and Beta dists 
	// Using cosmic muon start momentum 
	// there should be a better way of extracting 
	// beta and gamma
	const simb::MCTrajectory& muTrajectory = particle.Trajectory();	
	const size_t mu_last = particle.NumberTrajectoryPoints() - 1; 
	//const TLorentzVector& muPositionStart = particle.Position(0);  
	const TLorentzVector& muPositionEnd = particle.Position(mu_last);  
        const TVector3& muPositionEnd3 = muPositionEnd.Vect(); 
        const TLorentzVector& muP = particle.Momentum(0); 
  
	if (pdgId == 13) nMuonsMinus++; 
	if (pdgId == -13) nMuonsPlus++; 
	muEnergy->Fill(muP.E());  
	muBeta->Fill(muP.Beta());
	muGamma->Fill(muP.Gamma());
	// a global look at active track dist and muEInit 
	// fillEnergyLength(muTrajectory, muEoActiveDist, TEST_HIST, TEST_HIST1, TEST_HIST2);
        // Populate hists for Penetrating and Stopping Muons
        /* 
        if (inTPC(muPositionEnd)) // stopping 
           {
		fillEnergyLength(muTrajectory, muEInitStop, muActiveLenStop); 
           }
	else // penetrating 
	   {
		fillEnergyLength(muTrajectory, muEInitPen, muActiveLenPen); 
	   }
	*/ 
	// continue;  
        // To populate expected delta events for each muon, use a curried version of 
        // nDeltaEvents 

        //=======================Daughters of Mother=======================// 
        for (int i = 0 ; i < particle.NumberDaughters() ; i++) // access daughter trackIds 
           {
             const int dId = particle.Daughter(i); // daughter ID 
             auto search = particleMap.find(dId); // * to <int, MCParticle> pair 
             if (search == particleMap.end() ) // didn't find daughter in list 
               continue; 
                   
             // access daughter particle via search
             const simb::MCParticle& daughter = *(search->second);
             
             if (abs(daughter.PdgCode()) != 11) // abs necessary b/c u+ decay  
                  continue; // daughter not electron.
 
             const auto dProcess = daughter.Process(); 
	     //std::cout << "\n--\ne Parent Start Process: " << dProcess << std::endl; 
             //const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints(); 
             //const int last = numberTrajectoryPoints - 1;                    
             const TLorentzVector& positionStart = daughter.Position(0); 
             const TVector3& positionStart3 = positionStart.Vect(); 
             //const TLorentzVector& positionEnd = daughter.Position(last); 
             const TLorentzVector& momentumStart = daughter.Momentum(0); 
             //const TLorentzVector& momentumEnd = daughter.Momentum(last); 
             float eE = momentumStart.E();

	     const simb::MCTrajectory& dTrajectory = daughter.Trajectory();
	     
	     float dTrackLength = dTrajectory.TotalLength(); // track length for e 
	     // starting dist of all e w.r.t final mu position
	     // but ignoring time has the problem of picking up 
	     // wildily out of time but spatially coincident
	     
	     // are muons tracked after leaving the tpcs? yes. 
             float e_dist_mu = (muPositionEnd3 - positionStart3).Mag();  
	     //eDistMu->Fill(e_dist_mu); 
 	     nMuonElectrons++;  
             eEnergyMuon->Fill(eE);

	     // We can differentiate between drift and delta 
	     // by: track length (i.e reconstructable electrons)
	     // they have to move at least 3 wire hits (9mm) 
	     // which means they need at least that much energy 
	     // if they don't have that energy they are drift.
	     // That Minimum Energy == 2 [Mev] 
	     // Citation: https://physics.nist.gov/cgi-bin/Star/e_table.pl
	    	     
	     //Apply minimum energy cut (2Mev)
	     // or trackLength == 9 mm == 0.9 cm
	     // anything that travels less can not be reconstructed.
	     //
	     // Note: energy cut and dist length are not identical(?) 
	     //if (dTrackLength < 0.9) continue;
	     if (eE < 0.002) continue; // get rid of drifts 

	     float epsilon = 0.2; // [cm] from eDist plot
	     // this cut throws away tip deltas 
 	     if (e_dist_mu < epsilon) // michel 
               {
		nMichelElectrons++; 
		eEnergyMich->Fill(eE);

                // calculate muon life time
                // muLife ~= 2.19 [us] 
                // Then fit an exponential to it and find tau = 2.19 [us]
                float time_to_decay = positionStart.T() - muPositionEnd.T();
		// because decay time are different: 
		if (pdgId == 13)
		{
		    muMinusLife->Fill(time_to_decay);	
		}
		else
		{
		    muPlusLife->Fill(time_to_decay);  
 		} 
		continue; // get rid of michels  
               }

	     // plot POCA
	     // then apply a cut on pca ; no michels or drifts
	     // if what remains is just deltas then poca for all e == 0.
	     // if spread, then there is something else contaminating the
	     // delta hists   
	     double pca = poca(muTrajectory, positionStart); 
	     // ePocaMu->Fill(pca);

	     // Apply poca cut found from ePocaMU 
	     // pca <= 0.001 cm
	     // Now only delta's should survive
	     if (pca > 0.001)
		{
		 //std::cout << "Start Process for Poca < 0.02: " << dProcess << std::endl;  
		 nBremElectrons++; 
		 // eEnergyBrem->Fill(eE); 
		 if (dProcess == "muPairProd") eEnergyPP->Fill(eE); 
		 continue; 
		}
 
	     float cosTheta = momentumStart.CosTheta(); 
	     nDeltaElectrons++;
	     eCosDelta->Fill(cosTheta);
	     eTrackDist->Fill(dTrackLength);  
	     eEnergyDelta->Fill(eE);     

           } // finish u- daughter loop       
        } // finish mcparticle loop  
      } // finish event loop  

  /*
  // Perform various fits to hists
  muEInit->Fit("fMuEInit");
  muActiveLen->Fit("fMuLen"); 

  muEInitPen->Fit("fMuEInitPen");
  muActiveLenPen->Fit("fMuLenPen");
  
  muEInitStop->Fit("fMuEInitStop");
  muActiveLenStop->Fit("fMuLenStop");  
  */

  // Add to stackDelta
  deltaStack->Add(eEnergyDelta); 
 
  // Normalize delta plot to integral == 1
  // and re-draw
  Double_t scale = 1/eEnergyDelta->Integral();
  eEnergyDelta->Scale(scale);
  eEnergyDelta->DrawCopy(); 
 
  muMinusLife->Fit("fMinusLife");   
  muPlusLife->Fit("fPlusLife");  

  std::cout << "Total MuonsMinus: " << nMuonsMinus << std::endl;  
  std::cout << "Total MuonsPlus: " << nMuonsPlus << std::endl;  
  std::cout << "Total Muons: " << nMuonsMinus + nMuonsPlus << std::endl;  
  std::cout << "Total Muon Electrons: " << nMuonElectrons << std::endl;
  std::cout << "Total Michel Electrons: " << nMichelElectrons << std::endl;
  std::cout << "Total Delta Electrons: " << nDeltaElectrons << std::endl;

  // Save histogram (to file and PDF)
  TFile* fout = new TFile(outfile.c_str(), "recreate");

  TNtuple* allMuonStat = new TNtuple("allMuonStat", "All Muon Stats", "N:Eo:x");
  makeTNtupleFromHist(allMuonStat, muEoActiveDist);  

  // get delta_theory by exporting a calculating function from Deltas.cpp.
  // then add delta MC and delta theory 
  // define stack here 

  // Write the hists and plots 

  TEST_HIST->Write();  
  TEST_HIST1->Write(); 
  TEST_HIST2->Write();
 
  //eEnergy->Write();
  eEnergyMuon->Write();
  eEnergyMich->Write();
  eEnergyBrem->Write();
  eEnergyPP->Write();  
  
  eEnergyDelta->Write();
  //eCosDelta->Write();
  eDistMu->Write(); 
  ePocaMu->Write();  
  eTrackDist->Write();  
  //eStack->Write();  

  //mu_x->Write(); 

  muEInit->Write();
  muActiveLen->Write(); 

  muEInitPen->Write();
  muActiveLenPen->Write();  
  muEInitStop->Write();
  muActiveLenStop->Write();  

  allMuonStat->Write();  
  muEoActiveDist->Write();
  muEoActiveDistStop->Write();
  muEoActiveDistPen->Write();  
  muEnergy->Write();  
  muMinusLife->Write();  
  muPlusLife->Write(); 
  muBeta->Write();
  muGamma->Write();
  fout->Close();

  return 0;
}
  
