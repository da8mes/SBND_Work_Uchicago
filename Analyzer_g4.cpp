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
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
//#include "lardata/DetectorInfoServices/DetectorClocksService.h" 
//#include "larcore/Geometry/Geometry.h" 
//#include "larcorealg/Geometry/GeometryCore.h" 
//#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" 
//#include "larsim/Simulation/LArG4Parameters.h"

//===============GLOBAL VARIABLES==================================// 

/*
// Geometry and Clock Providers
geo::GeometryCore const* geoService; // Geometry provider
detinfo::DetectorClocks const* timeService; // Detector clock provider   
const detinfo::DetectorProperties* detprop; // Detector Properties
geoService = lar::providerFrom<geo::Geometry>(); 
timeService = lar::providerFrom<detinfo::DetectorClocksService>(); 
detProp = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
*/ 

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

// populate muEInit and muActiveLen given
// MCTrajectory of the Particle and the corresponding 
// two hists for track length and initial Energy 
void fillEnergyLength(const simb::MCTrajectory& track, TH1F* E, TH1F* len)
{

	simb::MCTrajectory::const_iterator pt; 
	TLorentzVector tPos;
	TLorentzVector tMom;
	TLorentzVector tPosInit; // first active point 
	TLorentzVector tPosFinal; // last active point  
	// advance on the track via step
	int step = 1; // step size  
	// flag used in determining initial Energy
	// of the muon's point of entry 
	int eInitUndetermined = 1;  
	for (pt = track.begin(); pt != track.end(); pt+=step)
	{ 
		tPos = (*pt).first;
		if (not inTPC(tPos)) continue; // wait till in Active Vol  
		tPosFinal = tPos; // get the final active point

		if (eInitUndetermined) // if initial energy undetermined so far 
		   {
			tPosInit = tPos; // save first active point
			tMom = (*pt).second; 
			E->Fill(tMom.E()); // get the initial energy	
			eInitUndetermined = 0; // indicate that init energy is determined
		   } 
	}

	len->Fill( (tPosFinal.Vect() - tPosInit.Vect()).Mag() ); // get the active track length        		
	return; 

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

/* 
// calculates the expected number of delta rays of energy T
// that come out of particle's track 
// Only valid for I << T <= W where I is the Ionization energy of Ar 
// I = 188.000000 ev
// All the enrgy terms have to be in [Mev] 
int nDeltaEvents(const simb::MCTrajectory& traject,  double T)
{
	// MCTrajectory object has type vector<Pair<TLVecPosition, TLVecMomentum>>; 
	// const simb::MCTrajectory& traject = particle.Trajectory();
	// std::cout << "\n--\nMuon Trajectory Info: " << traject << std::endl;
	
	double dN_dT = 0; // what is dN/dT_x0 = 0;   
	const double c2 = 1; // c^2, natural units 
	const double me = 0.511; // Mass of electron [Mev/c^2] 
	const double mu = 105.658; // Mass of Muon  [Mev/c^2]
	const double me_mu = me / mu; 
	const double K = 0.307075; // [MeV mol^-1 cm^2]
	const double z2 = 1; // charge number of muon
	const double Z = 18; // atmonic number of Ar 
	const double A =  39.948; // [g mol^-1] atomic mass of Ar 
	const double Ar_density = 1.3954; // [g / cm^3] density of liquid Ar  
	const double coeff = (0.5*K*z2*Z)/A; // 0.5*k*z^2*Z/A 
	TLorentzVector prevPos; // (0.,0.,0.,0.) in case I need this to find dx
 
	// using the const iterator for traject	
	simb::MCTrajectory::const_iterator step; 
	for (step = traject.begin(); step != traject.end(); step++)
	    {	 
		// Do stuff to trajectory steps	 
		// std::cout << "\n--\nStep Info: " << *step << std::endl; 

		auto const pt = *(step); // extract <Pos, Mom> pair
		auto const pos = pt.first;  
		auto const P = pt.second; 
		double dx = (pos.Vect() - prevPos.Vect()).Mag(); // [cm] get dx
		prevPos = pos; // save point for next step 

		// exludes points outside active volumes of TPCs
		if (not inTPC(pos)) continue;
 
		double beta = P.Beta(); 
		double gamma = P.Gamma(); 
		double eMu = 1000*P.E(); // [Mev] Energy of Muon  
		double W = (2*me*c2*pow(beta*gamma, 2)) / (1+ (2*gamma*me_mu) + pow(me_mu,2)); // max energy transfer 
		if (T > W) return 0; // if delta energy requested is more than max energy transfer 	
		double F = ( ((1/pow(T*beta, 2)) - (1/(T*W)) ) + 0.5*(1/pow(beta*(eMu + mu*c2), 2)) ); // spin dependent term
		// ... divided by beta^2 and T^2
		
		// integrate over dx
		// Multiplied by Ar_density because the 
		// formula in Phys Rev D is events per density
		// Citation: DOI: 10.1103/PhysRevD.98.030001 page 449. 
		dN_dT += Ar_density * coeff * F * dx;  		
	    } 
	return dN_dT; 
} 
*/ 

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
  gStyle->SetOptStat("nemr"); 
  gStyle->SetOptFit(1111);

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


//================= Muon Hists=======================================//

// ================= For Delta Theory ===========//
  // fMuEInit_ are uded for extracting the probability dist of cosmic muon 
  // energies (initial, inside active TPC) from CORSIKA. 
  // fMuLen_ is used for doing the same thing for length of muon Track 
  // in Active TPC.

  TF1* fMuEInit = new TF1("fMuEInit", "[0]*exp(-x/[1])", 0, 100);
  TH1F* muEInit = new TH1F("mu_energy_init", "mu_Energy_Initial", 100, 0, 100); 
  TF1* fMuLen = new TF1("fMuLen", "[0]*exp(-x/[1])", 0, 3000); // might not be exponential
  // max Active TPC Length is 4^2 + 4^2 + 5^2 == 57 [m] == 5700 [cm] 
  TH1F* muActiveLen = new TH1F("mu_len", "mu_ActiveTPC_TLength", 500, 0, 5700); 

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
  //THStack* eStack = new THStack("es", "Electron Energies"); 
 
  // All Electron E
  TH1F* eEnergy = new TH1F("eE", "e_Energy", 1000, 0, 1);
  
  eEnergy->SetFillColor(kBlue);
  // All Muon E
  TH1F* eEnergyMuon = new TH1F("eE_mu", "e_Energy_mu", 500, 0.001, 0.5); 
  eEnergyMuon->SetFillColor(kRed);
  eEnergyMuon->SetXTitle("e Energy [Gev]"); 
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
        // Populate hists for Penetrating and Stopping Muons 
        if (inTPC(muPositionEnd)) // stopping 
           {
		fillEnergyLength(muTrajectory, muEInitStop, muActiveLenStop); 
           }
	else // penetrating 
	   {
		fillEnergyLength(muTrajectory, muEInitPen, muActiveLenPen); 
	   }

	continue;  
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
	     // That Min Energy == 2 [Mev] 
	     // Citation: https://physics.nist.gov/cgi-bin/Star/e_table.pl
	    	     
	     //Apply minimum energy cut (2Mev)
	     // or trackLength == 9 mm == 0.9 cm
	     // anything that travels less can not be reconstructed.
	     //
	     // Note: energy cut and dist length are not identical(?) 
	     //if (dTrackLength < 0.9) continue;
	     if (eE < 0.002) continue; // get rid of drifts 

	     float epsilon = 0.2; // [cm] from eDist plot
	     // this cut throws away tip deltas that might 
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
		    //std::cout << "Mu Minus (?) pdg: " << pdgId << std::endl; 
		    muMinusLife->Fill(time_to_decay);	
		}
		else
		{
		    //std::cout << "Mu Plus (?) pdg: " << pdgId << std::endl; 
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
	     ePocaMu->Fill(pca);

	     // Apply poca cut found from ePocaMU 
	     // pca <= 0.001 cm
	     // Now only delta's should survive
	     if (pca > 0.001)
		{
		 //std::cout << "Start Process for Poca < 0.02: " << dProcess << std::endl;  
		 nBremElectrons++; 
		 eEnergyBrem->Fill(eE); 
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

  // Write the hists and plots 
 
  /* 
  eMomentum->Write(); 
  eMomentumX->Write();
  eMomentumY->Write();
  eMomentumZ->Write();  
  */ 

  //deltaTheory->Write(); 

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

  //muEInit->Write();
  //muActiveLen->Write(); 

  muEInitPen->Write();
  muActiveLenPen->Write();  
  muEInitStop->Write();
  muActiveLenStop->Write();  
 
  muEnergy->Write();  
  muMinusLife->Write();  
  muPlusLife->Write(); 
  muBeta->Write();
  muGamma->Write();
  fout->Close();

  return 0;
}
  
