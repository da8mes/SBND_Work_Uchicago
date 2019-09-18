/* Deltas.cpp
 *
 *  * Written by: D. Belayneh <dawit@uchicago.edu>, 2019/12/10
 *
 *
 */

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>

// ROOT classes
#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TGraph.h"
#include "TMath.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/GaussIntegrator.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "TStyle.h"

//===============GLOBAL VARIABLES==================================//
const double ERRORLIMIT = 1E-3;
// x =  [cm] and T = [Mev]
const double Ar_density = 1.3973;  
const double D = (Ar_density*18*0.307075)/(2*39.948);
const double Mmu = 105.6583745; // [Mev / c2]
const double Mel = 0.511; // [Mev / c2]
const double Mmu_Mel = Mmu / Mel; // for convinence 
const int c = 1; // unit_less
const double Emu = Mmu * pow(c,2); // [Mev]    
const double Eel = Mel * pow(c,2); // [Mev]
const double I = 0.000188; // [Mev] Mean Ionization Energy of Argon
std::map<double, double> Energy; // map of <double x, double E(x)>; for speed.   
// ----- //  
const double x_step = 0.01; // x step size for integration and Energy loss functions
const double TMIN = 0;
const double TMAX = 1; // [Gev]
const double BINWIDTH = 0.01;
const double NBINS = (TMAX - TMIN)/BINWIDTH; // [Mev] for now 
// deltas->SetXTitle("E [Mev]");  
//===============END OF GLOBAL VARIABLES===========================// 

// relativistic gamma
// gamma(x) via E_x
double gamma(const double E_x)
{
	double gamma_val = E_x / Emu; 
	return gamma_val; 
}

// relativistic beta2
// beta2(x) via E_x
double beta2(const double E_x)
{
	double beta_val = 1 - pow(1/gamma(E_x), 2); 
	return beta_val; 
}

// Calculates differential energy loss
// of a muon with energy E if it was to travel 
// a distance of x in Liquid Ar
// formula: Landau-Vavilov-Bischel most probable energy loss
// Ref: DOI: 10.1103/PhysRevD.98.030001, 447. 
double dE(const double E_x, const double x)
{
	if (x < x_step) return 0; 
	// omiting density correction terms for now 
	// off with accurate prediction by about 10% 
	double delta_E = pow(E_x,2) - pow(Emu,2); 
	double log_dep = TMath::Log((2*delta_E)/(Emu*I)) + TMath::Log(D*x/(beta2(E_x)*I)) + 0.200 - beta2(E_x);
	double dE_val = D*(1/beta2(E_x))*x*log_dep;  
	return dE_val; 
} 

double dETest(const double E_x, const double x)
{
	return 2.0; 
}

// Calculates E(x) = Eo - dE(x) 
// for a given initial energy Eo
// This should be interated to find the energy dissipation 
// of a muon over its track 
double E(const double xf) //const double Eo, const double xf) // going to be slow for the integral
{
	double E;
	double x;  
	if (!Energy.empty()) // current E(x)
	{
		x = (--Energy.end())->first; 
		E = (--Energy.end())->second;	
	}

	if (abs(x - xf) <= x_step) return E; 

	for (; x < xf ; x+=x_step)
	{
		// beta2(E) != 0 for dE to make sense
		if (beta2(E) < 0.001) return E; 
		double loss = dE(E, x); 
		E-=dE(E, x);
		Energy[x] = E; // update E(x) map  	
	}
	return E; 
}

// Calculates the Maximum Energy transfer to deltas
// W(x) via E_x 
double W(const double E_x)
{
	double delta_E = pow(E_x, 2) - pow(Emu,2); 
	double w_val = (2*delta_E) / (Emu * Mmu_Mel + 2*E_x + Eel);  
	return w_val; 
}

// Spin dependent term (for Muon Spin = 1/2)
// F(x, T) via E_x 
double F(const double E_x, const double T)
{
	double f_val = (1 - ((beta2(E_x)*T)/ W(E_x))) + 0.5*pow(T/(E_x + Mmu), 2); 	
	return f_val; 
}


// diff of distribution number N with respect to 
// x (distance travelled) and T (delta energy [Mev])
// Input R: {x, T} is an array of doubles 
// Output dNdTdx at R given Eo. 
double dN_dTdx(const double Eo, const double* R)
{
	double x = R[0]; 
	double T = R[1];

	double Tmax = W(Eo); // to protect from bad T values
	// protected from bad x values in double E();  

	// this should be inside the applyHist function. Temp; 
	if (T > Tmax) return 0; // b/c in this case prob == 0.
 
	// compute E(x) once and propagate it through the other funcitons
	// F(E_x, x, T); W(E_x); beta2(E_x); gamma(E_x);
	double E_x = E(x); 
	double diff_val = (D*F(E_x, T)) / (pow(T,2)*beta2(E_x)); 	
	return diff_val; 
}

// returns a lambda function that curries dN_dTdx with a given 
// Eo 
auto curried_dN(double Eo)
{
	return [=](const double* R){return dN_dTdx(Eo, R);}; 
}

// Integrator
// a[2] = {xi, xf} 
// b[2] = {Ti, Tf} 
double integrate(const double Eo, const double* a, const double* b)
{
	auto fn = curried_dN(Eo);
	// Wraps the function as a functor (whatever that is)
	// Define the rectangular region to integrate over
		
	ROOT::Math::Functor wf(fn, 2); 
	ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);  
	ig.SetFunction(wf); 
	double val = ig.Integral(a, b);

	/*
	std::cout << "Eo: " << Eo << std::endl;	
	std::cout << "Wmax: " << W(Eo) << std::endl;	
	std::cout << "xi: " << a[0] << std::endl;    	
	std::cout << "xf: " << b[0] << std::endl;    	
	std::cout << "Ti: " << a[1] << std::endl;    	
	std::cout << "Tf: " << b[1] << std::endl;    	
	std::cout << "# of Deltas expected: " << val << std::endl; 
	std::cout << "\n------\n" << std::endl;
 	*/ 
	return val; 
}

// applies integrator to bins THist
// checks against bins with bad values like T > Tmax = W(Eo) 
// are set to 0 probabbility
// can check against bad x values as well
// look at a given bin T1; look at T1min and T1max; if within
// permissible T range, find prob else set to 0. 
void makeHist(double Eo, double activeDist, double N, TH1F* hist)
{
	// Eo in Mev; hist in Gev 
	// Tmax is maximum energy transfer
	double Tmax = 1000 * W(Eo);  
	
	int nBin = NBINS;
	for (int i = 1; i <= nBin; i++)
	{
		// double T = hist->GetBinCenter(i);
		double Ti = 1000 * (hist->GetXaxis()->GetBinLowEdge(i));
		double Tf = 1000 * (hist->GetXaxis()->GetBinUpEdge(i)); 
		double val; 
		if (Tf > Tmax || Ti < 2) // outside valid range
		{
			val = 0;
			hist->SetBinContent(i, val); 
			continue; 	
		}

		double a[2] = {0, Ti};  
		double b[2] = {activeDist, Tf}; 
		val = N*integrate(Eo, a, b);
		hist->SetBinContent(i, val);  	 
	}
}

// The global delta distribution is the sum of the individual muon 
// delta distributions
int main(int argc, char* argv[]) 
{ 	
 
	TH1F* deltas = new TH1F("deltas", "Delta Distribution", NBINS, TMIN, TMAX);
	deltas->SetXTitle("E [Gev]");
	deltas->SetFillColor(kNeon);

	// open TNtuple; get N,Eo,activeDist
	const char* infile = argv[1];
	TFile* fileIn = TFile::Open(infile, "UPDATE");

	TTreeReader reader("allMuonStat", fileIn);
	TTreeReaderValue<Float_t> N(reader, "N");    
	TTreeReaderValue<Float_t> Eo(reader, "Eo"); 
	TTreeReaderValue<Float_t> activeDist(reader, "x"); 

	// Loop over all entries in TTree
	while (reader.Next()) 
	{		
		if (*N < 1) continue; 

		std::cout << "Eo [Gev]: " << *Eo << std::endl;    	
		std::cout << "x [cm]: " << *activeDist << std::endl;    	
		std::cout << "N: " << *N << std::endl; 
		std::cout << "\n------\n" << std::endl;

		Energy[0] = 1000*(double)*Eo;
		makeHist(Energy[0], (double)*activeDist, (double)*N, deltas);
	} 
	// normalize to inegral == 1
	Double_t scale = 1/deltas->Integral(); 
	deltas->Scale(scale); 
	deltas->DrawCopy(); 

	deltas->Write(); 
	fileIn->Close(); 

	return 0;  
} 
