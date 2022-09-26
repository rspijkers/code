// This is a very basic pythia script. 
// It generates events with settings specified in a separate config file.
// It calculates some properties on the fly, and creates some histograms
// which it writes to an output file. 
//
// First one should enter the required pythia and ROOT environments using
// alienv. Then one can make an executable via `make pythia_example`, and
// execute it with `./pythia_example "outputfilename.root"`. 

// Import the required libraries:
// std
#include <iostream> // std in and out
#include <cmath> // needed for modulo in deltaPhi calc
#include <vector> // needed for a variable array
#include <cstring> // string/char handling
#include <chrono> // for measuring performance of script
// ROOT and pythia
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TTree.h"

#define PI 3.14159265

using namespace Pythia8;
using std::cout; using std::endl;

bool IsStrange(Int_t particlepdg) {
	// Checks if given pdg code is that of a (anti)strange hadron, returns True/False
	// This also finds ssbar mesons such as the phi
	Int_t pdg = std::abs(particlepdg);
	pdg /= 10; // get rid of the last digit, is not important for quark content
	if (pdg % 10 == 3) return true; // 3rd quark
	pdg /= 10;
	if (pdg % 10 == 3) return true; // 2nd quark
	pdg /= 10;
	if (pdg % 10 == 3) return true; // 1st quark
	return false;
}

Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	// returns phi1 - phi2 in a range between -pi/2 and 3pi/2
	return std::fmod(phi1 - phi2 + 2.5*PI, 2*PI) - 0.5*PI;
}

int main(int argc, char** argv)
{
	if (argc != 2) { // main() counts as one argument, together with outputfile path makes 2
	cerr << " Unexpected number of command-line arguments. \n"
         << " You are expected to provide a file path for the output and nothing else. \n"
         << " Program stopped! " << endl;
    return 1;
	}

	// start keeping track of time
    auto start = std::chrono::high_resolution_clock::now();

	// create outputfile
	TFile* outFile = new TFile(argv[1], "CREATE"); // doesn't open file if it already exists 
	if(!outFile->IsOpen()) { // if output file isn't opened, abort program
		cerr << "Error: File " << argv[1] << " is not opened, perhaps because it already exists. Aborting script.";
		return 1;
	}

	const Double_t pTmin = 0.15; // minimum pT
	const Double_t etamax = 4.; // maximum value of |eta|

	Pythia pythia;

	// PYTHIA SETTINGS
	// put all settings in a separate file
	pythia.readFile("example_settings.cmnd");
	Int_t nEvents = pythia.mode("Main:numberOfEvents");

	// this is a rather ugly way of generating a random seed
	// pretty much only important if you run multiple instances of this script at the same time. 
	Int_t processid = getpid();
	string seedstr = "Random:seed = " + std::to_string((time(0) + processid)%900000000);
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedstr); // 0 means it uses the time to generate a seed

	pythia.init();

	// create output histo's
  	TH2D *hEtaPt = new TH2D("hEtaPt","p_{T} vs #eta for all particles;#eta;p_{T} (GeV/c)", 40, -4., 4., 50, 0, 10);
	TH1D *hPt = new TH1D("hPt", "Transverse momentum of all particles;p_{T} (GeV/c)", 50, 0, 10); 
	TH1D *hEta = new TH1D("hEta", "Pseudorapidity of all particles;#eta", 40, -4, 4); 
	// It is customary to report dphi between -pi/2 and 3pi/2, this way it is intuitive to see peaks at dphi = 0, pi
	TH1D *hXiDPhi = new TH1D("hXiDPhi", "Delta #varphi between a Xi-hadron other strange hadrons in the same event", 100, -PI/2, 3*PI/2); 
	// use Double_t (TH1D) to get around maximum bin content of Int_t:
	TH1D *hPDG = new TH1D("hPDG", "PDG codes of all particles", 12000, -6000, 6000); 
	TH1D *hPDGStrange = new TH1D("hPDGStrange", "PDG codes of strange hadrons", 12000, -6000, 6000); 

	cout << "Generating " << nEvents << " events..." << endl;

	// event loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++) {
		if(!pythia.next()) continue; // generate next event: skip the event if pythia.next() generates "false"

		int nPart = pythia.event.size(); // number of particles in this event
		
		// particle loop
		for(int iPart = 0; iPart < nPart; iPart++) {
      		const Particle &particle = pythia.event[iPart];
			if(!particle.isFinal()) continue; // skip if not a final state particle 

			Double_t pT = particle.pT();
			Double_t eta = particle.eta();

			if(pT < pTmin || eta > etamax) continue; // if particle fails kinematic checks skip it
			Int_t pdg = particle.id();

			hEtaPt->Fill(eta, pT); // fill 2D histo
			hPt->Fill(pT);
			hEta->Fill(eta);
			hPDG->Fill((Double_t) pdg); // cast pdg to double before filling
			if(IsStrange(pdg)) hPDGStrange->Fill((Double_t) pdg);

			// Maybe you are interested in the deltaphi between any Xi^{+/-} hadron you find and all other strange hadrons in the event
			// Then you could do something like this
			if(abs(pdg) == 3312){ // 3312 is the code for a Xi^{-}, -3312 is for Xi^{+}
				// we found a Xi^{+/-}! let's loop over all the other particles in the event
				Double_t phi = particle.phi();
				for(int jPart = 0; jPart < nPart; jPart++) {
					if(jPart == iPart) continue; // don't correlate particle with itself
					const Particle &part2 = pythia.event[jPart];
					Double_t part2pT = part2.pT();
					Double_t part2eta = part2.eta();
					Int_t part2pdg = part2.id();
					// See if it is strange, final state, and passes kinematic checks
					if(IsStrange(part2pdg) && part2.isFinal() && part2pT > pTmin && part2eta < etamax){
						Double_t dPhi = DeltaPhi(phi, part2.phi());
						hXiDPhi->Fill(dPhi);
					} // end check if 2nd particle passes checks
				} // end 2nd particle loop
			} // end if-statement checking for Xi
		} // end 1st particle loop
	} // end event loop

	// write output to file and close it
	outFile->Write();
	cout << "Output written to file " << outFile->GetName() << endl;
	outFile->Close();

	// stop keeping track of time, and calculate how long it took
    auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(end - start);
    cout << "This script took " << duration.count() << " minutes to run." << endl;
} // end main
