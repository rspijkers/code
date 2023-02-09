// This is a PYTHIA script that generates pp events and saves them to a root file. 
// It is taylored to save strange hadrons so we can analyze correlations between them after generating.
// This script should be shipped with a makefile, custom event and track classes, and PYTHIA config files.
//
// Usage:
// After making the executable, it can be run by `./ssbar_correlations "outputfile.root" runNumber`
// `runNumber` is an optional parameter that allows for easier parallel generating, This way you can have
// uniquely identified events across multiple instances of this script. Handy if you want to run on a grid
// or something similar. It is assumed that runNumber starts from 1, entering runNumber=0 will lead to 
// negative eventId's.
//
// Author: Rik Spijkers (rik.spijkers@nikhef.nl)

// std
#include <iostream>
#include <cmath> // needed for modulo in deltaPhi calc & abs of doubles
#include <vector> // needed for a variable array
#include <cstring> // string/char handling
#include <chrono> // for measuring performance of script
// ROOT/pythia
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH3.h"
#include "TTree.h"
// custom
#include "helperfunctions.h"
#include "SmallEvent.h"

#define PI 3.14159265

using namespace Pythia8;
using std::cout; using std::endl;

bool IsStrange(Int_t particlepdg) {
	// Checks if given pdg code is that of a (anti)strange hadron, returns True/False
	// This also find ssbar mesons such as the phi
	Int_t pdg = std::abs(particlepdg);
	pdg /= 10; // get rid of the last digit, is not important for quark content
	if (pdg % 10 == 3) return true; // 3rd quark
	pdg /= 10;
	if (pdg % 10 == 3) return true; // 2nd quark
	pdg /= 10;
	if (pdg % 10 == 3) return true; // 1st quark
	return false;
}

int main(int argc, char** argv)
{	
	// we expect the following arguments:
	// outputfilepath, cmnd filepath, runNr (optional)
	Int_t runNr = 0;
	// main() counts as one argument
	if (argc < 2) {
		cerr << "Error: Too few command-line arguments! You are expected to give at least a file path for the output." << endl;
		return 1;
	} else if (argc == 4) {
		try{
			runNr = std::stoi(argv[3]);
		} catch(std::invalid_argument) {
			cerr << "Error: Could not convert the run number to an integer! " << endl;
			return 1;
		}
	} else if (argc > 4) {
		cerr << "Error: Too many command-line arguments! You are expected to give at most a file path for the output and a run number. " << endl;
		return 1;
	}
	// TODO: check if argument is a valid path?		

	const char* outFilePath = argv[1];
	const char* pythiaOptions = argv[2];

	// start keeping track of time
    auto start = std::chrono::high_resolution_clock::now();

	// TODO: get kinematic options from json config file?

	// Analysis settings
	const Bool_t	doOmegaTrigger = true; 
	const Double_t 	pTmin = 0;
	const Double_t 	maxEta = 4.;

	Pythia pythia;

	// PYTHIA SETTINGS
	pythia.readFile(pythiaOptions);
	// pythia.readFile("pythia_settings/ssbar_monash.cmnd");
	Int_t nEvents = pythia.mode("Main:numberOfEvents");
	const Int_t eventIdOffset = nEvents*(runNr-1);

	Int_t processid = getpid();
	string seedstr = "Random:seed = " + std::to_string((time(0) + processid)%900000000);
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedstr); // 0 means it uses the time to generate a seed

	pythia.init();

	// OUTPUT INIT
	TFile* outFile = new TFile(outFilePath, "CREATE"); // doesn't open file if it already exists 
	if(!outFile->IsOpen()) { // if output file isn't opened, abort program
		cerr << "Error: File " << outFilePath << " is not opened, perhaps because it already exists. Aborting script.";
		return 1;
	}

	// Output histo's
	// Use THnSparse with PDG, pT, eta, since most bins are empty due to the PDG codes
	TH3D* hSpectra = new TH3D("hSpectra","PDG, p_{T}, #eta for all non-exotic particles (|PDG| < 4000);PDG;p_{T} (GeV/c);#eta", 8000, -4000, 4000, 100, 0, 20, 100, -maxEta, maxEta);

	Int_t pdg;
	Double_t pT;
	Double_t eta;
	Double_t phi;
	Double_t pTssbar;
	
	SmallEvent *event = new SmallEvent();
	SmallTrack track;
	TTree *tree = new TTree("tree", "tree with small events and small tracks");
	tree->Branch("event", &event); 

	cout << "Generating " << nEvents << " events..." << endl;
	for(int iEvent = 0; iEvent < nEvents; iEvent++) { // event loop
		if(!pythia.next()) continue;

		int nPart = pythia.event.size();

		if(doOmegaTrigger){
			bool trigger = false;
			for (int iPart = 0; iPart < nPart; iPart++){
				int abspdg = abs(pythia.event[iPart].id());
				if (abspdg == 3334){ // Omega pdg = 3334
					trigger = true;
					break;
				} 
			}
			if(!trigger) continue; // no trigger? go next!
			// If we do have a trigger, proceed as normal. 
		}

		int candidatesPerEvent = 0;
		int nFinalState = 0;
		pTssbar = -1.;
		
		for(int iPart = 0; iPart < nPart; iPart++) { // particle loop
      		const Particle &part = pythia.event[iPart];

			pdg = part.id();
			// in case of ssbar, save the highest pT. in case only one s(bar), pT always > -1 
			if(part.status()==-23 && abs(pdg)==3 && part.pT() > pTssbar) pTssbar = part.pT();
			
			if(!part.isFinal()) continue; // final state particle 
			nFinalState++;
			pT = part.pT();
			eta = part.eta();

			// Fill the UE histo
			hSpectra->Fill((Double_t) pdg, pT, eta);

			if(!IsStrange(pdg) || pT < pTmin || abs(eta) > maxEta) continue; // kine cuts & strangeness check

			phi = part.phi();

			track.setPDG(pdg);
			track.setpT(pT);
			track.setEta(eta);
			track.setPhi(phi);
			event->addCandidate(track);
			track.Clear();

			candidatesPerEvent++;
		} // end track loop
		event->setEventId(eventIdOffset+iEvent);
		event->setNtracks(nFinalState);
		event->setpTssbar(pTssbar);
		tree->Fill();
		event->Clear();
	} // end event loop
	outFile->Write();
	cout << "Output written to file " << outFile->GetName() << endl;
	outFile->Close();

	// stop keeping track of time, and calculate how long it took
    auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "This script took " << duration.count()/60. << " minutes to run." << endl;
	return 0;
} // end main
