// std
#include <iostream>
#include <cmath> // needed for modulo in deltaPhi calc & abs of doubles
#include <vector> // needed for a variable array
#include <cstring> // string/char handling
#include <chrono> // for measuring performance of script
// ROOT/pythia
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
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

// TODO: move this to helperfunctions.cpp?
Double_t DeltaPhi(Double_t phi1, Double_t phi2){
	// returns phi1 - phi2 in a range between -pi/2 and 3pi/2
	return std::fmod(phi1 - phi2 + 2.5*PI, 2*PI) - 0.5*PI;
}

int main(int argc, char** argv)
{	
	Int_t runNr = 0;
	// main() counts as one argument
	if (argc < 2) {
		cerr << "Error: Too few command-line arguments! You are expected to give at least a file path for the output." << endl;
		return 1;
	} else if (argc == 3) {
		try{
			runNr = std::stoi(argv[2]);
		} catch(std::invalid_argument) {
			cerr << "Error: Could not convert the run number to an integer! " << endl;
			return 1;
		}
	} else if (argc > 3) {
		cerr << "Error: Too many command-line arguments! You are expected to give at most a file path for the output and a run number. " << endl;
		return 1;
	}
	// TODO: check if argument is a valid path?		

	const char* outFilePath = argv[1];

	// start keeping track of time
    auto start = std::chrono::high_resolution_clock::now();

	// TODO: get kinematic options from json config file?

	// Analysis settings
	const Bool_t	doUnderlyingEvent = false;
	const Double_t 	pTmin = 0.15;
	const Double_t 	maxEta = 4.;

	Pythia pythia;

	// PYTHIA SETTINGS
	pythia.readFile("pythia_settings/ssbar_monash.cmnd");
	Int_t nEvents = pythia.mode("Main:numberOfEvents");
	const Int_t eventIdOffset = nEvents*runNr;

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

	// TODO structure histo's into QA, UE, and correlations?
	// maybe do this by initializing the directories before this, and cd() to the right directory in between the histo declarations?

	// Output histo's
  	TH2D *hEtaPt = new TH2D("hEtaPt","p_{T} vs #eta for all particles;#eta;p_{T} (GeV/c)", 40, -4., 4., 50, 0, 10);
	TH1D *hPDG = new TH1D("hPDG", "PDG code for trigger strange hadrons", 12000, -6000, 6000); 
	TH1D *hPDGAssoc = new TH1D("hPDGAssoc", "PDG code for associated strange hadrons", 12000, -6000, 6000); // use Double_t to get around maximum bin content of Int_t
	TH1I *hCandidatesPerEvent = new TH1I("hCandidatesPerEvent", "Number of candidates per event", 100, 0, 100);

	Int_t partpdg;
	Double_t partpT;
	Double_t parteta;
	Double_t pTssbar;
	
	SmallEvent *event = new SmallEvent();
	SmallTrack track;
	TTree *tree = new TTree("tree", "tree with small events and small tracks");
	tree->Branch("event", &event); 

	cout << "Generating " << nEvents << " events..." << endl;
	for(int iEvent = 0; iEvent < nEvents; iEvent++) { // event loop
		if(!pythia.next()) continue;

		int nPart = pythia.event.size();
		int candidatesPerEvent = 0;
		int nFinalState = 0;
		pTssbar = -1.;
		
		for(int iPart = 0; iPart < nPart; iPart++) { // particle loop
      		const Particle &part = pythia.event[iPart];

			partpdg = part.id();
			// in case of ssbar, save the highest pT. in case only one s(bar), pT always > -1 
			if(part.status()==-23 && abs(partpdg)==3 && part.pT() > pTssbar) pTssbar = part.pT();
			
			if(!part.isFinal()) continue; // final state particle 
			nFinalState++;
			partpT = part.pT();
			parteta = part.eta();

			if(!IsStrange(partpdg) || partpT < pTmin || abs(parteta) > maxEta) continue; // kine cuts & strangeness check
			// we have identified a strange trigger that satisfies the kinematic requirements
			// If we get this far with the trigger particle, we will correlate it with other strange hadrons
			// In order to be able to normalize, we need to keep track of how many triggers we have for each hadron

			Double_t partphi = part.phi();

			track.setPDG(partpdg);
			track.setpT(partpT);
			track.setEta(parteta);
			track.setPhi(partphi);
			event->addCandidate(track);
			track.Clear();

			candidatesPerEvent++;

			hEtaPt->Fill(parteta,partpT);
			hPDG->Fill((Double_t) partpdg);
		} // end track loop
		event->setEventId(eventIdOffset+iEvent);
		event->setNtracks(nFinalState);
		tree->Fill();
		event->Clear();
		hCandidatesPerEvent->Fill(candidatesPerEvent);
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
