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
	if (argc != 2) { // main() counts as one argument, together with outputfile path makes 2
	cerr << " Unexpected number of command-line arguments. \n"
         << " You are expected to provide a file path for the output and nothing else. \n"
         << " Program stopped! " << endl;
    return 1;
	}

	// start keeping track of time
    auto start = std::chrono::high_resolution_clock::now();

	// TODO: check if argument is a valid path?

	// TODO: get kinematic options from json config file?

	const Double_t pTminTrigg = 0.8;
	const Double_t pTminAssoc = 0.15;

	const Double_t maxEta = 4.;

	Pythia pythia;

	// PYTHIA SETTINGS
	pythia.readFile("pythia_settings/ssbar_monash.cmnd");
	Int_t nEvents = pythia.mode("Main:numberOfEvents");

	Int_t processid = getpid();
	string seedstr = "Random:seed = " + std::to_string((time(0) + processid)%900000000);
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedstr); // 0 means it uses the time to generate a seed

	pythia.init();

	// OUTPUT INIT
	TFile* outFile = new TFile(argv[1], "CREATE"); // doesn't open file if it already exists 
	if(!outFile->IsOpen()) { // if output file isn't opened, abort program
		cerr << "Error: File " << argv[1] << " is not opened, perhaps because it already exists. Aborting script.";
		return 1;
	}

	// Output histo's
  	TH2D *hEtaPt = new TH2D("hEtaPt","p_{T} vs #eta for all particles;#eta;p_{T} (GeV/c)", 40, -4., 4., 50, 0, 10);
	TH1D *hPDG = new TH1D("hPDG", "PDG code for trigger strange hadrons", 12000, -6000, 6000); 
	TH1D *hPDGAssoc = new TH1D("hPDGAssoc", "PDG code for associated strange hadrons", 12000, -6000, 6000); // use Double_t to get around maximum bin content of Int_t
	TH1I *hTriggersPerEvent = new TH1I("hTriggersPerEvent", "number of triggers per event", 20, 0, 20);

	Int_t partpdg;
	Double_t partpT;
	Double_t parteta;
	Double_t pTssbar;
	std::vector<Int_t> pdgAssoc;
	std::vector<Double_t> pTAssoc;
	std::vector<Double_t> etaAssoc;
	std::vector<Double_t> deltaPhi;
	std::vector<Double_t> deltaEta;
	TTree *tree = new TTree("tree", "tree with trigger/assoc strange hadrons");
	tree->Branch("pdgTrigger", &partpdg, "pdgTrigger/I");
	tree->Branch("pTTrigger", &partpT, "pTTrigger/D");
	tree->Branch("etaTrigger", &parteta, "etaTrigger/D");
	tree->Branch("pTssbar", &pTssbar, "pTssbar/D");
	tree->Branch("pdgAssoc", &pdgAssoc);
	tree->Branch("pTAssoc", &pTAssoc);
	tree->Branch("etaAssoc", &etaAssoc);
	tree->Branch("deltaPhi", &deltaPhi);
	tree->Branch("deltaEta", &deltaEta);

	cout << "Generating " << nEvents << " events..." << endl;
	// event loop
	for(int iEvent = 0; iEvent < nEvents; iEvent++) {
		if(!pythia.next()) continue;

		int nPart = pythia.event.size();
		int strangenessPerEvent = 0;
		int triggersPerEvent = 0;
		pTssbar = -1.;
		
		for(int iPart = 0; iPart < nPart; iPart++) {
      		const Particle &part = pythia.event[iPart];

			partpdg = part.id();
			// in case of ssbar, save the highest pT. in case only one s(bar), pT always > -1 
			if(part.status()==-23 && abs(partpdg)==3 && part.pT() > pTssbar) pTssbar = part.pT();
			
			if(!part.isFinal()) continue; // final state particle 
			partpT = part.pT();
			parteta = part.eta();

			if(!IsStrange(partpdg) || partpT < pTminTrigg || abs(parteta) > maxEta) continue; // kine cuts & strangeness check
			// we have identified a strange trigger that satisfies the kinematic requirements
			// If we get this far with the trigger particle, we will correlate it with other strange hadrons
			// In order to be able to normalize, we need to keep track of how many triggers we have for each hadron
			Double_t partphi = part.phi();
			triggersPerEvent++;

			// Clear the vectors with the associated/correlation variables
			pdgAssoc.clear(); pTAssoc.clear(); etaAssoc.clear(); deltaPhi.clear(); deltaEta.clear();
			for(int jPart = 0; jPart < nPart; jPart++) {
				if(jPart == iPart) continue; // don't correlate particle with itself
				const Particle &part2 = pythia.event[jPart];
				Double_t part2pT = part2.pT();
				Double_t part2eta = part2.eta();
				Double_t part2phi = part2.phi();
				Int_t part2pdg = part2.id();

				if(!IsStrange(part2pdg) || !part2.isFinal() || part2pT < pTminAssoc || abs(part2eta) > maxEta) continue; // all cuts at once

				Double_t dPhi = DeltaPhi(partphi, part2phi);
				Double_t dEta = parteta - part2eta;
				pdgAssoc.push_back(part2pdg);
				pTAssoc.push_back(part2pT);
				etaAssoc.push_back(part2eta);
				deltaPhi.push_back(dPhi);
				deltaEta.push_back(dEta);
			} // end assoc loop

			hEtaPt->Fill(parteta,partpT);
			hPDG->Fill((Double_t) partpdg);
			tree->Fill();

			// only fill the associated pdg's histo once per event
			if(triggersPerEvent == 1) for(int x : pdgAssoc) hPDGAssoc->Fill(x);
		} // end trigger loop
		hTriggersPerEvent->Fill(triggersPerEvent);
	} // end event loop
	cout << "test" << endl;
	outFile->Write();
	cout << "Output written to file " << outFile->GetName() << endl;
	outFile->Close();

	// stop keeping track of time, and calculate how long it took
    auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(end - start);
    cout << "This script took " << duration.count() << " minutes to run." << endl;
} // end main
