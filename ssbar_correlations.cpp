#include <iostream>
#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TTree.h"
#include "THnSparse.h" // obsolete?
#include <cmath> // needed for modulo in deltaPhi calc
#include <vector> // needed for a variable array
#include <cstring> // string/char handling

#define nEvents 5000	// 20000 // TODO: get this from cmnd settings file?
#define PI 3.14159265

using namespace Pythia8;

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
	if (argc != 2) { // main() counts as one argument, together with outputfile path makes 2
	cerr << " Unexpected number of command-line arguments. \n"
         << " You are expected to provide a file path for the output and nothing else. \n"
         << " Program stopped! " << endl;
    return 1;
	}

	// TODO: check if argument is a valid path?

	int mecorr=1;
	const Double_t pTminTrigg = 4.;
	const Double_t pTminAssoc = 0.;
	const Double_t maxEta = 4.;

	Pythia pythia;

	// PYTHIA SETTINGS
	pythia.readFile("ssbar_correlations.cmnd");

	Int_t processid = getpid();
	string seedstr = "Random:seed = " + std::to_string((time(0) + processid)%900000000);
	pythia.readString("Random:setSeed = on");
	pythia.readString(seedstr); // 0 means it uses the time to generate a seed

	//ME corrections
	//use of matrix corrections where available
	if(mecorr==0){
		pythia.readString("TimeShower:MECorrections=off");
	}
	pythia.init();

	TFile* outFile = new TFile(argv[1], "CREATE"); // doesn't open file if it already exists 
	if(!outFile->IsOpen()) { // if output file isn't opened, abort program
		cerr << "Error: File " << argv[1] << "is not opened, perhaps because it already exists. Aborting script.";
		return 1;
	}

	// Output histo's
  	TH2D *hEtaPt = new TH2D("hEtaPt","p_{T} vs #eta for all particles;#eta;p_{T} (GeV/c)", 40, -maxEta, maxEta, 50, 0, 10);
	TH1D *hPDG = new TH1D("hPDG", "PDG code for trigger strange hadrons", 8000, -4000, 4000); // use Double_t to get around maximum bin content of Int_t

	// TODO: create tree to save all the data of trigger/assoc pairs to
	// not sure how to handle multiple assocs per trigger
	Int_t partpdg;
	Double_t partpT;
	Double_t parteta;
	std::vector<Int_t> pdgAssoc;
	std::vector<Double_t> pTAssoc;
	std::vector<Double_t> etaAssoc;
	std::vector<Double_t> deltaPhi;
	std::vector<Double_t> deltaEta;
	TTree *tree = new TTree("tree", "tree with trigger/assoc strange hadrons");
	tree->Branch("pdgTrigger", &partpdg, "partpdg/I");
	tree->Branch("pTTrigger", &partpT, "partpT/D");
	tree->Branch("etaTrigger", &parteta, "parteta/D");
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
		
		for(int iPart = 0; iPart < nPart; iPart++) {
      		const Particle &part = pythia.event[iPart];
			partpT = part.pT();
			parteta = part.eta();
			hEtaPt->Fill(parteta,partpT);
			if(partpT > pTminTrigg && part.isFinal() && std::abs(parteta) < maxEta) { // final state particle + kine cuts
				// we have identified a potential trigger that satisfies the kinematic requirements
				partpdg = part.id();
				if(IsStrange(partpdg)) {
					// If we get this far with the trigger particle, we will correlate it with other strange hadrons
					// In order to be able to normalize, we need to keep track of how many triggers we have for each hadron, because in the next part we will fill the trigger particle info for each pair, which could be more than one.
					hPDG->Fill((Double_t) partpdg);

					// Clear the vectors with the associated/correlation variables
					pdgAssoc.clear(); pTAssoc.clear(); etaAssoc.clear(); deltaPhi.clear(); deltaEta.clear();
					for(int jPart = 0; jPart < nPart; jPart++) {
						const Particle &part2 = pythia.event[jPart];
						Double_t part2pT = part2.pT();
						Double_t part2eta = part2.eta();
						Int_t part2pdg = part2.id();
						if(IsStrange(part2pdg) && part2.isFinal() && partpT > part2pT && part2pT > pTminAssoc && std::abs(part2eta) < maxEta){
							Double_t dPhi = std::fmod(part.phi() - part2.phi() + 2.5*PI, 2*PI) - 0.5*PI; 
							Double_t dEta = parteta - part2eta;
							pdgAssoc.push_back(part2pdg);
							pTAssoc.push_back(part2pT);
							etaAssoc.push_back(part2eta);
							deltaPhi.push_back(dPhi);
							deltaEta.push_back(dEta);
						}
					} // end assoc loop
					tree->Fill();
				}
			}
		} // end trigger loop
	} // end event loop

	outFile->Write();
	cout << "Tree written to file " << outFile->GetName() << endl;
	outFile->Close();
} // end main
