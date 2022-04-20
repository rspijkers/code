#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TPad.h"
#include <vector>
#include <unordered_map>
#include "hadrons.h"

void tree_handler_v2() {
    
    TFile *inputfile = new TFile("output/ssbar_5M_14TeV_skands_mode2_noEta.root", "READ");
    TTree *tree = (TTree*) inputfile->Get("tree");
    // tree->Print();
    // tree->Show(0);

    TFile *outputfile = new TFile("output/plots.root", "RECREATE");

    // Set branch address, so we can use the variables when we GetEntry()
    Int_t pdgTrigger;
    std::vector<int>* pdgAssoc = 0;
    Double_t pTTrigger, etaTrigger;
    std::vector<double>* pTAssoc = 0;
    std::vector<double>* etaAssoc = 0;
    std::vector<double>* deltaPhi = 0;
    std::vector<double>* deltaEta = 0;
    tree->SetBranchAddress("pdgTrigger", &pdgTrigger);
    tree->SetBranchAddress("pdgAssoc", &pdgAssoc);
    tree->SetBranchAddress("pTTrigger", &pTTrigger);
    tree->SetBranchAddress("etaTrigger", &etaTrigger);
    tree->SetBranchAddress("pTAssoc", &pTAssoc);
    tree->SetBranchAddress("etaAssoc", &etaAssoc);
    tree->SetBranchAddress("deltaPhi", &deltaPhi);
    tree->SetBranchAddress("deltaEta", &deltaEta);

    std::vector<Hadron*> hadron_vec = {Kminus, Lambda, Sigmaminus, Sigmazero, Sigmaplus, Ximinus, Xizero, Omegaminus};
    // std::vector<Hadron*> hadron_vec = {Kminus, Lambda}; // less particles for testing purposes
    Int_t ndim = hadron_vec.size(); 
    
    struct mapStruct
    {
        TString name;
        TH1D* histogram; // same for strange and antistrange
        Int_t strangeness;

        Int_t ntriggers = 0;
        Double_t chargedKaonBkg = 0;
        Double_t chargedKaonSig = 0;
    };

    TH1D *hInclusiveTrigger = new TH1D("hInclusiveTrigger", "Inclusve transverse momentum spectrum for trigger hadrons", 100, 4, 50);
    TH1D *hInclusiveAssoc = new TH1D("hInclusiveAssoc", "Inclusve transverse momentum spectrum for associated hadrons", 100, 0, 20);
    TH1D *hTemp = new TH1D("template", "template", ndim + 1, 0, ndim + 1); // + 1 for K0_S/L
    TAxis *ax = hTemp->GetXaxis();
    ax->SetBinLabel(1, "K0_S/L");
    for (Int_t i = 0; i < ndim; i++){
        // 0th bin is underflow, 1st bin is K0_S/L, so start with i + 2
        ax->SetBinLabel(i + 2, hadron_vec[i]->getAntiName());
    } 
    std::unordered_map<Int_t, mapStruct> map; 
    for (Hadron *hadron : hadron_vec) {
        TH1D *hp = (TH1D*) hTemp->Clone();
        hp->SetName(hadron->getName());
        hp->SetOption("HIST");
        Int_t pdg = hadron->getPDG();
        Int_t abspdg = abs(pdg);
        Int_t strange;
        if(abspdg < 3300) strange = 1;
        else if(abspdg < 3330) strange = 2;
        else strange = 3;

        // normal pdg
        mapStruct normal;
        normal.name = hadron->getAntiName();
        normal.histogram = hp;
        normal.strangeness = strange;
        map[pdg] = normal;
        
        // anti pdg
        mapStruct anti;
        anti.name = hadron->getAntiName();
        anti.histogram = hp;
        anti.strangeness = -1*strange;
        map[-1*pdg] = anti;
    }
    
    Long64_t nentries = tree->GetEntries(); 
    Double_t pTDuplicateCheck = -1;
    for(Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        // inclusive spectra before kine cuts
        hInclusiveTrigger->Fill(pTTrigger);
        // check that the pT of the first assoc isn't exactly a match for the previous one
        // if it is, then we most likely have associated from the same event, resulting in duplicates in the inclusive spectrum
        // so we skip them
        if((*pTAssoc)[0] != pTDuplicateCheck){ 
            for(Double_t pT : *pTAssoc) hInclusiveAssoc->Fill(pT);
            pTDuplicateCheck = (*pTAssoc)[0]; 
        }

        // kine cuts before proceeding, for perfomance
        if(false) continue;
        // exception for K0_S/L
        if(pdgTrigger == Kzerolong->getPDG() || pdgTrigger == Kzeroshort->getPDG()) continue; 

        TH1D *h = map[pdgTrigger].histogram;
        Int_t tStrangeness = map[pdgTrigger].strangeness;
        map[pdgTrigger].ntriggers++;
        // loop over assoc
        for(Int_t pdg : *pdgAssoc){  
            // kine cuts
            if(false) continue;
            // exception for K0_S/L
            if(pdg == Kzerolong->getPDG() || pdg == Kzeroshort->getPDG()){
                h->Fill("K0_S/L", 0.5); // K0_S/L counts as half strange
                continue; // don't do anything else
            }

            Int_t aStrangeness = map[pdg].strangeness;
            // check if pair is ss or os and fill relevant histo
            if((tStrangeness > 0) != (aStrangeness > 0)){ // opposite sign
                h->Fill(map[pdg].name, abs(aStrangeness));
                // if K +/-, keep track so we can postprocess K0_S/L bkg
                if (abs(pdg) == Kminus->getAntiPDG()) map[pdgTrigger].chargedKaonSig++;
            } else if ((tStrangeness > 0) == (aStrangeness > 0)){ // same sign
                h->Fill(map[pdg].name, -1*abs(aStrangeness));
                // if K +/-, keep track so we can postprocess K0_S/L bkg
                if (abs(pdg) == Kminus->getAntiPDG()) map[pdgTrigger].chargedKaonBkg++;
            } else {
                std::cout << "wtf did you do???" << std::endl;
            }
        }
    }

    for (Hadron* hadron : hadron_vec){
        Int_t pdg = hadron->getPDG();
        Int_t antipdg = hadron->getAntiPDG();
        TH1D* h = map[pdg].histogram;
        // treat the neutral kaon bkg here:
        Double_t totBkg = map[pdg].chargedKaonBkg + map[antipdg].chargedKaonBkg;
        // std::cout << totBkg << std::endl;
        Double_t totSig = map[pdg].chargedKaonSig + map[antipdg].chargedKaonSig;
        // std::cout << totSig << std::endl;
        Double_t Bkg = totBkg*(h->GetBinContent(1)/totSig);
        // std::cout << Bkg << std::endl;
        h->Fill("K0_S/L", -1*Bkg);
        Int_t Ntotal = map[pdg].ntriggers + map[antipdg].ntriggers;
        Double_t ratio = h->Integral()/Ntotal;
        std::cout << hadron->getName() << " ratio: " << ratio << " found with ntriggers = " << Ntotal << std::endl;
    }
    
    inputfile->Close();
    outputfile->Write();
    outputfile->Close();

    return;
}