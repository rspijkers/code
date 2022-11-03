// std
#include <iostream>
#include <vector>
#include <set>
#include <map>
// ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TSystem.h"
// custom
#include "include/hadrons.h"
#include "include/helperfunctions.h"
#include "include/SmallEvent.h" // also includes SmallTrack 

// take a csv file with the lower bin edges and efficiency values. Note that the last entry should be the upper bin edge of the last bin. 
TH1D makeEfficiency(const char* filepath, Double_t BR = 1){
    // use vectors to get efficiency info from csv
    std::vector<Double_t> _BinEdges, _Efficiency;
    CSVtoXYArrays(filepath, &_BinEdges, &_Efficiency);
    const Int_t nbins = _Efficiency.size();
    // copy info to arrays so we can pass to TH1D
    Double_t BinEdges[nbins], Efficiency[nbins]; 
    std::copy(_BinEdges.begin(), _BinEdges.end(), BinEdges);
    std::copy(_Efficiency.begin(), _Efficiency.end(), Efficiency);
    TH1D hEfficiency = TH1D("hXiEff", "Efficiency for Xi", nbins - 1, BinEdges); 
    for(int i = 0; i < nbins; i++) hEfficiency.SetBinContent(i+1, BR*Efficiency[i]);
    return hEfficiency;
}

int main() {
    gSystem->Load("lib/libEvent.so");
    // TH1::SetDefaultSumw2(); // make sure errors are propagated in histo'ss
    
    // stbc input filepaths:
    // "/data/alice/rspijker/output/220926_Monash_pp_10M_14TeV_minpT0p15/"
    // "/data/alice/rspijker/output/220926_Monash_ppbar_10M_14TeV_minpT0p15/"
    // "/data/alice/rspijker/output/220926_Ropes_pp_10M_14TeV_minpT0p15/"
    // "/data/alice/rspijker/output/220926_SkandsMode2_pp_10M_14TeV_minpT0p15/"

    TChain* chain = new TChain("tree");
    chain->Add("PythiaEventGen/test.root");

    // Create list of files 
    TObjArray* fileElements = chain->GetListOfFiles(); // not actually a list of files, hence the TChainElement fuckery
    const Int_t nfiles = fileElements->GetEntries();
    TFile* fileList[nfiles];
    for(Int_t i = 0; i < nfiles; i++) { // create an ACTUAL fileList
        TChainElement* element = dynamic_cast<TChainElement*>((*fileElements)[i]);
        TFile* file;
        // check for nullpointers every step of the way
        if(element) file = new TFile(element->GetTitle());
        if(file) fileList[i] = file;
    }
    if(!fileList[0]){
        std::cout << "file is nullptr, something went wrong" << std::endl;
        return 1;
    }

    TFile* outputfile = new TFile("treev2eff.root", "RECREATE");

    struct XiStruct{
        TH1D* hSS;
        TH1D* hOS;
        Double_t nSS = 0;
        Double_t nOS = 0; 
    };

    // make efficiency
    std::vector<Double_t> _XiBinEdges, _XiEff;
    CSVtoXYArrays("efficiencies/XiMin.csv", &_XiBinEdges, &_XiEff);
    const Int_t effbins = _XiEff.size();
    Double_t XiBinEdges[effbins], XiEff[effbins]; // one more entry for edges: also the upper edge of last bin
    std::copy(_XiBinEdges.begin(), _XiBinEdges.end(), XiBinEdges);
    std::copy(_XiEff.begin(), _XiEff.end(), XiEff);
    TH1D* hXiEff = new TH1D("hXiEff", "Efficiency for Xi", effbins - 1, XiBinEdges); 
    for(int i = 0; i < effbins; i++) hXiEff->SetBinContent(i+1, 0.64*XiEff[i]); // 64% is the BR for Lambda to p + pion
    // TODO: handle under- and overflow: underflow is 0 probably good, but overflow > 0.

    // Make vector with all the strange hadrons
    std::vector<Hadron*> posStrangeHadrons;
    for(Hadron* h: StrangeHadrons) if(h->getStrangeness() > 0) posStrangeHadrons.push_back(h);
    Int_t nbins = posStrangeHadrons.size();

    std::map<Int_t, XiStruct> XiMinMap; 
    std::map<Int_t, XiStruct> XiZeroMap; 
    TH1D* hTempDPhi = new TH1D("hTempDPhi", "Overwrite this title", 32, -0.5*PI, 1.5*PI); // 32 bins is easy to rebin
    TH1D* hTempRatio = new TH1D("hTempRatio", "template", nbins, 0, nbins); 
    TAxis* ax = hTempRatio->GetXaxis();
    hTempRatio->SetXTitle("Associate hadron");
    Int_t binnr = 0;
    for(Hadron* assoc : posStrangeHadrons){
        Hadron* antiassoc = assoc->getAntiParticle();
        TString latex = assoc->getLatex();
        TString antilatex = antiassoc->getLatex();
        ax->SetBinLabel(binnr + 1, antilatex);
        binnr++;

        Int_t pdg = assoc->getPDG();
        XiStruct XiMin;
        XiStruct XiZero;

        TH1D* hSS = (TH1D*) hTempDPhi->Clone();
        hSS->SetName("Xi-"+assoc->getName()+"Dphi");
        hSS->SetTitle("#Xi^{-}(#Xi^{+}) - " + latex + "(" + antilatex + ") correlations");
        TH1D* hOS = (TH1D*) hTempDPhi->Clone();
        hOS->SetName("Xi-"+antiassoc->getName()+"Dphi");
        hOS->SetTitle("#Xi^{-}(#Xi^{+}) - " + antilatex + "(" + latex + ") correlations");

        XiMin.hSS = hSS;
        XiMin.hOS = hOS;
        XiMinMap[pdg] = XiMin;

        TH1D* hSS1 = (TH1D*) hTempDPhi->Clone();
        hSS1->SetName("Xi0"+assoc->getName()+"Dphi");
        hSS1->SetTitle("#Xi^{0}(#bar{#Xi^{0}}) - " + latex + "(" + antilatex + ") correlations");
        TH1D* hOS1 = (TH1D*) hTempDPhi->Clone();
        hOS1->SetName("Xi0"+antiassoc->getName()+"Dphi");
        hOS1->SetTitle("#Xi^{0}(#bar{#Xi^{0}}) - " + antilatex + "(" + latex + ") correlations");

        XiZero.hSS = hSS1;
        XiZero.hOS = hOS1;
        XiZeroMap[pdg] = XiZero;
    }

    // Analysis options
    const Bool_t doEff = true;
    const Double_t minpT = 0.15;
    const Double_t maxEtaTrigger = 2.0;
    const Double_t maxEtaAssoc = 3.0;
    const Double_t maxRapidity = 0.5;
    const Bool_t BkgScaling =  false; // do bkg scaling, needed to account for asymmetry in pp collisions
    std::set<Int_t> uniqueTriggerPDGs, uniqueAssocPDGs; // to keep track of possible PDG

    // Set branch address, so we can use the variables when we GetEntry()
    SmallEvent *event = new SmallEvent();
    chain->SetBranchAddress("event", &event);

    // initialize expensive objects
    SmallTrack XiCand;
    SmallTrack Assoc;
    std::vector<SmallTrack> cands = {};
    XiStruct* fillXi;
    std::map<Int_t, XiStruct>* XiMap;

    const Long64_t nEvents = chain->GetEntries(); 
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
        chain->GetEntry(iEvent);
        cands = event->getCandidates();
        Int_t nCands = cands.size();
        for(Int_t i = 0; i < nCands; i++){
            XiCand = cands[i];
            Int_t pdg = XiCand.getPDG();
            Int_t abspdg = abs(pdg);

            // look for Xi
            if(abspdg == Ximinus.getPDG()){
                XiMap = &XiMinMap;
            } else if (abspdg == Xizero.getPDG()){
                XiMap = &XiZeroMap;
            } else continue;
            // found a Xi, kine cuts
            // if (XiCand.getpT() < minpT || XiCand.geteta() > maxEtaTrigger) continue; 
            Double_t y = rapidityFromEta(XiCand.geteta(), XiCand.getpT(), 1.320);
            if (XiCand.getpT() < minpT || abs(y) > maxRapidity) continue; // Xi's are roughly 1.320 GeV
            // do eff cut here:
            Double_t eff = 1;
            if(doEff){
                eff *= hXiEff->GetBinContent(hXiEff->GetBin(XiCand.getpT()));
                if(0.5 < abs(y)) eff*= 2*(1-abs(y)); // we already know y < ymax
            }
            // keep track of strangeness
            Int_t SS;
            try{
                SS = StrangeHadronPDGMap.at(pdg)->getStrangeness();
            } catch (std::out_of_range){
                std::cout << "unknown pdg, skipping. pdg = " << pdg << std::endl;
                continue;
            }
            
            for(Int_t j = 0; j < nCands; j++){
                if(i == j) continue; // don't correlate with self
                Assoc = cands[j];
                // kine cuts 
                Double_t mass;
                try{
                    mass = StrangeHadronPDGMap.at(Assoc.getPDG())->getMass();
                } catch (std::out_of_range){
                    std::cout << "unknown pdg, skipping. pdg = " << Assoc.getPDG() << std::endl;
                    continue;
                }

                Double_t y2 = rapidityFromEta(Assoc.geteta(), Assoc.getpT(), mass);
                if(Assoc.getpT() < minpT || abs(y2) > maxRapidity) continue; 

                Int_t pdgAssoc = Assoc.getPDG();
                // skip K0's
                if(pdgAssoc == Kzeroshort.getPDG() || pdgAssoc == Kzerolong.getPDG()) continue;
                //efficiency in case of Ximin
                if(doEff && abs(pdgAssoc) == Ximinus.getPDG()) {
                    eff *= hXiEff->GetBinContent(hXiEff->GetBin(Assoc.getpT()));
                    if(0.5 < abs(y2)) eff*= 2*(1-abs(y2)); // we already know y < ymax
                }
                Double_t dphi = DeltaPhi(XiCand.getphi(), Assoc.getphi());
                try{
                    fillXi = &XiMap->at(pdgAssoc);
                } catch (std::out_of_range){
                    // if it's out of range, it's an anti strange assoc
                    fillXi = &XiMap->at(-1*pdgAssoc);
                    SS*=-1;
                }
                // fill correct dphi hist and update yield.
                if(SS >= 1){
                    fillXi->hSS->Fill(dphi, eff);
                    fillXi->nSS++;
                } else {
                    fillXi->hOS->Fill(dphi, eff);
                    fillXi->nOS++;
                }
            }
        }
    }
    // do ratio plots here:
    // Don't forget to handle the K0's
    TH1D* hXiMinRatio = (TH1D*) hTempRatio->Clone();
    hXiMinRatio->SetName("hXiMinRatio");
    hXiMinRatio->SetTitle("test");
    for(auto bla : XiMinMap){
        TH1D* OS = bla.second.hOS;
        TH1D* SS = bla.second.hSS;
        // scale with expected N events before resetting sumw2
        if(doEff){
            OS->Scale(100);
            SS->Scale(100);
        }
        // reset errors to be equal to sqrt(N), needed when filling with weights due to efficiency
        OS->GetSumw2()->Set(0);
        OS->Sumw2();
        SS->GetSumw2()->Set(0);
        SS->Sumw2();

        // TODO: make OS - SS plot (perhaps also compute integral?)
        TH1D* hSignal = (TH1D*) hTempDPhi->Clone();
        TString name = OS->GetName();
        hSignal->SetName(name + "_subtracted");
        hSignal->Add(OS, SS, 1, -1);

        // make ratio plot
        TString binname = StrangeHadronPDGMap.at(bla.first)->getAntiParticle()->getLatex();
        Double_t nOS = bla.second.nOS;
        Double_t nSS = bla.second.nSS;
        Double_t bincontent = nOS - nSS;
        Double_t error = sqrt(nOS + nSS);

        binnr = hXiMinRatio->Fill(binname, bincontent);
        hXiMinRatio->SetBinError(binnr, error);
    }

    outputfile->Write();
    outputfile->Close();
    return 0;
}