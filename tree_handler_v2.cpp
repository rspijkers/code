// std
#include <iostream>
#include <vector>
#include <set>
#include <map>
// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TSystem.h"
// custom
#include "include/hadrons.h"
#include "include/helperfunctions.h"
#include "include/SmallEvent.h" // also includes SmallTrack 

using std::cout; using std::endl;

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
    // TODO: handle under- and overflow: underflow is 0 probably good, but overflow > 0.
    return hEfficiency;
}

int main() {

    // make inputfiles an input argument

    gSystem->Load("lib/libEvent.so");
    // TH1::SetDefaultSumw2(); // make sure errors are propagated in histo'ss
    
    // stbc input filepaths:
    // "/data/alice/rspijker/output/220926_Monash_pp_10M_14TeV_minpT0p15/"
    // "/data/alice/rspijker/output/220926_Monash_ppbar_10M_14TeV_minpT0p15/"
    // "/data/alice/rspijker/output/220926_Ropes_pp_10M_14TeV_minpT0p15/"
    // "/data/alice/rspijker/output/220926_SkandsMode2_pp_10M_14TeV_minpT0p15/"

    TChain* chain = new TChain("tree");
    chain->Add("~/alice/data/ModelStudyJan/Monash_pp_100M_14TeV/*.root");

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

    // outputfile structure
    TFile* outputfile = new TFile("treev2eff.root", "RECREATE");
    TDirectory* XiMinusdir = outputfile->mkdir(Ximinus.getName());
    TDirectory* XiZerodir = outputfile->mkdir(Xizero.getName());
    TDirectory* Omegadir = outputfile->mkdir(Omegaminus.getName());
    TDirectory* ratiodir = outputfile->mkdir("ratios");
    TDirectory* XiMinusSubdir = outputfile->mkdir(Ximinus.getName()+"/subtracted");
    TDirectory* XiZeroSubdir = outputfile->mkdir(Xizero.getName()+"/subtracted");
    TDirectory* OmegaSubdir = outputfile->mkdir(Omegaminus.getName()+"/subtracted");

    // Do QA from spectra histo (PDG, pT, eta):
    TH3D* hSpectra = (TH3D*) fileList[0]->Get<TH3D>("hSpectra")->Clone();
    for(Int_t i = 1; i < nfiles; i++){
        TH3D* h_ = fileList[i]->Get<TH3D>("hSpectra");
        if(!h_){ // check h_ is not nullpointer
            std::cout << "histogram (or file) is nullptr, something went wrong" << std::endl;
            return 1;
        } 
        hSpectra->Add(h_);
    }
    TH1D* hPDG = hSpectra->ProjectionX();
    hPDG->SetName("hPDG");

    // TODO setnames
    TH1D hXiEffobj = makeEfficiency("efficiencies/XiMin.csv", 0.64);
    TH1D* hXiEff = &hXiEffobj;

    TH1D hOmegaEffobj = makeEfficiency("efficiencies/Omegamin.csv", 0.68*0.64); //Omega -> Lambda K- ~ 68%
    TH1D* hOmegaEff = &hOmegaEffobj;

    // Make vector with all the strange hadrons
    std::vector<Hadron*> posStrangeHadrons;
    for(Hadron* h: StrangeHadrons) if(h->getStrangeness() > 0) posStrangeHadrons.push_back(h);
    Int_t nbins = posStrangeHadrons.size();

    // Make analysis structures that hold the relevant data
    struct DataStruct{
        TH2D* hSS;
        TH2D* hOS;
        Double_t nSS = 0;
        Double_t nOS = 0; 
    };
    using assocmap = std::map<Int_t, DataStruct>;
    assocmap XiMinMap; 
    assocmap XiZeroMap; 
    assocmap OmegaMap;
    
    // array of interesting triggers
    std::map<Hadron*, assocmap> triggermap = {{&Ximinus, XiMinMap}, {&Xizero, XiZeroMap}, {&Omegaminus, OmegaMap}};

    Double_t multiplicityBinning[7] = {0, 50, 200, 350, 500, 650, 1000};
    TH2D* hTempDPhi = new TH2D("hTempDPhi", "Overwrite this title;#Delta#varphi;Ntracks", 32, -0.5*PI, 1.5*PI, 6, multiplicityBinning); // 32 bins is easy to rebin
    // TH2D* hTempRatio = new TH2D("hTempRatio", "Overwrite this title;Associate hadron;Multiplicity", nbins, 0, nbins, 6, multiplicityBinning); 
    TH1D* hTempRatio = new TH1D("hTempRatio", "Overwrite this title;Associate hadron;Ntracks", nbins, 0, nbins); 
    TAxis* ax = hTempRatio->GetXaxis();
    hTempRatio->SetXTitle("Associate hadron");
    assocmap assocmapdummy; // needs to be defined before the loop, otherwise it won't survive the loop
    for (auto& trig : triggermap){
        Hadron* trigger = trig.first; 
        Int_t binnr = 0;
        outputfile->cd(trigger->getName());
        DataStruct datastruct;
        for(Hadron* assoc : posStrangeHadrons){
            // names and titles
            Hadron* antiassoc = assoc->getAntiParticle();
            Hadron* antitrigger = trigger->getAntiParticle();
            TString assoclatex = assoc->getLatex();
            TString antiassoclatex = antiassoc->getLatex();
            TString triglatex = trigger->getLatex();
            TString antitriglatex = antitrigger->getLatex();

            ax->SetBinLabel(binnr + 1, antiassoclatex);
            binnr++;
            Int_t pdg = assoc->getPDG();

            TH2D* hSS = (TH2D*) hTempDPhi->Clone();
            hSS->SetName(trigger->getName()+assoc->getName()+"Dphi");
            hSS->SetTitle(triglatex + "(" + antitriglatex + ")" + assoclatex + "(" + antiassoclatex + ") correlations");
            TH2D* hOS = (TH2D*) hTempDPhi->Clone();
            hOS->SetName(trigger->getName()+antiassoc->getName()+"Dphi");
            hOS->SetTitle(triglatex + "(" + antitriglatex + ")" + antiassoclatex + "(" + assoclatex + ") correlations");

            datastruct.hSS = hSS;
            datastruct.hOS = hOS;
            assocmapdummy[pdg] = datastruct;
        }
        trig.second = assocmapdummy; // set the map
    }

    // Analysis options
    const Bool_t doEff = false;
    const Double_t minpT = 1.2;
    const Double_t maxEtaTrigger = 2.;
    const Double_t maxEtaAssoc = 2.;
    const Double_t maxRapidity = 2.;
    // TODO: const Bool_t BkgScaling =  false; // do bkg scaling, needed to account for asymmetry in pp collisions
    std::set<Int_t> uniqueTriggerPDGs, uniqueAssocPDGs; // to keep track of possible PDG

    // Set branch address, so we can use the variables when we GetEntry()
    SmallEvent *event = new SmallEvent();
    chain->SetBranchAddress("event", &event);

    // initialize expensive objects
    SmallTrack XiCand;
    SmallTrack Assoc;
    std::vector<SmallTrack> cands = {};
    DataStruct* fillXi;
    std::map<Int_t, DataStruct>* XiMap;

    // TODO: make nice, this is to debug
    TH1I* hncands = new TH1I("hncands", "N strange particles per event", 100, 0, 100);

    const Long64_t nEvents = chain->GetEntries(); 
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
        chain->GetEntry(iEvent);
        Int_t multiplicity = event->getNtracks();
        cands = event->getCandidates();
        Int_t nCands = cands.size();
        hncands->Fill(nCands);
        for(Int_t i = 0; i < nCands; i++){
            XiCand = cands[i];
            Int_t pdg = XiCand.getPDG();
            Int_t abspdg = abs(pdg);

            // ad hoc check for xi or omega
            Bool_t isXi = true; 
            if(abspdg == Ximinus.getPDG()){
                // XiMap = &XiMinMap;
                XiMap = &triggermap[&Ximinus];
            } else if (abspdg == Xizero.getPDG()){
                // XiMap = &XiZeroMap;
                XiMap = &triggermap[&Xizero];
            } else if (abspdg == Omegaminus.getPDG()){
                isXi = false;
                // XiMap = &OmegaMap;
                XiMap = &triggermap[&Omegaminus];
            } else continue;

            // Kine cuts: do either rapidity or pseudorapidity cut ( + pt cut)
            // if (XiCand.getpT() < minpT || XiCand.geteta() > maxEtaTrigger) continue; 
            Double_t y = rapidityFromEta(XiCand.getEta(), XiCand.getpT(), 1.320);
            if (XiCand.getpT() < minpT || abs(y) > maxRapidity) continue; // Xi's are roughly 1.320 GeV
            // do eff cut here (optional):
            Double_t eff = 1;
            if(doEff){
                if(isXi) {
                    eff *= hXiEff->GetBinContent(hXiEff->GetBin(XiCand.getpT()));
                } else {
                    eff *= hOmegaEff->GetBinContent(hOmegaEff->GetBin(XiCand.getpT()));
                }
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

                Double_t y2 = rapidityFromEta(Assoc.getEta(), Assoc.getpT(), mass);
                if(Assoc.getpT() < minpT || abs(y2) > maxRapidity) continue; 

                Int_t pdgAssoc = Assoc.getPDG();
                // skip K0's
                if(pdgAssoc == Kzeroshort.getPDG() || pdgAssoc == Kzerolong.getPDG()) continue;
                //efficiency in case of Ximin or Omega
                if(doEff) {
                    if (abs(pdgAssoc) == Ximinus.getPDG()){
                        eff *= hXiEff->GetBinContent(hXiEff->GetBin(Assoc.getpT()));
                    } else if (abs(pdgAssoc) == Omegaminus.getPDG()) {
                        eff *= hOmegaEff->GetBinContent(hOmegaEff->GetBin(Assoc.getpT()));
                    } else continue;
                    if(0.5 < abs(y2)) eff*= 2*(1-abs(y2)); // we already know y < ymax
                }
                Double_t dphi = DeltaPhi(XiCand.getPhi(), Assoc.getPhi());
                try{
                    fillXi = &XiMap->at(pdgAssoc);
                } catch (std::out_of_range){
                    // if it's out of range, it's an anti strange assoc
                    fillXi = &XiMap->at(-1*pdgAssoc);
                    SS*=-1;
                }
                // fill correct dphi hist and update yield.
                if(SS >= 1){
                    fillXi->hSS->Fill(dphi, multiplicity, eff);
                    fillXi->nSS++;
                } else {
                    fillXi->hOS->Fill(dphi, multiplicity, eff);
                    fillXi->nOS++;
                }
            } // 2nd cand loop
        } // 1st cand loop
    } // event loop

    // do ratio plots here:
    // Let's try to do it with a loop over the triggers
    Int_t nMultiplicityBins = hTempDPhi->GetNbinsY();
    cout << nMultiplicityBins << endl;
    for(auto& trig : triggermap){
        // Don't forget to handle the K0's
        Hadron* trigger = trig.first;
        outputfile->cd(trigger->getName()+"/subtracted");
        TH1D* hRatio = (TH1D*) hTempRatio->Clone();
        hRatio->SetName("h"+trigger->getName()+"_ratio");
        hRatio->SetTitle("Correlations between " + trigger->getLatex() + " and assocates (number of pairs)");
        hRatio->SetDirectory(ratiodir);
        std::vector<TH1D*> vMultiplicityRatios;
        for(Int_t i = 0; i < nMultiplicityBins; i++){
            TString loweredge = std::to_string((int) multiplicityBinning[i]);
            TString upperedge = std::to_string((int) multiplicityBinning[i+1]);
            TH1D* h = (TH1D*) hTempRatio->Clone();
            TString name = TString(loweredge+"-"+upperedge+" N_part");
            h->SetName(name);
            vMultiplicityRatios.push_back(h);
        }
        for(auto bla : trig.second){
            TH2D* OS = bla.second.hOS;
            TH2D* SS = bla.second.hSS;
            // scale with expected N events before resetting sumw2
            if(doEff){ // FIXME
                OS->Scale(100);
                SS->Scale(100);
            }

            // reset errors to be equal to sqrt(N), needed when filling with weights due to efficiency
            // doesn't change anything when running without efficiency
            OS->GetSumw2()->Set(0);
            OS->Sumw2();
            SS->GetSumw2()->Set(0);
            SS->Sumw2();
            
            // TODO: make OS - SS plot (perhaps also compute integral?)
            TH2D* hSignal = (TH2D*) hTempDPhi->Clone();
            TString name = OS->GetName();
            TString title = OS->GetTitle();
            hSignal->SetName(name + "_subtracted");
            hSignal->SetTitle(title + " (subtracted)");
            hSignal->Add(OS, SS, 1, -1);

            // make total ratio plot
            TString binname = StrangeHadronPDGMap.at(bla.first)->getAntiParticle()->getLatex();
            Double_t nOS = bla.second.nOS;
            Double_t nSS = bla.second.nSS;
            Double_t bincontent = nOS - nSS;
            Double_t error = sqrt(nOS + nSS);
            // perhaps error should be sqrt(2*nSS+nOS-nSS)

            Int_t binnr = hRatio->Fill(binname, bincontent);
            hRatio->SetBinError(binnr, error);

            // So with this 2D signal histo, we can make ratioplots per multiplicitybin
            for(Int_t ybin = 0; ybin < nMultiplicityBins; ybin++){
                Double_t integral, error;
                Int_t nxbins = hSignal->GetNbinsX();
                integral = hSignal->IntegralAndError(0, nxbins, ybin, ybin+1, error);
                binnr = vMultiplicityRatios[ybin]->Fill(binname, integral);
                vMultiplicityRatios[ybin]->SetBinError(binnr, error);
            }
        }
    }

    outputfile->Write();
    outputfile->Close();
    return 0;
}
