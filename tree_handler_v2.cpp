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
    TH1D hEfficiency = TH1D("hEff", "Efficiency", nbins - 1, BinEdges); 
    for(int i = 0; i < nbins; i++) hEfficiency.SetBinContent(i+1, BR*Efficiency[i]);
    // TODO: handle under- and overflow: underflow is 0 probably good, but overflow > 0.
    return hEfficiency;
}

int main() {

    // TODO: make inputfiles an input argument

    gSystem->Load("lib/libEvent.so");
    TH1::SetDefaultSumw2(); // make sure errors are propagated in histo'ss

    TChain* chain = new TChain("tree");
    chain->Add("~/alice/data/ModelStudyFeb/Monash_pp_100M_14TeV/*.root");

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
    TDirectory* yielddir = outputfile->mkdir("yield");
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

    const Int_t nMultiplicityBins = 7;
    Double_t multiplicityBinning[nMultiplicityBins+1] = {0, 50, 100, 200, 300, 400, 500, 1000};
    TH1D* hPDG = hSpectra->ProjectionX();
    hPDG->SetName("hPDG");
    TH1I* hncands = new TH1I("hNStrangeHadrons", "N strange particles per event", 100, 0, 100); // to be filled in event loop
    TH2I* hPDGMultiplicity = new TH2I("hPDGMultiplicity", "Yield of hadrons that pass kinematic selections;PDG;Ntracks", 8000, -4000, 4000, nMultiplicityBins, multiplicityBinning);
    TH1D* hStrangenessPerStrangeTrigger = new TH1D("hStrangenessPerStrangeTrigger", "Percent strangeness per trigger (+ = OS)", 37, -3-(1./12.), 3+(1./12.));
    TH1D* hStrangenessPerAntiStrangeTrigger = new TH1D("hStrangenessPerAntiStrangeTrigger", "Percent strangeness per anti trigger (+ = OS)", 37, -3-(1./12.), 3+(1./12.));


    // TODO setnames
    TH1D hXiEffobj = makeEfficiency("efficiencies/XiMin.csv", 0.64); // Lambda -> p+ K- ~ 64%
    TH1D* hXiEff = &hXiEffobj;
    hXiEff->SetName("hXiEff");
    TH1D hOmegaEffobj = makeEfficiency("efficiencies/Omegamin.csv", 0.68*0.64); //Omega -> Lambda K- ~ 68%
    TH1D* hOmegaEff = &hOmegaEffobj;
    hOmegaEff->SetName("hOmegaEff");

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

    TH2D* hTempDPhi = new TH2D("hTempDPhi", "Overwrite this title;#Delta#varphi;Ntracks", 64, -0.5*PI, 1.5*PI, nMultiplicityBins, multiplicityBinning); // 64 bins is easy to rebin
    TH1D* hTempRatio = new TH1D("hTempRatio", "Overwrite this title;Associate hadron;pairs", nbins, 0, nbins); 
    TAxis* ax = hTempRatio->GetXaxis();
    hTempRatio->SetXTitle("Associate hadron");
    // TODO: draw dotted gray line at 0 for the ratio plot
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
            hSS->SetTitle(triglatex + "(" + antitriglatex + ") - " + assoclatex + "(" + antiassoclatex + ") correlations");
            TH2D* hOS = (TH2D*) hTempDPhi->Clone();
            hOS->SetName(trigger->getName()+antiassoc->getName()+"Dphi");
            hOS->SetTitle(triglatex + "(" + antitriglatex + ") - " + antiassoclatex + "(" + assoclatex + ") correlations");

            datastruct.hSS = hSS;
            datastruct.hOS = hOS;
            assocmapdummy[pdg] = datastruct;
        }
        trig.second = assocmapdummy; // set the map
    }

    // Kinematics
    const Double_t minpT = 1.2;
    const Double_t maxEta = 4.;
    const Double_t maxY = 2.;

    // allow for different kinematic cuts between trigger and assoc
    const Double_t minpTTrigger = minpT;
    const Double_t minpTAssoc = minpT;
    const Double_t maxEtaTrigger = maxEta;
    const Double_t maxEtaAssoc = maxEta;
    const Double_t maxYTrigger = maxY;
    const Double_t maxYAssoc = maxY;

    // Analysis options
    const Bool_t doEff = false;
    const Bool_t doEta = true;
    const Bool_t doRapidity = false;
    // do either eta or rapidity cut, not both
    if(doEta && doRapidity){
        cout << "Error: you cannot both cut on eta and rapidity, please only choose one! Exiting..." << endl;
        return 1;
    }

    // TODO: const Bool_t BkgScaling =  false; // do bkg scaling, needed to account for asymmetry in pp collisions

    // Set branch address, so we can use the variables when we GetEntry()
    SmallEvent *event = new SmallEvent();
    chain->SetBranchAddress("event", &event);

    // initialize objects
    std::map<Hadron*, Double_t> nTriggers;
    for (auto t : triggermap) nTriggers[t.first] = 0;        
    SmallTrack trigger;
    SmallTrack assoc;
    std::vector<SmallTrack> cands = {};
    DataStruct* fillXi;
    std::map<Int_t, DataStruct>* XiMap;
    Double_t y,y2;

    const Long64_t nEvents = chain->GetEntries(); 
    for(Int_t iEvent = 0; iEvent < nEvents; iEvent++){
        chain->GetEntry(iEvent);
        Int_t multiplicity = event->getNtracks4p0();
        cands = event->getCandidates();
        Int_t nCands = cands.size();
        hncands->Fill(nCands);
        for(Int_t i = 0; i < nCands; i++){
            trigger = cands[i];
            Int_t pdg = trigger.getPDG();
            Int_t abspdg = abs(pdg);
            Hadron* triggerHadron;
            try{
                triggerHadron = StrangeHadronPDGMap.at(abspdg);
            } catch (std::out_of_range){
                std::cout << "unknown trigger pdg, skipping. pdg = " << abspdg << std::endl;
                continue;
            }

            // ad hoc check for xi or omega
            Bool_t isXi = true; // keep track of Xi or Omega
            if(triggerHadron == &Ximinus){
                XiMap = &triggermap[&Ximinus];
            } else if (triggerHadron == &Xizero){
                XiMap = &triggermap[&Xizero];
            } else if (triggerHadron == &Omegaminus){
                isXi = false;
                XiMap = &triggermap[&Omegaminus];
            } else continue;


            // Psuedorapidity
            if (doEta && (trigger.getpT() < minpTTrigger || trigger.getEta() > maxEtaTrigger)) continue; 
            // Rapidity
            if (doRapidity){
                Double_t massTrigger;
                try{
                    massTrigger = StrangeHadronPDGMap.at(pdg)->getMass();
                } catch (std::out_of_range){
                    std::cout << "unknown pdg, skipping. pdg = " << trigger.getPDG() << std::endl;
                    continue;
                }
                y = rapidityFromEta(trigger.getEta(), trigger.getpT(), massTrigger);
                if (trigger.getpT() < minpT || abs(y) > maxYTrigger) continue; 
            }

            // keep track of n_triggers for normalization
            // nTriggers.at(triggerHadron)++;
            hPDGMultiplicity->Fill(pdg, multiplicity);

            // do eff cut here (optional):
            Double_t eff = 1;
            if(doEff){
                if(isXi) {
                    eff *= hXiEff->GetBinContent(hXiEff->GetBin(trigger.getpT()));
                } else {
                    eff *= hOmegaEff->GetBinContent(hOmegaEff->GetBin(trigger.getpT()));
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

            //debug
            Double_t StrangeTrigger = (Double_t) SS;
            Double_t sPerT = 0;
            
            for(Int_t j = 0; j < nCands; j++){
                if(i == j) continue; // don't correlate with self
                assoc = cands[j];

                // Pseudorapidity
                if (doEta && (assoc.getpT() < minpTAssoc || assoc.getEta() > maxEtaAssoc)) continue; 
                // Rapidity
                if (doRapidity){
                    Double_t massAssoc;
                    try{
                        massAssoc = StrangeHadronPDGMap.at(assoc.getPDG())->getMass();
                    } catch (std::out_of_range){
                        std::cout << "unknown pdg, skipping. pdg = " << assoc.getPDG() << std::endl;
                        continue;
                    }
                    y2 = rapidityFromEta(assoc.getEta(), assoc.getpT(), massAssoc);
                    if(assoc.getpT() < minpT || abs(y2) > maxYAssoc) continue; 
                }

                Int_t pdgAssoc = assoc.getPDG();
                // skip K0's
                if(pdgAssoc == Kzeroshort.getPDG() || pdgAssoc == Kzerolong.getPDG()) continue;
                
                if(doEff) { // Efficiency 
                    if (abs(pdgAssoc) == Ximinus.getPDG()){
                        eff *= hXiEff->GetBinContent(hXiEff->GetBin(assoc.getpT()));
                    } else if (abs(pdgAssoc) == Omegaminus.getPDG()) {
                        eff *= hOmegaEff->GetBinContent(hOmegaEff->GetBin(assoc.getpT()));
                    } else continue;
                    if(0.5 < abs(y2)) eff*= 2*(1-abs(y2)); // we already know y < ymax
                }
                // calculate dphi
                Double_t dphi = DeltaPhi(trigger.getPhi(), assoc.getPhi());

                // duplicate trigger check: incase of for example Omega-Omega correlations, we only want to fill the histo once. 
                // However, we will encounter this combination twice: one time for each Omega. 
                Double_t doublecount = 1.;
                if(abspdg==abs(pdgAssoc)) doublecount = 0.5;

                // figure out which histogram to fill
                try{
                    fillXi = &XiMap->at(pdgAssoc);
                } catch (std::out_of_range){ // if it's out of range, it's an anti strange assoc
                    fillXi = &XiMap->at(-1*pdgAssoc);
                    SS*=-1;
                }
                // fill correct dphi hist and update yield.
                if(SS >= 1){
                    fillXi->hSS->Fill(dphi, multiplicity, eff*doublecount);
                    fillXi->nSS+=eff*doublecount;
                } else {
                    fillXi->hOS->Fill(dphi, multiplicity, eff*doublecount);
                    fillXi->nOS+=eff*doublecount;
                }
                sPerT += StrangeHadronPDGMap.at(pdgAssoc)->getStrangeness();
            } // assoc loop

            if(StrangeTrigger>0) {
                hStrangenessPerStrangeTrigger->Fill(-sPerT/StrangeTrigger);
            } else {
                hStrangenessPerAntiStrangeTrigger->Fill(-sPerT/StrangeTrigger);
            }
        } // trigger loop
    } // event loop

    // postprocessing
    for(auto& trig : triggermap){
        Hadron* trigger = trig.first;
        ratiodir->cd();
        TH1D* hRatio = (TH1D*) hTempRatio->Clone();
        hRatio->SetName("h"+trigger->getName()+"_ratio");
        hRatio->SetTitle("Correlations between " + trigger->getLatex() + " and assocates (number of pairs)");
        // Int_t pdg = trigger->getPDG(); 
        // Int_t minpdgbin = hPDGMultiplicity->GetXaxis()->FindBin(-pdg);
        // Int_t pospdgbin = hPDGMultiplicity->GetXaxis()->FindBin(pdg);
        // TH1D* projection = hPDGMultiplicity->ProjectionX();
        // Double_t totalnorm = projection->GetBinContent(minpdgbin) + projection->GetBinContent(pospdgbin);
        // cout << totalnorm << endl;
        // std::vector<Double_t> norm;
        std::vector<TH1D*> vMultiplicityRatios;
        for(Int_t i = 0; i < nMultiplicityBins; i++){
            // norm.push_back(hPDGMultiplicity->GetBinContent(minpdgbin,i+1) + hPDGMultiplicity->GetBinContent(pospdgbin,i+1));
            TString loweredge = std::to_string((int) multiplicityBinning[i]);
            TString upperedge = std::to_string((int) multiplicityBinning[i+1]);
            TH1D* h = (TH1D*) hTempRatio->Clone();
            TString name = TString(trigger->getName()+"_"+loweredge+"-"+upperedge+"_N_part");
            TString title = TString("Correlations between " + trigger->getLatex() + " and assocates, "+loweredge+" < N_parts < "+upperedge);
            h->SetName(name);
            h->SetTitle(title);
            vMultiplicityRatios.push_back(h);
        }
        outputfile->cd(trigger->getName()+"/subtracted");
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
            
            // make signal plots
            TH2D* hSignal = (TH2D*) hTempDPhi->Clone();
            TString name = OS->GetName();
            TString title = OS->GetTitle();
            hSignal->SetName(name + "_subtracted");
            hSignal->SetTitle(title + " (subtracted)");
            hSignal->Add(OS, SS, 1, -1);

            // make total ratio plot
            Hadron* antiAssoc = StrangeHadronPDGMap.at(bla.first)->getAntiParticle();
            TString binname = antiAssoc->getLatex();
            Double_t nOS = bla.second.nOS;
            Double_t nSS = bla.second.nSS;
            Double_t yield = nOS - nSS;
            Double_t error = sqrt(nOS + nSS);
            Int_t binnr = hRatio->Fill(binname, yield);
            hRatio->SetBinError(binnr, error);

            // yields per multiplicity
            TH1D* hYield = new TH1D("hYield", "Yield", nMultiplicityBins, multiplicityBinning);
            hYield->SetTitle("Yields of "+trigger->getLatex()+" - "+antiAssoc->getLatex()+" pairs");
            hYield->SetName(trigger->getName()+"-"+antiAssoc->getName());
            hYield->SetDirectory(yielddir);
            // So with this 2D signal histo, we can make ratioplots per multiplicitybin
            for(Int_t i = 0; i < nMultiplicityBins; i++){
                Int_t nxbins = hSignal->GetNbinsX();
                yield = hSignal->IntegralAndError(0, nxbins, i, i+1, error);
                binnr = vMultiplicityRatios[i]->Fill(binname, yield);
                vMultiplicityRatios[i]->SetBinError(binnr, error);
                hYield->SetBinContent(i+1, yield);
                hYield->SetBinError(i+1, error);
            } // end multiplicity loop
        } // end assoc loop
        // normalize histograms
        hRatio->Scale(1/hRatio->Integral()); 
        for (TH1D* h : vMultiplicityRatios){
            h->Scale(1/h->Integral()); 
        }
    } // end trigger loop

    outputfile->Write();
    outputfile->Close();
    return 0;
}
