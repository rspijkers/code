#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TPad.h"
#include <vector>
#include <set>
#include <unordered_map>
#include "hadrons.h"
#include "helperfunctions.h"

#define PI 3.14159265

void tree_handler() {
    TH1::SetDefaultSumw2(); // make sure errors are propagated in histo's
    
    TFile* inputfile = new TFile("500ktest.root", "READ");
    TTree* tree = (TTree*) inputfile->Get("tree");
    // tree->Print();
    // tree->Show(0);

    // Check to see if the under/overflow bins in the pdg are empty. If not, smt funky is going on
    TH1D* hPDGTrigger = (TH1D*) inputfile->Get("hPDG");
    Int_t nPDGbins = hPDGTrigger->GetNbinsX();
    if(hPDGTrigger->GetBinContent(0) > 0 || hPDGTrigger->GetBinContent(nPDGbins + 1) > 0){
        std::cout << "WARNING: non-empty underflow and/or overflow bins in the PDG histogram!!! This means we are not accounting for all particle species. Please investigate" << std::endl;
    }

    TFile* outputfile = new TFile("test_handler.root", "RECREATE");

    // Kinematic cuts
    const Double_t maxEtaTrigger = 2.0;
    const Double_t maxEtaAssoc = 3.0;
    // std::cout << "We are using eta cuts: maxEtaTrigger = " << maxEtaTrigger << " and maxEtaAssoc = " << maxEtaAssoc << std::endl;
    std::set<Int_t> uniqueTriggerPDGs, uniqueAssocPDGs; // to keep track of possible PDG

    // Set branch address, so we can use the variables when we GetEntry()
    Int_t pdgTrigger;
    std::vector<Int_t>* pdgAssoc = 0;
    Double_t pTTrigger, etaTrigger;
    std::vector<Double_t>* pTAssoc = 0;
    std::vector<Double_t>* etaAssoc = 0;
    std::vector<Double_t>* deltaPhi = 0;
    std::vector<Double_t>* deltaEta = 0;
    tree->SetBranchAddress("pdgTrigger", &pdgTrigger);
    tree->SetBranchAddress("pdgAssoc", &pdgAssoc);
    tree->SetBranchAddress("pTTrigger", &pTTrigger);
    tree->SetBranchAddress("etaTrigger", &etaTrigger);
    tree->SetBranchAddress("pTAssoc", &pTAssoc);
    tree->SetBranchAddress("etaAssoc", &etaAssoc);
    tree->SetBranchAddress("deltaPhi", &deltaPhi);
    tree->SetBranchAddress("deltaEta", &deltaEta);

    std::vector<Hadron> hadron_vec = {Kminus, Lambda, Sigmaminus, Sigmazero, Sigmaplus, Ximinus, Xizero, Omegaminus, Dsubs, Bsubs, Xicplus, Xiczero, Xibmin, Xibzero, Omegac, Omegab, Omegacc, Omegabb, Omegabc};
    const Int_t nbins = hadron_vec.size() + 1; // + 1 for K0_S/L
    
    struct mapStruct
    {
        Hadron hadron;
        TH1D* hSig; 
        TH1D* hBkg; 
        Int_t strangeness;

        Int_t ntriggers = 0;
        Double_t chargedKaonBkg = 0;
        Double_t chargedKaonSig = 0;
    };

    // QA & Kinematic plots
    TH1D* hInclusiveTrigger = new TH1D("hInclusiveTrigger", "Inclusive transverse momentum spectrum for trigger hadrons", 100, 0, 25);
    TH1D* hInclusiveAssoc = new TH1D("hInclusiveAssoc", "Inclusive transverse momentum spectrum for associated hadrons", 100, 0, 10);
    TH1D* hEtaTrigger = new TH1D("hEtaTrigger", "Pseudorapidity spectrum for trigger hadrons", 100, -10, 10);
    TH1D* hEtaAssoc = new TH1D("hEtaAssoc", "Pseudorapidity spectrum for associated hadrons", 100, -10, 10);
    TH1D* hYTrigger = new TH1D("hYTrigger", "Rapidity spectrum for trigger hadrons", 100, -10, 10);
    TH1D* hYAssoc = new TH1D("hYAssoc", "Rapidity spectrum for associated hadrons", 100, -10, 10);
    TH1D* hEtaNormal = new TH1D("hEtaNormal", "Pseudorapidity spectrum for 'normal' baryons", 100, -10, 10);
    TH1D* hEtaAnti = new TH1D("hEtaAnti", "Pseudorapidity spectrum for anti baryons", 100, -10, 10);

    // DeltaPhi plots
    TH1D* hLLdphi = new TH1D("hLLdphi", "deltaphi for L-L", 50, -0.5*PI, 1.5*PI);
    TH1D* hLLbardphi = new TH1D("hLLbardphi", "deltaphi for L-Lbar", 50, -0.5*PI, 1.5*PI);
    TH1D* hLbarLdphi = new TH1D("hLbarLdphi", "deltaphi for Lbar-L", 50, -0.5*PI, 1.5*PI);
    TH1D* hLbarLbardphi = new TH1D("hLbarLbardphi", "deltaphi for Lbar-Lbar", 50, -0.5*PI, 1.5*PI);

    TH1D* hLKdphi = new TH1D("hLKdphi", "deltaphi for L-Kmin", 50, -0.5*PI, 1.5*PI);
    TH1D* hLKbardphi = new TH1D("hLKbardphi", "deltaphi for L-Kplus", 50, -0.5*PI, 1.5*PI);
    TH1D* hLbarKdphi = new TH1D("hLbarKdphi", "deltaphi for Lbar-Kmin", 50, -0.5*PI, 1.5*PI);
    TH1D* hLbarKbardphi = new TH1D("hLbarKbardphi", "deltaphi for Lbar-Kplus", 50, -0.5*PI, 1.5*PI);

    TH1D* hLSdphi = new TH1D("hLSdphi", "deltaphi for L-Sigmamin", 50, -0.5*PI, 1.5*PI);
    TH1D* hLSbardphi = new TH1D("hLSbardphi", "deltaphi for L-Sigmaminbar", 50, -0.5*PI, 1.5*PI);
    TH1D* hLbarSdphi = new TH1D("hLbarSdphi", "deltaphi for Lbar-Sigmamin", 50, -0.5*PI, 1.5*PI);
    TH1D* hLbarSbardphi = new TH1D("hLbarSbardphi", "deltaphi for Lbar-Sigmaminbar", 50, -0.5*PI, 1.5*PI);

    TH1D* hKKdphi = new TH1D("hKKdphi", "deltaphi for Kmin-Kmin", 50, -0.5*PI, 1.5*PI);
    TH1D* hKKbardphi = new TH1D("hKKbardphi", "deltaphi for Kmin-Kplus", 50, -0.5*PI, 1.5*PI);
    TH1D* hKbarKdphi = new TH1D("hKbarKdphi", "deltaphi for Kplus-Kmin", 50, -0.5*PI, 1.5*PI);
    TH1D* hKbarKbardphi = new TH1D("hKbarKbardphi", "deltaphi for Kplus-Kplus", 50, -0.5*PI, 1.5*PI);


    TH1D* hTemp = new TH1D("template", "template", nbins, 0, nbins); 
    // hTemp->SetOption("HIST E");
    TAxis* ax = hTemp->GetXaxis();
    // TODO: use LaTeX for bin labels? and should we do the normal, or antiname? or smt else?
    ax->SetBinLabel(1, "K^{0}_{S/L}");
    for (Int_t i = 0; i < nbins - 1; i++){
        // 0th bin is underflow, 1st bin is K0_S/L, so start with i + 2
        ax->SetBinLabel(i + 2, hadron_vec[i].getAntiLatex());
    } 
    std::unordered_map<Int_t, mapStruct> map; 
    for (Hadron hadron : hadron_vec) {
        TString name = hadron.getName();
        TString antiname = hadron.getAntiName();

        TH1D* hss = (TH1D*) hTemp->Clone();
        hss->SetName(name + "_bkg");
        hss->SetTitle("strange strange pairs");
        TH1D* hssbar = (TH1D*) hTemp->Clone();
        hssbar->SetName(name + "_sig");
        hssbar->SetTitle("strange anti-strange pairs");
        TH1D* hsbars = (TH1D*) hTemp->Clone();
        hsbars->SetName(antiname + "_sig");
        hsbars->SetTitle("anti-strange strange pairs");
        TH1D* hsbarsbar = (TH1D*) hTemp->Clone();
        hsbarsbar->SetName(antiname + "_bkg");
        hsbarsbar->SetTitle("anti-strange anti-strange pairs");
        
        Int_t pdg = hadron.getPDG();
        Int_t abspdg = abs(pdg);
        Int_t strange = strangenessFromPDG(pdg);

        // normal pdg
        mapStruct normal;
        normal.hadron = hadron;
        normal.hSig = hssbar;
        normal.hBkg = hss;
        normal.strangeness = strange;
        map[pdg] = normal;

        // anti pdg
        mapStruct anti;
        anti.hadron = hadron;
        anti.hSig = hsbars;
        anti.hBkg = hsbarsbar;
        anti.strangeness = -1*strange;
        map[-1*pdg] = anti;
    }
    // hTemp->Delete();
    mapStruct k0shortstruct, k0longstruct;
    k0shortstruct.hadron = Kzeroshort;
    k0longstruct.hadron = Kzerolong;
    map[Kzeroshort.getPDG()] = k0shortstruct;
    map[Kzerolong.getPDG()] = k0longstruct;

    const Long64_t nentries = tree->GetEntries(); 
    Double_t pTDuplicateCheck = -1;
    for(Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        // inclusive spectra before kine cuts
        uniqueTriggerPDGs.insert(pdgTrigger); // will only insert if unique
        Hadron trigger;
        // try to access element of map. if doesn't exist, skip this entire iteration
        try{
            trigger = map.at(pdgTrigger).hadron;
        } catch (std::out_of_range){
            std::cout << "Unknown trigger hadron! pdg = " << pdgTrigger << ". Skipping this iteration..." << std::endl;
            continue;
        }
        hInclusiveTrigger->Fill(pTTrigger);
        hEtaTrigger->Fill(etaTrigger);
        Double_t y = rapidityFromEta(etaTrigger, pTTrigger, trigger.getMass()/1000.); // convert mass to GeV
        hYTrigger->Fill(y);
        // check that the pT of the first assoc isn't exactly a match for the previous one
        // if it is, then we most likely have associated from the same event, resulting in duplicates in the inclusive spectrum
        // so we skip them
        Bool_t isDuplicate = true;
        if((*pTAssoc)[0] != pTDuplicateCheck){ 
            isDuplicate = false;
            pTDuplicateCheck = (*pTAssoc)[0]; 
        }

        // kine cuts before proceeding, for perfomance
        // if(etaTrigger > maxEtaTrigger) continue;
        // exception for K0_S/L
        if(pdgTrigger == Kzerolong.getPDG() || pdgTrigger == Kzeroshort.getPDG()) continue; 
        // || pdgTrigger == Kminus.getPDG()  || pdgTrigger == Kminus.getAntiPDG()
        if(pdgTrigger > 1000) {hEtaNormal->Fill(etaTrigger);}
        else if(pdgTrigger < -1000) {hEtaAnti->Fill(etaTrigger);}

        TH1D* hSig = map.at(pdgTrigger).hSig;
        TH1D* hBkg = map.at(pdgTrigger).hBkg;
        Int_t tStrangeness = map.at(pdgTrigger).strangeness;
        map.at(pdgTrigger).ntriggers++;
        // loop over assoc
        for(Int_t j = 0; j < pdgAssoc->size(); j++){
            Int_t pdg = (*pdgAssoc)[j];
            Double_t eta = (*etaAssoc)[j];
            Double_t pt = (*pTAssoc)[j];
            Hadron assoc;
            try{
                assoc = map.at(pdg).hadron;
            } catch (std::out_of_range){
                std::cout << "Unknown associated hadron! pdg = " << pdg << ". Skipping this iteration..." << std::endl;
                continue;
            }
            if(!isDuplicate){
                uniqueAssocPDGs.insert(pdg); // will only insert if unique
                hInclusiveAssoc->Fill(pt);
                hEtaAssoc->Fill(eta);
                y = rapidityFromEta(eta, pt, assoc.getMass()/1000.); // convert mass to GeV
                hYAssoc->Fill(y);
            }

            // kine cuts
            // if(eta > maxEtaAssoc) continue;
            // exception for K0_S/L
            if(pdg == Kzerolong.getPDG() || pdg == Kzeroshort.getPDG()){
                // we use the number of pairs with K0_S/L to estimate the error therefor we fill the bkg histo
                // the signal of K0_S/L is calculated by using conservation of strangeness, we don't need to do fill any signal histo here
                hBkg->Fill("K^{0}_{S/L}", 0.5); // K0_S/L counts as half strange
                continue; // don't do anything else
            }

            Int_t aStrangeness = map.at(pdg).strangeness;
            // check if pair is ss or os and fill relevant histo
            if((tStrangeness > 0) != (aStrangeness > 0)){ // opposite sign
                hSig->Fill(assoc.getAntiName(), abs(aStrangeness));
                // if K +/-, keep track so we can postprocess K0_S/L bkg
                if (abs(pdg) == Kminus.getAntiPDG()) map.at(pdgTrigger).chargedKaonSig++;
            } else if ((tStrangeness > 0) == (aStrangeness > 0)){ // same sign
                hBkg->Fill(assoc.getAntiName(), abs(aStrangeness));
                // if K +/-, keep track so we can postprocess K0_S/L bkg
                if (abs(pdg) == Kminus.getAntiPDG()) map.at(pdgTrigger).chargedKaonBkg++;
            } else {
                std::cout << "wtf did you do???" << std::endl;
            }

            // make some deltaphi plots for lambda correlations
            if(pdgTrigger == Lambda.getPDG()){
                if(pdg == Lambda.getPDG()){hLLdphi->Fill((*deltaPhi)[j]);} // LL
                else if(pdg == Lambda.getAntiPDG()){hLLbardphi->Fill((*deltaPhi)[j]);} // LLbar
                else if(pdg == Kminus.getPDG()){hLKdphi->Fill((*deltaPhi)[j]);} // LK
                else if(pdg == Kminus.getAntiPDG()){hLKbardphi->Fill((*deltaPhi)[j]);} // LKbar
                else if(pdg == Sigmaminus.getPDG()){hLSdphi->Fill((*deltaPhi)[j]);} // LS
                else if(pdg == Sigmaminus.getAntiPDG()){hLSbardphi->Fill((*deltaPhi)[j]);} // LSbar
            } else if(pdgTrigger == Lambda.getAntiPDG()){
                if(pdg == Lambda.getPDG()){hLbarLdphi->Fill((*deltaPhi)[j]);} // LbarL
                else if(pdg == Lambda.getAntiPDG()){hLbarLbardphi->Fill((*deltaPhi)[j]);} // LbarLbar
                else if(pdg == Kminus.getPDG()){hLbarKdphi->Fill((*deltaPhi)[j]);} // LbarK
                else if(pdg == Kminus.getAntiPDG()){hLbarKbardphi->Fill((*deltaPhi)[j]);} // LbarKbar
                else if(pdg == Sigmaminus.getPDG()){hLbarSdphi->Fill((*deltaPhi)[j]);} // LbarS
                else if(pdg == Sigmaminus.getAntiPDG()){hLbarSbardphi->Fill((*deltaPhi)[j]);} // LbarSbar
            }
            // Kaon corr.
            else if(pdgTrigger == Kminus.getPDG()){
                if(pdg == Kminus.getPDG()){hKKdphi->Fill((*deltaPhi)[j]);} // KminKmin
                else if(pdg == Kminus.getAntiPDG()){hKKbardphi->Fill((*deltaPhi)[j]);} // KminKplus
            } else if(pdgTrigger == Kminus.getAntiPDG()){
                if(pdg == Kminus.getPDG()){hKbarKdphi->Fill((*deltaPhi)[j]);} // KplusKmin
                else if(pdg == Kminus.getAntiPDG()){hKbarKbardphi->Fill((*deltaPhi)[j]);} // KplusKplus
            }
        }
    }

    // Bkg, normalization, and errors
    for (Hadron hadron : hadron_vec){
        Int_t pdg = hadron.getPDG();
        Int_t antipdg = hadron.getAntiPDG();
        mapStruct normal = map.at(pdg);
        mapStruct anti = map.at(antipdg);
        TH1D* hNormalSig = normal.hSig;
        TH1D* hNormalBkg = normal.hBkg;
        TH1D* hAntiSig = anti.hSig;
        TH1D* hAntiBkg = anti.hBkg;

        Double_t NNormalTriggers = normal.ntriggers;
        Double_t NAntiTriggers = anti.ntriggers;

        // so now that we use strangeness conservation to determine the K0_S/L bins, how would we do this if we scale the background to account for L/Lbar asymmetry?
        // first scaling then using conservation, or vice versa?
        // Also, how do we propagate the uncertainties that come from estimating the scaling factor?

        hNormalSig->Add(hNormalBkg, -1.);
        hAntiSig->Add(hAntiBkg, -1.);

        // so, instead of this funky K0_S/L handling, we can try to use a missing energy technique: assigning all strangeness we miss to the K0_S/L
        // compare integral of hSig (after bkg subtraction) to ntriggers, and assign the difference to the K0_S/L bin. 
        // as for the statistical error, we can use the number of K0_S/L pairs we find maybe?
        int NormalMissingStrange = NNormalTriggers*normal.strangeness - hNormalSig->Integral(2, nbins); // skip the K0_S/L bin in integration
        int AntiMissingStrange = -NAntiTriggers*anti.strangeness - hAntiSig->Integral(2, nbins); // skip the K0_S/L bin in integration
        hNormalSig->SetBinContent(1, NormalMissingStrange); 
        Double_t K0Error = sqrt(NormalMissingStrange + 2*hNormalBkg->GetBinContent(1)); // calculate errors as sqrt(signal + 2*bkg)
        hNormalSig->SetBinError(1, K0Error); 
        hAntiSig->SetBinContent(1, AntiMissingStrange); 
        K0Error = sqrt(NormalMissingStrange + 2*hNormalBkg->GetBinContent(1)); // calculate errors as sqrt(signal + 2*bkg)
        hAntiSig->SetBinError(1, K0Error); 

        if(hNormalSig->GetBinContent(0) > 0 || hAntiSig->GetBinContent(0) > 0) {
            std::cout << "Warning: non-zero underflow bin detected in trigger" << hadron.getName() << std::endl;
        }

        if(hNormalSig->GetBinContent(nbins + 1) > 0 || hAntiSig->GetBinContent(nbins + 1) > 0) {
            std::cout << "Warning: non-zero overflow bin detected in trigger" << hadron.getName() << std::endl;
        }

        Double_t ratio, ratio1, errNormal, errAnti;
        // divide by zero check
        if (NNormalTriggers > 0 && NAntiTriggers > 0){
            ratio = hNormalSig->IntegralAndError(0, nbins, errNormal)/NNormalTriggers;
            ratio1 = hAntiSig->IntegralAndError(0, nbins, errAnti)/NAntiTriggers;
            errNormal /= NNormalTriggers; errAnti /= NAntiTriggers;
            std::cout << hadron.getName() << ": " << ratio << " +/- " << errNormal << " found with ntriggers = " << NNormalTriggers << std::endl;
            std::cout << hadron.getAntiName() << ": " << ratio1 << " +/- " << errAnti << " found with ntriggers = " << NAntiTriggers << std::endl;
        } else { 
            std::cout << "No triggers found for " << hadron.getName() << " and/or " << hadron.getAntiName() << std::endl;
        }
    }

    // list all the possible pdgs for both trigger and assoc:
    std::cout << "number of unique strange trigger hadrons found, should be at most " << 2*nbins << ": " << uniqueTriggerPDGs.size() << std::endl;
    std::cout << "number of unique strange assoc hadrons found, should be at most " << 2*nbins << ": " << uniqueAssocPDGs.size() << std::endl;
    // for(Int_t pdg : uniqueTriggerPDGs){std:: cout << pdg << std::endl;}
    // for(Int_t pdg : uniqueAssocPDGs){std:: cout << pdg << std::endl;}
    
    inputfile->Close();
    outputfile->Write();
    outputfile->Close();

    return;
}