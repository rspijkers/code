// std
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
// ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
// custom
#include "hadrons.h"
#include "helperfunctions.h"

const double PI = 3.14159265358979323846; 

void tree_handler() {
    TH1::SetDefaultSumw2(); // make sure errors are propagated in histo's
    
    // stbc input filepath = "/data/alice/rspijker/output/monash_pp_50M_14TeV_Trigger4GeV/"

    TChain* chain = new TChain("tree");
    chain->Add("output/ssbar_monash_ppbar_5M_14TeV.root");

    // Create list of files 
    TObjArray* fileElements = chain->GetListOfFiles(); // not actually a list of files, hence the TChainElement fuckery
    Int_t nfiles = fileElements->GetEntries();
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
        return;
    }
    
    // handle the QA histo's: 
    // from the "first" file, get the relevant QA's and clone them so we can add the others
    TH1D* hPDGTrigger = (TH1D*) fileList[0]->Get<TH1D>("hPDG")->Clone();
    TH1D* hPDGAssoc = (TH1D*) fileList[0]->Get<TH1D>("hPDGAssoc")->Clone();

    for(Int_t i = 1; i < nfiles; i++){
        TFile* file = fileList[i];
        for(TH1D* h : {hPDGTrigger, hPDGAssoc}){
            TH1D* h_ = file->Get<TH1D>(h->GetName());
            if(!h_){ // check h_ is not nullpointer
                std::cout << "file is nullptr, something went wrong" << std::endl;
                return;
            } 
            h->Add(h_);
        }
    }

    // Check to see if the under/overflow bins in the pdg are empty. If not, smt funky is going on
    Int_t nPDGbins = hPDGTrigger->GetNbinsX();
    if(hPDGTrigger->GetBinContent(0) > 0 || hPDGTrigger->GetBinContent(nPDGbins + 1) > 0){
        std::cout << "WARNING: non-empty underflow and/or overflow bins in the PDG histogram! Please investigate" << std::endl;
    }

    TFile* outputfile = new TFile("ppbar_scaled.root", "RECREATE");

    // Analysis options
    const Double_t maxEtaTrigger = 2.0;
    const Double_t maxEtaAssoc = 3.0;
    const Bool_t BkgScaling =  true; // do bkg scaling, needed to account for asymmetry in pp collisions
    std::set<Int_t> uniqueTriggerPDGs, uniqueAssocPDGs; // to keep track of possible PDG

    // Set branch address, so we can use the variables when we GetEntry()
    Double_t pTssbar;
    Int_t pdgTrigger;
    std::vector<Int_t>* pdgAssoc = 0;
    Double_t pTTrigger, etaTrigger;
    std::vector<Double_t>* pTAssoc = 0;
    std::vector<Double_t>* etaAssoc = 0;
    std::vector<Double_t>* deltaPhi = 0;
    std::vector<Double_t>* deltaEta = 0;
    chain->SetBranchAddress("pTssbar", &pTssbar);
    chain->SetBranchAddress("pdgTrigger", &pdgTrigger);
    chain->SetBranchAddress("pdgAssoc", &pdgAssoc);
    chain->SetBranchAddress("pTTrigger", &pTTrigger);
    chain->SetBranchAddress("etaTrigger", &etaTrigger);
    chain->SetBranchAddress("pTAssoc", &pTAssoc);
    chain->SetBranchAddress("etaAssoc", &etaAssoc);
    chain->SetBranchAddress("deltaPhi", &deltaPhi);
    chain->SetBranchAddress("deltaEta", &deltaEta);

    std::vector<Hadron> hadron_vec = {Kminus, Lambda, Sigmaminus, Sigmazero, Sigmaplus, Ximinus, Xizero, Omegaminus, Dsubs, Bsubs, Xicplus, Xiczero, Xibmin, Xibzero, Omegac, Omegab, Omegacc, Omegabb, Omegabc};
    const Int_t nbins = hadron_vec.size() + 1; // + 1 for K0_S/L
    
    struct mapStruct
    {
        Hadron hadron;
        TH1D* hSig; 
        TH1D* hBkg; 
        Int_t strangeness;
        Int_t ntriggers = 0;
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

    // "Ratio" plots
    TH1D* hTemp = new TH1D("template", "template", nbins, 0, nbins); 
    TH1D* hTempBar = new TH1D("templatebar", "templatebar", nbins, 0, nbins); 
    TAxis* ax = hTemp->GetXaxis();
    ax->SetBinLabel(1, "K^{0}_{S/L}");
    hTemp->SetXTitle("Associate hadron");
    TAxis* ax1 = hTempBar->GetXaxis();
    ax1->SetBinLabel(1, "K^{0}_{S/L}");
    hTempBar->SetXTitle("Associate hadron");
    for (Int_t i = 1; i < nbins; i++){
        // 0th bin is underflow, 1st bin is K0_S/L, so start with i + 2
        ax->SetBinLabel(i + 1, hadron_vec[i-1].getLatex());
        ax1->SetBinLabel(i + 1, hadron_vec[i-1].getAntiLatex());
    } 
    std::unordered_map<Int_t, mapStruct> map; 
    for (Hadron hadron : hadron_vec) {
        TString name = hadron.getName();
        TString antiname = hadron.getAntiName();

        TH1D* hss = (TH1D*) hTempBar->Clone();
        hss->SetName(name + "_bkg");
        hss->SetTitle(hadron.getLatex() + " trigger strange correlations");
        TH1D* hssbar = (TH1D*) hTempBar->Clone();
        hssbar->SetName(name + "_sig");
        hssbar->SetTitle(hadron.getLatex() + " trigger anti-strange correlations");
        TH1D* hsbars = (TH1D*) hTemp->Clone();
        hsbars->SetName(antiname + "_sig");
        hsbars->SetTitle(hadron.getAntiLatex() + " trigger strange correlations");
        TH1D* hsbarsbar = (TH1D*) hTemp->Clone();
        hsbarsbar->SetName(antiname + "_bkg");
        hsbarsbar->SetTitle(hadron.getAntiLatex() + " trigger anti-strange correlations");
        
        Int_t pdg = hadron.getPDG();
        Int_t strange = strangenessFromPDG(pdg);

        mapStruct normal; // normal pdg
        normal.hadron = hadron;
        normal.hSig = hssbar;
        normal.hBkg = hss;
        normal.strangeness = strange;
        map[pdg] = normal;

        mapStruct anti; // anti pdg
        anti.hadron = hadron;
        anti.hSig = hsbars;
        anti.hBkg = hsbarsbar;
        anti.strangeness = -1*strange;
        map[-1*pdg] = anti;
    }

    mapStruct k0shortstruct, k0longstruct;
    k0shortstruct.hadron = Kzeroshort;
    k0longstruct.hadron = Kzerolong;
    map[Kzeroshort.getPDG()] = k0shortstruct;
    map[Kzerolong.getPDG()] = k0longstruct;

    const Long64_t nentries = chain->GetEntries(); 
    Double_t pTDuplicateCheck = -1;
    for(Long64_t i = 0; i < nentries; i++) {
        chain->GetEntry(i);

        uniqueTriggerPDGs.insert(pdgTrigger); // will only insert if unique
        Hadron trigger;
        // try to access element of map. if doesn't exist, skip this entire iteration
        try{
            trigger = map.at(pdgTrigger).hadron;
        } catch (std::out_of_range){
            std::cout << "Unknown trigger hadron! pdg = " << pdgTrigger << ". Skipping this iteration..." << std::endl;
            continue;
        }
        // inclusive spectra before kine cuts
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
                if(aStrangeness > 0){
                    hSig->Fill(assoc.getLatex(), abs(aStrangeness));
                } else {
                    hSig->Fill(assoc.getAntiLatex(), abs(aStrangeness));
                }
            } else if ((tStrangeness > 0) == (aStrangeness > 0)){ // same sign
                if(aStrangeness > 0){
                    hBkg->Fill(assoc.getAntiLatex(), abs(aStrangeness));
                } else {
                    hBkg->Fill(assoc.getLatex(), abs(aStrangeness));
                }
            } else {
                std::cout << "wtf did you do???" << std::endl;
            }

            // fill deltaphi plots
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
            else if(pdgTrigger == Kminus.getPDG()){
                if(pdg == Kminus.getPDG()){hKKdphi->Fill((*deltaPhi)[j]);} // KminKmin
                else if(pdg == Kminus.getAntiPDG()){hKKbardphi->Fill((*deltaPhi)[j]);} // KminKplus
            } else if(pdgTrigger == Kminus.getAntiPDG()){
                if(pdg == Kminus.getPDG()){hKbarKdphi->Fill((*deltaPhi)[j]);} // KplusKmin
                else if(pdg == Kminus.getAntiPDG()){hKbarKbardphi->Fill((*deltaPhi)[j]);} // KplusKplus
            }
        }
    }

    // make the bkg scaling histogram so we can scale by multiplying/dividing by a histogram
    TH1D* hBkgScaling = (TH1D*) hTemp->Clone();
    if(BkgScaling){ // allow for switching the scaling on/off
        hBkgScaling->SetName("hBkgScaling");
        hBkgScaling->SetTitle("BkgScaling title");
        TAxis* PDGAxis = hPDGAssoc->GetXaxis();
        for(Hadron hadron : hadron_vec){
            Double_t Nanti = hPDGAssoc->GetBinContent(PDGAxis->FindBin(hadron.getAntiPDG()));
            Double_t Nnorm =  hPDGAssoc->GetBinContent(PDGAxis->FindBin(hadron.getPDG()));
            Double_t scale, error;
            if(Nanti == 0 || Nnorm == 0){ // if there are no normal and/or anti assocs, define scale to be 1. (i.e. no scaling)
                scale = 1;
                error = 0;
            } else {
                scale = Nanti/Nnorm;
                // calculate error as follows (see wiki of propagation of uncertainty)
                error = scale*sqrt(1/Nnorm + 1/Nanti);
                // it makes no sense to define a scaling with a large error, this will only make it more complicated
                if(error > 0.01){ // cut at 0.01 practically means that we implement the scaling only for Kaons, Lambda's, Sigma's, and Xi's, note that this may differ between productions
                    scale = 1;
                    error = 0;
                }
            }
            Int_t binnr = hBkgScaling->Fill(hadron.getLatex(), scale); // Filling this way returns the bin number which has been filled
            hBkgScaling->SetBinError(binnr, error);
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

        // Do bkg scaling first, then use strangeness conservation to calculate K0_S/L. This worsens the error in K0, but is more realistic
        if(BkgScaling){ // allow for switching the scaling on/off
            hNormalBkg->Multiply(hBkgScaling);
            hAntiBkg->Divide(hBkgScaling); 
        }

        hNormalSig->Add(hNormalBkg, -1.);
        hAntiSig->Add(hAntiBkg, -1.);

        // We can try to use a missing energy technique: assigning all strangeness we miss to the K0_S/L
        // compare integral of hSig (after bkg subtraction) to ntriggers (total strangeness), and assign the difference to the K0_S/L bin. 
        // The statistical errors are automatically calculated through IntegralAndError(), only need to manually catch and fill them
        Double_t NormK0, NormError, AntiK0, AntiError;
        NormK0 = hNormalSig->IntegralAndError(2, nbins, NormError);
        AntiK0 = hAntiSig->IntegralAndError(2, nbins, AntiError);
        hNormalSig->SetBinContent(1, NNormalTriggers*normal.strangeness - NormK0); 
        hNormalSig->SetBinError(1, NormError); 
        hAntiSig->SetBinContent(1, -NAntiTriggers*anti.strangeness - AntiK0); 
        hAntiSig->SetBinError(1, AntiError); 

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

    std::cout << "number of unique strange trigger hadrons found, should be at most " << 2*nbins << ": " << uniqueTriggerPDGs.size() << std::endl;
    std::cout << "number of unique strange assoc hadrons found, should be at most " << 2*nbins << ": " << uniqueAssocPDGs.size() << std::endl;

    // cleanup
    hTemp->Delete(); hTempBar->Delete();
    outputfile->Write();
    outputfile->Close();
}