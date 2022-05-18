#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "THStack.h"
#include <vector>
#include <unordered_map>
using std::cout;
using std::endl;

#include "hadrons.h"
#include "helperfunctions.h"
#include "myStyle.h" // Load TStyle myStyle
using namespace customStyle; // custom style objects are stored in this namespace

void compare() {
    myStyle->cd(); // set custom style
    gROOT->ForceStyle(); // force custom style on objects created with a different style

    std::vector<NamedFile*> infiles;

    // TString inputdir = "output/ssbar_results_split_by_particle-antiparticle";
    // TFile *monash = new TFile(inputdir + "/plots_monash.root", "READ");
    // infiles.push_back(monash);
    // TFile *skands_mode2 = new TFile(inputdir + "/plots_skands_mode2.root", "READ");
    // infiles.push_back(skands_mode2);
    // TFile *ropes = new TFile(inputdir + "/plots_ropes.root", "READ");
    // infiles.push_back(ropes);
    // TFile *monashppbar = new TFile(inputdir + "/plots_monash_ppbar.root", "READ");
    // infiles.push_back(monashppbar);

    NamedFile *monash = (NamedFile*) new TFile("output/plots_monash_test.root", "READ");
    monash->SetCustomName("Monash");
    infiles.push_back(monash);
    NamedFile *monashppbar = (NamedFile*) new TFile("output/plots_monash_ppbar_test.root", "READ");
    monashppbar->SetCustomName("Monash p#bar{p}");
    infiles.push_back(monashppbar);

    TFile *outfile = new TFile("output/output_compare.root", "RECREATE");
    
    Int_t nfiles = infiles.size();

    TList *list = infiles[0]->GetListOfKeys();
    for(auto *key : *list){
        cout << key->GetName() << endl;
        TString name = (TString) key->GetName();
        TString ctitle = TString("title of ") + name;
        TCanvas *c = new TCanvas(name, ctitle, 800, 800);
        TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9, "PYTHIA Tunes");
        TString plottitle = key->GetTitle();
        int i = 0;
        THStack *hStack = new THStack(name, plottitle);
        for(NamedFile *file : infiles){
            TH1D *h = (TH1D*) file->Get(name);
            h->SetLineColor(colors[i]);
            h->SetOption("HIST"); // plot histogram style
            hStack->Add(h);
            legend->AddEntry(h, file->GetCustomName());
            i++;
        }
        hStack->Draw("nostack E1"); // don't stack histo's, draw errors
        legend->Draw();
        c->Write();
    }
    for(TFile *file : infiles) file->Close();
    outfile->Close();
}