#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TTree.h"
#include "TPad.h"
#include <cmath> // needed for modulo in deltaPhi calc
#include <vector> // needed for a variable array
#include <cstring> // string/char handling

void tree_handler() {
    
    TFile *inputfile = new TFile("output/treetest25k.root", "READ");
    TTree *tree = (TTree*) inputfile->Get("tree");
    // tree->Print();
    // tree->Show();

    TFile *outputfile = new TFile("output/plots.root", "RECREATE");

    Int_t Nss = tree->Draw("pdgAssoc", "pdgTrigger == 3312 && pdgAssoc > 0", "goff");
    Nss += tree->Draw("pdgAssoc", "pdgTrigger == -3312 && (pdgAssoc < 0 || pdgAssoc == 310 || pdgAssoc == 130)", "goff");
    Int_t Nos = tree->Draw("pdgAssoc", "pdgTrigger == 3312 && (pdgAssoc < 0 || pdgAssoc == 310 || pdgAssoc == 130)", "goff");
    Nos += tree->Draw("pdgAssoc", "pdgTrigger == -3312 && pdgAssoc > 0", "goff");
    Int_t Nt = tree->Draw("pdgTrigger", "pdgTrigger == 3312", "goff");
    Double_t ratio = ((Double_t) Nos - (Double_t) Nss)/((Double_t) Nt);
    // gPad->Update();
    // tree->Draw("pdgAssoc", "pdgTrigger == 3312");
    // Int_t entries = tree->GetEntries("pdgTrigger");
    std::cout << Nss << std::endl;
    std::cout << Nos << std::endl;
    std::cout << Nt << std::endl;
    std::cout << ratio << std::endl;

    inputfile->Close();
    outputfile->Close();

    // try to get the total strangeness correlated with Xi?

    return;
}