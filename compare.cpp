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
#include "TKey.h"
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
    NamedFile *pp_scaling = new NamedFile("output/scalingtest/pp_scaling.root", "pp_scaling", "READ");
    infiles.push_back(pp_scaling);
    NamedFile *ppbar_scaling = new NamedFile("output/scalingtest/ppbar_scaling.root", "ppbar_scaling", "READ");
    infiles.push_back(ppbar_scaling);
    // NamedFile *pp_noscaling = new NamedFile("output/scalingtest/pp_noscaling.root", "pp_noscaling", "READ");
    // infiles.push_back(pp_noscaling);
    // NamedFile *ppbar_noscaling = new NamedFile("output/scalingtest/ppbar_noscaling.root", "ppbar_noscaling", "READ");
    // infiles.push_back(ppbar_noscaling);

    TFile *outfile = new TFile("output/output_compare.root", "RECREATE");
    
    for(auto *key : *infiles[0]->GetListOfKeys()){
        cout << key->GetName() << endl;
        TString name = (TString) key->GetName();
        TString ctitle = TString("title of "+ name);
        TCanvas *c = new TCanvas(name, ctitle, 800, 600);
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9, "PYTHIA Tunes");
        TString plottitle = key->GetTitle();
        int i = 0;
        THStack *hStack = new THStack(name, plottitle);
        for(NamedFile *file : infiles){
            TH1D *h = file->Get<TH1D>(name); // TODO: implement safeguard against trying to get something that doesnt exist
            h->SetLineColor(colors[i]);
            // h->SetOption("HIST"); // plot histogram style
            hStack->Add(h);
            legend->AddEntry(h, file->GetCustomName());
            i++;
        }
        hStack->Draw("nostack E1"); // don't stack histo's, draw errors
        hStack->GetXaxis()->SetLabelSize(0.06); // has to be here, root is dumb
        hStack->GetXaxis(); // has to be here, root is dumb
        legend->Draw();
        c->Write();
    }

    // generalize the following code to include these combi's: Lambda-Kaon, Lambda-Sigma, Kaon-Kaon (Lambda-Lambda)
    TString correlations[4][2] = {{"Lambda", "Lambda"}, {"Lambda", "Kaon"}, {"Lambda", "Sigma"}, {"Kaon", "Kaon"}};
    for (auto corr : correlations){
        // close your eyes, this is ugly
        TString trigger = TString(corr[0][0]); TString assoc = TString(corr[1][0]);
        TString signame = TString("h"+trigger+assoc+"bardphi");
        TString sigbarname = TString("h"+trigger+"bar"+assoc+"dphi");
        TString bkgname = TString("h"+trigger+assoc+"dphi");
        TString bkgbarname = TString("h"+trigger+"bar"+assoc+"bardphi");
        // okay you can open your eyes again
        for(auto production : infiles){
            TString prodname = production->GetCustomName();
            // nevermind the string fuckery continues
            TCanvas *c = new TCanvas(corr[0]+corr[1]+"bar_"+prodname, corr[0]+" - "+corr[1]+"(bar) correlations in "+prodname, 800, 800);
            TCanvas *cbar = new TCanvas(corr[0]+"bar"+corr[1]+"_"+prodname, corr[0]+"bar - "+corr[1]+" correlations in "+prodname, 800, 800);
            THStack *stack = new THStack("stack", corr[0]+" - "+corr[1]+"(bar) correlations in "+prodname);
            THStack *stackbar = new THStack("stackbar", corr[0]+"bar - "+corr[1]+" correlations in "+prodname);
            TH1D* hbkg = (TH1D*) production->Get(bkgname); hbkg->SetLineColor(kBlack);
            TH1D* hsig = (TH1D*) production->Get(signame); hsig->SetLineColor(kRed+1);
            TH1D* hsigbar = (TH1D*) production->Get(sigbarname); hsigbar->SetLineColor(kRed+1);
            TH1D* hbkgbar = (TH1D*) production->Get(bkgbarname); hbkgbar->SetLineColor(kBlack);

            stack->Add(hsig); stack->Add(hbkg);
            stackbar->Add(hsigbar); stackbar->Add(hbkgbar);
            c->cd();
            stack->Draw("nostack E1");
            c->Write();
            cbar->cd();
            stackbar->Draw("nostack E1");
            cbar->Write();
        }
    }

    for(TFile *file : infiles) file->Close();
    outfile->Close();
}