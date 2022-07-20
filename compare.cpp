
// std
#include <iostream>
#include <vector>
#include <unordered_map>
// ROOT
#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TString.h"
#include "TTree.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TKey.h"
// custom
#include "hadrons.h"
#include "helperfunctions.h"
#include "myStyle.h" 

using std::cout;
using std::endl;
using namespace customStyle; // custom style objects are stored in this namespace

void compare() {
    TH1::SetDefaultSumw2(); // make sure errors are propagated in histo's
    myStyle->cd(); // set custom style
    gROOT->ForceStyle(); // force custom style on objects created with a different style

    std::vector<NamedFile*> infiles;
    NamedFile *pp_scaled = new NamedFile("output/scalingtest/pp_scaled.root", "pp_scaled", "READ");
    infiles.push_back(pp_scaled);
    NamedFile *ppbar_scaled = new NamedFile("output/scalingtest/ppbar_scaled.root", "ppbar_scaled", "READ");
    infiles.push_back(ppbar_scaled);
    NamedFile *pp_unscaled = new NamedFile("output/scalingtest/pp_unscaled.root", "pp_unscaled", "READ");
    infiles.push_back(pp_unscaled);
    NamedFile *ppbar_unscaled = new NamedFile("output/scalingtest/ppbar_unscaled.root", "ppbar_unscaled", "READ");
    infiles.push_back(ppbar_unscaled);
    Int_t nfiles = infiles.size();

    // simply loop over all the histo's and make comparisons between different files(= productions)
    TFile *outfile = new TFile("output/output_compare.root", "RECREATE");
    for(auto *key : *infiles[0]->GetListOfKeys()){
        TString name = (TString) key->GetName();
        cout << name << endl;
        TString title = (TString) key->GetTitle();
        TCanvas *c = new TCanvas(name, title, 800, 600);
        TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9, "PYTHIA Tunes");

        TH1D* h0 = infiles[0]->Get<TH1D>(name); // get first histogram, we can use this to handle ranges and stuff
        if(!h0) { // check nullptr
            cout << "Warning: histogram is nullptr, skipping!" << endl;
            continue;
        }
        h0->SetLineColor(colors[0]);
        h0->Draw();
        legend->AddEntry(h0, infiles[0]->GetCustomName());
        for(Int_t i = 1; i < nfiles; i++){
            TH1D *h = infiles[i]->Get<TH1D>(name); 
            if(h == h0) continue; // skip h0, we already have it
            if(!h) { // check nullptr
                cout << "Warning: histogram is nullptr, skipping!" << endl;
                continue;
            }
            updateRanges(h0, h); // make sure the ranges include everything from all plots
            h->SetLineColor(colors[i]);
            h->Draw("SAME");
            legend->AddEntry(h, infiles[i]->GetCustomName());
        }
        // Fix labels
        TAxis* XAxis = h0->GetXaxis();
        for(int j = 0; j < XAxis->GetNbins(); j++) XAxis->ChangeLabel(j, -1, 0.06, 11, -1, -1, "");

        legend->Draw();
        c->Write();
    }

    // deltaPhi signal vs background plots
    TString correlations[4][2] = {{"Lambda", "Lambda"}, {"Lambda", "Kaon"}, {"Lambda", "Sigma"}, {"Kaon", "Kaon"}};
    for (auto corr : correlations){
        // Sorry for all the string fuckery
        TString trigger = TString(corr[0][0]); TString assoc = TString(corr[1][0]);
        TString signame = TString("h"+trigger+assoc+"bardphi");
        TString sigbarname = TString("h"+trigger+"bar"+assoc+"dphi");
        TString bkgname = TString("h"+trigger+assoc+"dphi");
        TString bkgbarname = TString("h"+trigger+"bar"+assoc+"bardphi");

        for(auto production : infiles){
            TString prodname = production->GetCustomName();
            TCanvas *c = new TCanvas(corr[0]+corr[1]+"bar_"+prodname, corr[0]+" - "+corr[1]+"(bar) correlations in "+prodname, 800, 800);
            TCanvas *cbar = new TCanvas(corr[0]+"bar"+corr[1]+"_"+prodname, corr[0]+"bar - "+corr[1]+"(bar) correlations in "+prodname, 800, 800);
            TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
            TLegend *legendbar = new TLegend(0.7, 0.7, 0.9, 0.9);

            TH1D* hsig = production->Get<TH1D>(signame); hsig->SetLineColor(kRed+1); legend->AddEntry(hsig, corr[0]+"-"+corr[1]+"bar");
            TH1D* hbkg = production->Get<TH1D>(bkgname); hbkg->SetLineColor(kBlack); legend->AddEntry(hbkg, corr[0]+"-"+corr[1]);
            TH1D* hsigbar = production->Get<TH1D>(sigbarname); hsigbar->SetLineColor(kRed+1); legendbar->AddEntry(hsigbar, corr[0]+"bar-"+corr[1]);
            TH1D* hbkgbar = production->Get<TH1D>(bkgbarname); hbkgbar->SetLineColor(kBlack); legendbar->AddEntry(hbkgbar, corr[0]+"bar-"+corr[1]+"bar");

            if (prodname == "pp_scaled" || prodname == "ppbar_scaled"){
                // do scaling
                TH1D* hBkgScaling = production->Get<TH1D>("hBkgScaling");
                // bin nrs are: K 2, L 3, S 4.
                Double_t scaler, error;
                if(corr[1] == "Kaon"){
                    scaler = hBkgScaling->GetBinContent(2); error = hBkgScaling->GetBinError(2);
                } else if(corr[1] == "Lambda"){
                    scaler = hBkgScaling->GetBinContent(3); error = hBkgScaling->GetBinError(3);
                } else if(corr[1] == "Sigma"){
                    scaler = hBkgScaling->GetBinContent(4); error = hBkgScaling->GetBinError(4);
                }
                TH1D* hDPhiScaler = (TH1D*) hbkg->Clone();
                for(Int_t i = 0; i < hDPhiScaler->GetNbinsX(); i++){
                    hDPhiScaler->SetBinContent(i + 1, scaler); 
                    hDPhiScaler->SetBinError(i + 1, error);
                }
                hbkg->Multiply(hDPhiScaler);
                hbkgbar->Divide(hDPhiScaler);
            }
            updateRanges(hsig, hbkg);
            updateRanges(hsigbar, hbkgbar);

            c->cd();
            hsig->SetTitle(corr[0]+" Trigger, "+prodname);
            hsig->SetXTitle("#Delta#varphi");
            hsig->Draw();
            hbkg->Draw("SAME");
            legend->Draw();
            c->Write();

            cbar->cd();
            hsigbar->Draw();
            hsigbar->SetTitle(corr[0]+"bar Trigger, "+prodname);
            hsigbar->SetXTitle("#Delta#varphi");
            hbkgbar->Draw("SAME");
            legendbar->Draw();
            cbar->Write();
            
            // Okay, let's do some calculations. We want to calculate the integral of the diff between sig and bkg from 1/2 PI to 3/2 Pi.
            // This should be compatible with 0, but ofcourse it wont be in the unscaled pp case.
            TH1D *hDiff = (TH1D*) hsig->Clone();
            hDiff->SetName("Diff"+prodname+corr[0]+corr[1]);
            TH1D *hDiffbar = (TH1D*) hsigbar->Clone();
            hDiffbar->SetName("Diff"+prodname+corr[0]+"bar"+corr[1]);

            hDiff->Add(hbkg, -1);
            hDiffbar->Add(hbkgbar, -1);
            Int_t nbins = hDiff->GetNbinsX();
            Double_t diff, diffbar, err, errbar;
            diff = hDiff->IntegralAndError(nbins/2, nbins, err);
            diffbar = hDiffbar->IntegralAndError(nbins/2, nbins, errbar);
            cout << prodname << ", " << corr[0] << " Trigger: " << diff << " +/- " << err << endl;
            cout << prodname << ", " << corr[0] << "bar Trigger: " << diffbar << " +/- " << errbar << endl;
            hDiff->Write(); hDiffbar->Write();
            // We should also make another histogram + calculation, namely the (sig+sigbar - (bkg+bkgbar))/2 one. 
            // A lot of things will average out then, but will the unscaled pp case be as good as the scaled one?
            // Quantify the difference. 

        }
    }

    

    for(TFile *file : infiles) file->Close();
    outfile->Close();
}