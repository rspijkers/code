// std
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <cmath>
#include <cassert>
// // json parsing
// #include <nlohmann/json.hpp>
// ROOT
#include "TROOT.h"
#include "TLatex.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TString.h"
#include "TStyle.h"
#include "TGrid.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRatioPlot.h"

#include "postprocessingTools.h"

std::vector<TString> plotnames = {"hV0Radius", "hCascRadius", "hV0CosPA", "hCascCosPA", "hDCAPosToPV",
   "hDCANegToPV", "hDCABachToPV", "hDCAV0ToPV", "hDCAV0Dau", "hDCACascDau", "hLambdaMass", "hITSnClustersPos", 
   "hITSnClustersNeg", "hITSnClustersBach", "hTPCnCrossedRowsPos", "hTPCnCrossedRowsNeg", "hTPCnCrossedRowsBach"};

int topocheck(){

  TFile *fDefault = new TFile("results/576504/AnalysisResults.root", "READ");
  TFile *fLoose = new TFile("results/576505/AnalysisResults.root", "READ");
  TFile *fTight = new TFile("results/576506/AnalysisResults.root", "READ");
  TFile *fDataDefault = new TFile("results/577711/AnalysisResults.root", "READ");
  TFile *fDataLoose = new TFile("results/577712/AnalysisResults.root", "READ");
  TFile *fDataTight = new TFile("results/578138/AnalysisResults.root", "READ");

  TFile *fOutput = new TFile("topocheck.root", "RECREATE");
  fOutput->cd();
  TDirectory* dSelvar = fOutput->mkdir("selvar");
  TDirectory* dLoose = fOutput->mkdir("loose");
  TDirectory* dTight = fOutput->mkdir("tight");
  TDirectory* dDefault = fOutput->mkdir("default");
  TDirectory* dDataDefault = fOutput->mkdir("data_default");
  TDirectory* dDataLoose = fOutput->mkdir("data_loose");
  TDirectory* dDataTight = fOutput->mkdir("data_tight");
  TDirectory* dBefore = fOutput->mkdir("before");

  TCanvas *c1 = new TCanvas();
  c1->cd();
  TLegend *leg = new TLegend(0.4,0.8,0.7,0.9);

  for(TString name : plotnames){
    dSelvar->cd();
    for(TFile* f : {fDataLoose, fDataDefault, fDataTight}){
      TString label;
      if(f == fDataLoose) label = "Loose";
      else if(f == fDataDefault) label = "Default";
      else if(f == fDataTight) label = "Tight";
      TH1F *h;
      f->GetObject("cascade-correlations/" + name, h);
      if(!h){
        std::cout << "Could not retrieve histogram " << name << " from file " << f->GetName() << std::endl;
        continue;
      }
      h->SetStats(kFALSE);
      h->SetLineColor(label == "Loose" ? kRed : (label == "Default" ? kBlack : kBlue));
      leg->AddEntry(h, label.Data(), "l");
      h->Draw(label == "Loose" ? "" : "SAME");
      c1->SetName(h->GetName());
    }
    leg->Draw();
    c1->Write();
    c1->Clear();
    leg->Clear();

    for(TFile* f : {fDataLoose, fDataDefault, fDataTight}){
      if(f == fDataLoose) dDataLoose->cd();
      else if(f == fDataDefault) dDataDefault->cd();
      else if(f == fDataTight) dDataTight->cd();

      TH3F *hBefore3D;
      f->GetObject("cascade-selector/" + name, hBefore3D);
      if(!hBefore3D){
        std::cout << "Could not retrieve histogram " << name << " from file " << f->GetName() << std::endl;
        continue;
      }
      TH1F* hBefore = (TH1F*)hBefore3D->Project3D("x");
      hBefore->SetTitle(name);

      TH1F* hAfter;
      f->GetObject("cascade-correlations/" + name, hAfter);
      if(!hAfter){
        std::cout << "Could not retrieve after histogram " << name << " from file " << f->GetName() << std::endl;
        continue;
      }
      // ad hoc rebinning in case of DCA to PV
      for (TString check : {"hDCAPosToPV", "hDCANegToPV", "hDCABachToPV", "hDCAV0ToPV"})
        if(name == check) hBefore->Rebin(5);

      hBefore->SetStats(kFALSE);
      hBefore->SetLineColor(kBlack);
      hAfter->SetLineColor(kRed);
      leg->AddEntry(hBefore, "before", "l");
      leg->AddEntry(hAfter, "after", "l");
      hBefore->Draw();
      hAfter->Draw("SAME");
      leg->Draw();
      c1->SetName(hBefore->GetName());
      c1->Write();
      c1->Clear();
      leg->Clear();
    }
  }

  for (TFile* f : {fLoose, fDefault, fTight}){
    TDirectory* d;
    if(f == fLoose) d = dLoose;
    else if(f == fDefault) d = dDefault;
    else if(f == fTight) d = dTight;

    d->cd();
    for(TString name : plotnames){
      TH3F *hBefore3D;
      f->GetObject("cascade-selector/" + name, hBefore3D);
      if(!hBefore3D){
        std::cout << "Could not retrieve histogram " << name << " from file " << f->GetName() << std::endl;
        continue;
      }
      TH1F* hBefore = (TH1F*)hBefore3D->Project3D("x");
      hBefore->SetTitle(name);

      TH1F* hTrueRec;
      f->GetObject("cascade-selector/truerec/" + name, hTrueRec);
      if(!hTrueRec){
        std::cout << "Could not retrieve true-rec histogram " << name << " from file " << f->GetName() << std::endl;
        continue;
      } 

      hBefore->SetStats(kFALSE);
      hBefore->SetLineColor(kBlack);
      hTrueRec->SetLineColor(kRed);
      leg->AddEntry(hBefore, "before", "l");
      leg->AddEntry(hTrueRec, "true-rec", "l");
      hBefore->Draw();
      hTrueRec->Draw("SAME");
      leg->Draw();
      c1->SetName(hBefore->GetName());
      c1->Write();
      c1->Clear();
      leg->Clear();
    }
  }

  dBefore->cd();
  for(TString name : plotnames){
    TH3F *hBefore3DMC, *hBefore3DData;
    fLoose->GetObject("cascade-selector/" + name, hBefore3DMC);
    hBefore3DMC->SetName("BeforeMC_"+name);
    fDataLoose->GetObject("cascade-selector/" + name, hBefore3DData);
    hBefore3DData->SetName("BeforeData_"+name);
    if(!hBefore3DMC || !hBefore3DData){
      std::cout << "Could not retrieve histogram " << name << " from files " << fLoose->GetName() << " and " << fDataLoose->GetName() << std::endl;
      continue;
    }
    // cout << hBefore3DData->GetZaxis()->GetXmin() << endl;
    // hBefore3DData->GetZaxis()->SetRangeUser(4,6);
    hBefore3DData->GetYaxis()->SetRangeUser(1.31, 1.33);
    // hBefore3DMC->GetZaxis()->SetRangeUser(4,6);
    hBefore3DMC->GetYaxis()->SetRangeUser(1.31, 1.33);
    TH1F* hBeforeData = (TH1F*)hBefore3DData->Project3D("x");
    TH1F* hBeforeMC = (TH1F*)hBefore3DMC->Project3D("x");
    hBeforeData->Scale(1./hBeforeData->Integral());
    hBeforeMC->Scale(1./hBeforeMC->Integral());

    hBeforeData->SetStats(kFALSE);
    hBeforeData->SetLineColor(kBlue);
    hBeforeMC->SetLineColor(kBlack);
    leg->AddEntry(hBeforeData, "data", "l");
    leg->AddEntry(hBeforeMC, "MC", "l");
    // hBeforeData->SetTitle(name);
    leg->AddEntry((TObject*)nullptr, "min pT > 1 GeV/c", "");
    hBeforeData->Draw("HIST");
    hBeforeMC->Draw("SAME HIST");
    leg->Draw();
    c1->SetName(hBeforeData->GetName());
    c1->Write();
    c1->Clear();
    leg->Clear();
  }
  // fOutput->Write();
  fOutput->Close();
  for(TFile* f : {fLoose, fDefault, fTight}){
    f->Close();
  }

  return 0;
}