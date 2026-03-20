#include <vector>
#include <cstdarg>

#include "postprocessingTools.h"
#include "TRatioPlot.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"

// colours to cycle through:
const int nColours = 6;
int colours[nColours] = {kRed, kBlue, kGreen, kBlack, kMagenta, kCyan};

TFile *outFile;

std::vector<TString> CascRadLabels = {"> 0.5 cm", "> 1.5 cm", "> 2.5 cm"};
std::vector<TString> DCAPosNegToPVLabels = {"> 0.01 cm", "> 0.20 cm", "> 0.30 cm"};
std::vector<TString> varsLabels = {"Loose", "Default", "Tight"};

void plotMCClosure(TString runnr) {
  TFile *fData = new TFile("plots/"+runnr+".root", "READ");
  TFile *fMC = new TFile("plots/Closure"+runnr+".root", "READ");
  
  if(!fData){
    cout << "data file is null! aborting..." << endl;
    return;
  }
  if(!fMC){
    cout << "MC file is null! aborting..." << endl;
    return;
  }

  int nplots = 3;
  // TString dataHistNames[3] = {"XiXi/XiXidphiCorrectedOS_pT_1.0_8.0", "XiXi/XiXidphiCorrectedSS_pT_1.0_8.0", "XiXi/XiXidphiCorrectedSub_pT_1.0_8.0"};
  TString dataHistNames[3] = {"hXiXiOS", "hXiXiSS", "hXiXisubtracted"};

  TString mcHistNames[3] = {"hOS", "hSS", "hSubtracted"};

  TCanvas *c = new TCanvas("c");
  c->cd();

  for (int i = 0; i < nplots; i++){
    TH1D *hData, *hMC;
    fData->GetObject(dataHistNames[i], hData);
    fMC->GetObject(mcHistNames[i], hMC);

    hData->SetLineColor(kRed);
    hData->SetStats(false);

    TRatioPlot *hRatio = new TRatioPlot(hData, hMC, "divsym");
    hRatio->SetH1DrawOpt("E");
    hRatio->Draw("nogrid");
    hRatio->GetUpperRefYaxis()->SetRangeUser(0., 0.035);
    hRatio->GetLowerRefYaxis()->SetRangeUser(0., 2.);
    TLegend *l = new TLegend(0.55, 0.75, 0.75, 0.85);
    l->AddEntry(hData, "data");
    l->AddEntry(hMC, "MC truth");
    l->Draw();
    c->Print("figures/closure/"+mcHistNames[i]+".pdf");
    c->Clear();
  }

}

int plotEfficiencies(std::vector<TString> runnumbers, TString name){
  TCanvas *c = new TCanvas("cEff"+name);
  c->cd();
  TArrayI cols = TColor::GetPalette();
  TLegend *l = new TLegend(0.55, 0.75, 0.75, 0.85);

  // TH1D *hEffDiv;
  // int i = 0;
  // for (auto runnr : runnumbers){
  //   TFile *f = new TFile("plots/MC"+runnr+".root", "READ");
  //   if(!f){
  //     cout << "file " << runnr << " is null! skipping..." << endl;
  //     continue;
  //   }
  //   TH1D *hEffXiMinus;
  //   f->GetObject("hPtXiMinEff", hEffXiMinus);
  //   hEffXiMinus->SetStats(false);
  //   hEffXiMinus->SetTitle("Xi- Efficiency");
  //   hEffXiMinus->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //   hEffXiMinus->GetYaxis()->SetTitle("Efficiency");
  //   hEffXiMinus->SetLineColor(colours[i]);
  //   // hEffXiMinus->SetMarkerColor(kRed + runnr.Atoi()%10);
  //   // hEffXiMinus->SetMarkerStyle(20 + runnr.Atoi()%5);
  //   if (runnr == runnumbers[0]){
  //     hEffXiMinus->Draw("E");
  //     hEffDiv = (TH1D*) hEffXiMinus->Clone("hEffDiv");
  //   } else {
  //     hEffXiMinus->Draw("Esame");
  //     hEffDiv->Divide(hEffXiMinus);
  //   }
  //   l->AddEntry(hEffXiMinus, "Run "+runnr);
  //   i++;
  // }
  int n = runnumbers.size();
  for(int i = 0; i < n; i++){
    TString runnr = runnumbers[i];
    TFile *f = new TFile("plots/MC"+runnr+".root", "READ");
    if(!f){
      cout << "file " << runnr << " is null! skipping..." << endl;
      continue;
    }
    TH1D *hEff;
    f->GetObject("hPtXiMinEff", hEff);
    hEff->SetStats(false);
    hEff->SetTitle("Xi- Efficiency");
    hEff->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hEff->GetYaxis()->SetTitle("Efficiency");
    hEff->SetLineColor(colours[i % nColours]);
    if (i == 0){
      hEff->Draw("E");
    } else {
      hEff->Draw("Esame");
    }
    l->AddEntry(hEff, "Run "+runnr);
  }
  outFile->cd();
  l->Draw();
  c->Write();
  c->Clear();
  l->Clear();

  return 0;
}

int plotSpectra(std::vector<TString> runnumbers, std::vector<TString> labels, TString name){
  TCanvas *c = new TCanvas("cSpectra"+name);
  c->cd();
  TArrayI cols = TColor::GetPalette();
  TLegend *l = new TLegend(0.55, 0.75, 0.75, 0.85);

  int n = runnumbers.size();
  for(int i = 0; i < n; i++){
    TString runnr = runnumbers[i];
    TFile *f = new TFile("spectra"+runnr+".root", "READ");
    if(!f){
      cout << "file " << runnr << " is null! skipping..." << endl;
      continue;
    }
    TH1D *hSpectra;
    f->GetObject("hSpectraRebin", hSpectra);
    hSpectra->SetStats(false);
    hSpectra->SetTitle("Xi- Spectra");
    hSpectra->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hSpectra->GetYaxis()->SetTitle("1/N_{evt} dN/dp_{T} (GeV/c)^{-1}");
    hSpectra->SetLineColor(colours[i % nColours]);
    if (i == 0){
      hSpectra->Draw("E");
    } else {
      hSpectra->Draw("Esame");
    }
    l->AddEntry(hSpectra, labels[i]);
  }
  outFile->cd();
  l->Draw();
  c->Write();
  c->Clear();
  l->Clear();

  TH1D* def;
  TFile *f = new TFile("spectra"+runnumbers[1]+".root", "READ");
  f->GetObject("hSpectraRebin", def);
  for (int i = 0; i < n; i++){
    if (i == 1) continue;
    TString runnr = runnumbers[i];
    TFile *f = new TFile("spectra"+runnr+".root", "READ");
    if(!f){
      cout << "file " << runnr << " is null! skipping..." << endl;
      continue;
    }
    TH1D *hSpectra;
    f->GetObject("hSpectraRebin", hSpectra);
    hSpectra->Divide(def);
    hSpectra->SetStats(false);
    hSpectra->SetTitle("Xi- Spectra Ratio to Run "+runnumbers[1]);
    hSpectra->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hSpectra->GetYaxis()->SetTitle("Ratio");
    hSpectra->SetAxisRange(0.9, 1.1, "Y");
    hSpectra->SetLineColor(colours[i % nColours]);
    if (i == 0){
      hSpectra->Draw("E");
    } else {
      hSpectra->Draw("Esame");
    }
    l->AddEntry(hSpectra, labels[i]);
  }
  outFile->cd();
  l->Draw();
  c->Write("ratio"+name);
  c->Clear(); 
  l->Clear();
  return 0;
}

int plottingMacro(){

  // plotMCClosure(runnr);
  outFile = new TFile("plots.root", "RECREATE");
  outFile->cd();

  // plotEfficiencies({"586917", "586918", "586919", "586920", "586921", "588051"});

  // start with loose variation for automatic y-axis
  plotEfficiencies({"607020", "607019", "607021"}, "spectra"); // 24f4d - spectra
  plotEfficiencies({"607164", "607163", "607165"}, "correlations"); // 24j7 - correlations

  plotSpectra({"603231", "603230", "603229"}, CascRadLabels, "CascRad");
  plotSpectra({"603228", "603227", "603226"}, DCAPosNegToPVLabels, "DCAPosNegToPV");
  plotSpectra({"608045", "608044", "608046"}, varsLabels, "sysvars");


	return 0;
}
