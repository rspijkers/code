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

TFile *fOutput;
const int nColours = 5;
int colours[nColours] = {kRed, kBlue, kGreen, kBlack, kMagenta};

void topo(TString nrDefault, TString nrLoose, TString nrTight){
  cout << "doing topological variations..." << endl;

  TFile* fLoose = new TFile("plots/"+nrLoose+".root", "READ");
  TFile* fDefault = new TFile("plots/"+nrDefault+".root", "READ");
  TFile* fTight = new TFile("plots/"+nrTight+".root", "READ");
  TDirectory* topo = fOutput->mkdir("topo");
  topo->cd();

  for (TString s : {"OS", "SS", "Sub"}) {
    TString hname = "XiXidphiCorrected"+s+"_pT_1.0_8.0";

    TH1D *hLoose, *hDefault, *hTight;
    fLoose->GetObject("XiXi/"+hname, hLoose);
    fDefault->GetObject("XiXi/"+hname, hDefault);
    fTight->GetObject("XiXi/"+hname, hTight);

    if(!hLoose || !hDefault || !hTight){
      std::cout << "Could not retrieve all histograms for " << hname << std::endl;
      continue;
    }
    hLoose->SetName("Loose"+hname);
    hTight->SetName("Tight"+hname);

    if(s == "Sub"){ // in case of subtracted histogram, do subtraction of variations instead of ratio (ratio --> 0/0)
      SubtractAndComputeRogerBarlow(hLoose, hDefault);
      hLoose->Write();
      SubtractAndComputeRogerBarlow(hTight, hDefault);
      hTight->Write();
    } else {
      DivideAndComputeRogerBarlow(hLoose, hDefault);
      hLoose->Write();
      DivideAndComputeRogerBarlow(hTight, hDefault);
      hTight->Write();
    }
    // Do systematic uncertainty calculations here
  }
}

void SBcheck(TString run){
  cout << "doing SBcheck..." << endl;

  TFile* f = new TFile("plots/"+run+".root", "READ");
  TDirectory* dSBcheck = fOutput->mkdir("SBcheck");
  dSBcheck->cd();

  TCanvas *c1 = new TCanvas();
  c1->cd();
  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  for (TString s : {"OS", "SS", "Sub"}){
    TH1D *hNoSBcorrection, *hSBcorrection;
    f->GetObject("XiXi/XiXidphi"+s+"_pT_1.0_8.0", hNoSBcorrection);
    f->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", hSBcorrection);

    if(!hNoSBcorrection || !hSBcorrection){
      std::cout << "Could not retrieve all histograms for " << s << std::endl;
      continue;
    }
    hNoSBcorrection->SetName("NoSBcorrection"+s);
    hSBcorrection->SetName("SBcorrection"+s);

    hNoSBcorrection->SetLineColor(colours[0]);
    hSBcorrection->SetLineColor(colours[1]);
    leg->AddEntry(hNoSBcorrection, "No SB correction", "l");
    leg->AddEntry(hSBcorrection, "SB correction", "l");
    hNoSBcorrection->Draw();
    hSBcorrection->Draw("SAME");
    leg->Draw();
    c1->SetName("SBcheck_"+s);
    c1->Write();
    c1->Clear();
    leg->Clear();

    DivideAndComputeRogerBarlow(hSBcorrection, hNoSBcorrection);
    hSBcorrection->Write();
    double error;
    double avg = hSBcorrection->IntegralAndError(1, hSBcorrection->GetNbinsX(), error);
    avg /= hSBcorrection->GetNbinsX(); error /= hSBcorrection->GetNbinsX();
    std::cout << "Average ratio for " << s << ": " << avg << " +/- " << error << std::endl;
    std::cout << "Significance: " << std::abs(1.-avg)/error << " sigma" << std::endl;
  }

  fOutput->cd(); // go back to root of output file
}

void SBvariations(TString run4_10, TString run4_8, TString run6_10){
  cout << "doing SB variations..." << endl;
  TFile* f4_10 = new TFile("plots/"+run4_10+".root", "READ");
  TFile* f4_8 = new TFile("plots/"+run4_8+".root", "READ");
  TFile* f6_10 = new TFile("plots/"+run6_10+".root", "READ");
  TDirectory* dSBvariations = fOutput->mkdir("SBvariations");
  dSBvariations->cd();

  TCanvas *c1 = new TCanvas();
  c1->cd();
  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  for (TString s : {"OS", "SS", "Sub"}){
    TH1D *h4_10, *h4_8, *h6_10;
    f4_10->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", h4_10);
    f4_8->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", h4_8);
    f6_10->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", h6_10);

    if(!h4_10 || !h4_8 || !h6_10){
      std::cout << "Could not retrieve all histograms for " << s << std::endl;
      continue;
    }
    h4_10->SetName("SB4_10_"+s);
    h4_8->SetName("SB4_8_"+s);
    h6_10->SetName("SB6_10_"+s);
    h4_10->SetLineColor(colours[0]);
    h4_8->SetLineColor(colours[1]);
    h6_10->SetLineColor(colours[2]);
    leg->AddEntry(h4_10, "SB 4-10 sigma", "l");
    leg->AddEntry(h4_8, "SB 4-8 sigma", "l");
    leg->AddEntry(h6_10, "SB 6-10 sigma", "l");
    h4_10->Draw();
    h4_8->Draw("SAME");
    h6_10->Draw("SAME");
    leg->Draw();
    c1->SetName("SBvariation_"+s);
    c1->Write();
    c1->Clear();
    leg->Clear();

    if (s == "Sub") { // in case of subtracted histogram, do subtraction of variations instead of ratio (ratio --> 0/0)
      SubtractAndComputeRogerBarlow(h4_8, h4_10);
      h4_8->Write();
      SubtractAndComputeRogerBarlow(h6_10, h4_10);
      h6_10->Write();
    } else {
      DivideAndComputeRogerBarlow(h4_8, h4_10);;
      h4_8->Write();
      DivideAndComputeRogerBarlow(h6_10, h4_10);
      h6_10->Write();
    }
  }
}

void MCClosure(TString runnr){
  cout << "doing MC closure..." << endl;
  TFile* fData = new TFile("plots/"+runnr+".root", "READ");
  TFile* fMC = new TFile("plots/Closure"+runnr+".root", "READ");
  
  TDirectory* dMCClosure = fOutput->mkdir("MCClosure");
  dMCClosure->cd();

  TCanvas *c1 = new TCanvas();
  c1->cd();
  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  for (TString s : {"OS", "SS", "Sub"}){
    TH1D *hData, *hMC;
    fData->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", hData);
    fMC->GetObject("h"+s, hMC);

    if(!hData || !hMC){
      std::cout << "Could not retrieve all histograms for " << s << std::endl;
      continue;
    }

    // plot data vs MC
    hData->SetName("Data"+s);
    hMC->SetName("MC"+s);
    hData->SetLineColor(colours[0]);
    hMC->SetLineColor(colours[1]);
    leg->AddEntry(hData, "Reco", "l");
    leg->AddEntry(hMC, "Gen", "l");
    hData->Draw();
    hMC->Draw("SAME");
    leg->Draw();
    c1->SetName("MCClosure_"+s);
    c1->Write();
    c1->Clear();
    leg->Clear();

    if(s == "Sub"){ // in case of subtracted histogram, do subtraction of variations instead of ratio (ratio --> 0/0)
    // if(true){ // in case of subtracted histogram, do subtraction of variations instead of ratio (ratio --> 0/0)
      SubtractAndComputeRogerBarlow(hMC, hData);
      hMC->Write();
      checkSignificance(hMC, "difference");
    } else {
      DivideAndComputeRogerBarlow(hMC, hData);
      hMC->Write();
      checkSignificance(hMC, "ratio");
    }
  }
}

int systematics(){
  fOutput = new TFile("systematics.root", "RECREATE");
  fOutput->cd();

  // this is the run with the default settings. In case of variations in postprocessing, simply use this one. 
  TString standard = "577711";

  // // topological variations
  // TString topoLoose = "541739";
  // TString topoDefault = standard;
  // TString topoTight = "541740";

  TString topoLoose = "577712";
  TString topoDefault = "577711";
  TString topoTight = "578138";  

  topo(topoDefault, topoLoose, topoTight);

  // check whether sideband correction is significant
  // as the comparison between with and without correction uses the same dataset, the statistical uncertainty is fully correlated.
  // Hence, we use the Roger Barlow method to determine whether the difference is significant
  SBcheck(standard);

  // vary the mass window for the sideband to estimate the systematic uncertainty from the sideband subtraction
  TString SB_4_10 = "541738_SB_4_10";
  TString SB_4_8 = "541738_SB_4_8";
  TString SB_6_10 = "541738_SB_6_10";
  // SBvariations(SB_4_10, SB_4_8, SB_6_10);

  TString closureRun = "577062";
  MCClosure(closureRun);

  fOutput->Write();
  fOutput->Close();
  return 0;
} 