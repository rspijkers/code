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

// systematics: we want to have 20 bins of systematics. let's have each function output 3 histograms (OS SS Sub) so we can plot them together, and calculate the final systemaic at the end
const int bins = 20;
const bool makePDFs = true;  // Set to true to save PDF plots to figures/systematics/
const TString xytitle = "#Delta#varphi (rad);1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})";

TFile *fOutput;
TDirectory* spectraDir;
const int nColours = 5;
int colours[nColours] = {kRed, kBlue, kGreen, kBlack, kMagenta};

std::vector<TString> CascRadLabels = {"> 0.5 cm", "> 1.5 cm", "> 2.5 cm"};
std::vector<TString> DCAPosNegToPVLabels = {"> 0.01 cm", "> 0.20 cm", "> 0.30 cm"};
std::vector<TString> varsLabels = {"Loose", "Default", "Tight"};
std::vector<TString> signTypes = {"OS", "SS", "Sub"};

std::array<TH1D*, 3> topo(TString nrDefault, TString nrLoose, TString nrTight){
  cout << "--- Doing topological variations... ---" << endl;

  TFile* fLoose = new TFile("plots/"+nrLoose+".root", "READ");
  TFile* fDefault = new TFile("plots/"+nrDefault+".root", "READ");
  TFile* fTight = new TFile("plots/"+nrTight+".root", "READ");
  TDirectory* topo = fOutput->mkdir("topo");
  topo->cd();

  std::array<TH1D*, 3> res;

  for(int i = 0; i < signTypes.size(); i++){
    TString s = signTypes[i];
    TString hname = "XiXidphiCorrected"+s+"_pT_1.0_8.0";

    TH1D *hLoose, *hDefault, *hTight;
    fLoose->GetObject("XiXi/"+hname, hLoose); fDefault->GetObject("XiXi/"+hname, hDefault); fTight->GetObject("XiXi/"+hname, hTight);

    if(!hLoose || !hDefault || !hTight){
      std::cout << "Could not retrieve all histograms for " << hname << std::endl;
      continue;
    }
    hLoose->SetName("Loose"+hname); hTight->SetName("Tight"+hname);

    bool doSys = false;
    cout << "--- Processing " << s << " ---" << endl;
    
    // Store results for both variations to compute final systematic at the end
    struct VariationResults {
      double avgArithmetic, errArithmetic;
      double avgWeighted, errWeighted;
      double avgFromRatio, errFromRatio;
      bool significant;
    };
    std::map<TString, VariationResults> varResults;
    
    // Process both Loose and Tight variations
    std::vector<TString> variations = {"Loose", "Tight"};
    std::map<TString, TH1D*> varHistos = {{"Loose", hLoose}, {"Tight", hTight}};
    
    for(const auto& var : variations) {
      TH1D* hVar = varHistos[var];
      
      // Create canvas for comparing different averaging methods
      TCanvas *cComparison = new TCanvas(("cAverageComparison_"+var+"_"+s).Data(), 
                                         ("Average Comparison "+var+" "+s).Data(), 800, 800);
      cComparison->cd();
      
      // Compute averages BEFORE the RogerBarlow functions modify the histograms
      TH1D *hRatio;
      if(s == "Sub"){ // in case of subtracted histogram, do subtraction
        hRatio = SubtractAndComputeRogerBarlow(hVar, hDefault);
      } else { // ratio case
        hRatio = DivideAndComputeRogerBarlow(hVar, hDefault);
      }
      hRatio->SetName(("hRatio_"+var+"_"+s).Data());
      
      // Compute averages using different methods
      VariationResults& vr = varResults[var];
      vr.avgArithmetic = computeArithmeticMean(hRatio, vr.errArithmetic);
      vr.avgWeighted = computeWeightedAverage(hRatio, vr.errWeighted);
      if(s == "Sub") {
        vr.avgFromRatio = computeAverageDifference(hVar, hDefault, vr.errFromRatio);
      } else {
        vr.avgFromRatio = computeAverageRatio(hVar, hDefault, vr.errFromRatio);
      }
      
      cout << var << " variation:" << endl;
      cout << "  Arithmetic mean: " << vr.avgArithmetic << " +/- " << vr.errArithmetic << endl;
      cout << "  Weighted average: " << vr.avgWeighted << " +/- " << vr.errWeighted << endl;
      cout << "  Average " << (s == "Sub" ? "difference" : "ratio") << ": " << vr.avgFromRatio << " +/- " << vr.errFromRatio << endl;
      
      // Draw the ratio histogram first
      hRatio->SetLineColor(kBlack); hRatio->SetMarkerStyle(20); hRatio->SetMarkerColor(kBlack);
      hRatio->SetTitle((s+" Topological Variations ("+var+");#Delta#varphi (rad);"+TString(s == "Sub" ? "Difference" : "Ratio")).Data());
      
      // Calculate y-axis range to accommodate legend
      double yMin = 1e10, yMax = -1e10;
      for (int i = 1; i <= hRatio->GetNbinsX(); i++) {
        double val = hRatio->GetBinContent(i);
        double err = hRatio->GetBinError(i);
        if (val != 0) {
          yMin = TMath::Min(yMin, val - err);
          yMax = TMath::Max(yMax, val + err);
        }
      }
      // Include the horizontal line values
      yMin = TMath::Min(yMin, TMath::Min(vr.avgArithmetic - vr.errArithmetic, TMath::Min(vr.avgWeighted - vr.errWeighted, vr.avgFromRatio - vr.errFromRatio)));
      yMax = TMath::Max(yMax, TMath::Max(vr.avgArithmetic + vr.errArithmetic, TMath::Max(vr.avgWeighted + vr.errWeighted, vr.avgFromRatio + vr.errFromRatio)));
      
      // Set range with extra padding at top for legend
      if (yMax > yMin) {
        double range = yMax - yMin;
        hRatio->GetYaxis()->SetRangeUser(yMin - 0.05 * range, yMax + 0.35 * range);
      }
      
      hRatio->Draw("E");
      
      // Draw horizontal lines with error bands for each method
      double xMin = hRatio->GetXaxis()->GetXmin();
      double xMax = hRatio->GetXaxis()->GetXmax();
      
      // Arithmetic mean (red)
      TGraphErrors *gArithmetic = new TGraphErrors(2);
      gArithmetic->SetPoint(0, xMin, vr.avgArithmetic); gArithmetic->SetPoint(1, xMax, vr.avgArithmetic);
      gArithmetic->SetPointError(0, 0, vr.errArithmetic); gArithmetic->SetPointError(1, 0, vr.errArithmetic);
      gArithmetic->SetLineColor(kRed); gArithmetic->SetLineWidth(2); gArithmetic->SetFillColorAlpha(kRed, 0.2);
      gArithmetic->Draw("E3 L SAME");
      
      // Weighted average (blue)
      TGraphErrors *gWeighted = new TGraphErrors(2);
      gWeighted->SetPoint(0, xMin, vr.avgWeighted); gWeighted->SetPoint(1, xMax, vr.avgWeighted);
      gWeighted->SetPointError(0, 0, vr.errWeighted); gWeighted->SetPointError(1, 0, vr.errWeighted);
      gWeighted->SetLineColor(kBlue); gWeighted->SetLineWidth(2); gWeighted->SetFillColorAlpha(kBlue, 0.2);
      gWeighted->Draw("E3 L SAME");
      
      // Average ratio/difference (green)
      TGraphErrors *gRatio = new TGraphErrors(2);
      gRatio->SetPoint(0, xMin, vr.avgFromRatio); gRatio->SetPoint(1, xMax, vr.avgFromRatio);
      gRatio->SetPointError(0, 0, vr.errFromRatio); gRatio->SetPointError(1, 0, vr.errFromRatio);
      gRatio->SetLineColor(kGreen+2); gRatio->SetLineWidth(2); gRatio->SetFillColorAlpha(kGreen+2, 0.2);
      gRatio->Draw("E3 L SAME");
      
      // Add legend (top-right for OS/Sub, top-left for SS)
      bool topLeft = (s == "SS");
      TLegend *legComp = new TLegend(0.18, 0.72, 0.82, 0.87);
      legComp->AddEntry(hRatio, var+"/Default", "lep");
      if(s == "Sub") {
        legComp->AddEntry(gArithmetic, Form("Arithmetic: %.4g #pm %.4g", vr.avgArithmetic, vr.errArithmetic), "l");
        legComp->AddEntry(gWeighted, Form("Weighted: %.4g #pm %.4g", vr.avgWeighted, vr.errWeighted), "l");
        legComp->AddEntry(gRatio, Form("Diff of averages: %.4g #pm %.4g", vr.avgFromRatio, vr.errFromRatio), "l");
      } else {
        legComp->AddEntry(gArithmetic, Form("Arithmetic: %.2f%% #pm %.2f%%", std::abs((vr.avgArithmetic - 1.0) * 100.0), vr.errArithmetic * 100.0), "l");
        legComp->AddEntry(gWeighted, Form("Weighted: %.2f%% #pm %.2f%%", std::abs((vr.avgWeighted - 1.0) * 100.0), vr.errWeighted * 100.0), "l");
        legComp->AddEntry(gRatio, Form("Ratio of sums: %.2f%% #pm %.2f%%", std::abs((vr.avgFromRatio - 1.0) * 100.0), vr.errFromRatio * 100.0), "l");
      }
      legComp->Draw();
      
      cComparison->Write();
      if(makePDFs) {
        cComparison->SetTitle("");
        cComparison->Print("figures/systematics/topo_" + var + "_" + s + ".pdf");
      }
      
      // Check significance for this variation
      vr.significant = checkSignificance(vr.avgFromRatio, vr.errFromRatio, (s == "Sub" ? "difference" : "ratio"));
      
      // Create ratio plot using plotWithRatio function with proper uncertainties
      // Set colors and styles
      if(var == "Loose") {
        hVar->SetLineColor(kBlue); hVar->SetMarkerColor(kBlue); hVar->SetMarkerStyle(20);
      } else {
        hVar->SetLineColor(kRed); hVar->SetMarkerColor(kRed); hVar->SetMarkerStyle(22);
      }
      hDefault->SetLineColor(kBlack); hDefault->SetMarkerColor(kBlack); hDefault->SetMarkerStyle(21);
      
      hDefault->SetTitle((";"+xytitle).Data());
      hVar->SetTitle((";"+xytitle).Data());
      hRatio->GetYaxis()->SetTitle(s == "Sub" ? "Diff." : ("#frac{" + var + "}{Default}").Data());
      
      // Create legend (top-right for OS/Sub, top-left for SS)
      TLegend *legRatioPlot = topLeft ? new TLegend(0.17, 0.72, 0.47, 0.85) : new TLegend(0.55, 0.72, 0.85, 0.85);
      legRatioPlot->AddEntry(hVar, var, "lep");
      legRatioPlot->AddEntry(hDefault, "Default", "lep");
      
      // Create ratio plot
      TCanvas *cRatioPlot = plotWithRatio(("cRatioPlot_"+var+"_"+s).Data(), hDefault, hVar, hRatio, 
                                          "E", "E", "E", legRatioPlot, 800, 800, nullptr, nullptr, s == "Sub");
      cRatioPlot->Write();
      if(makePDFs) {
        hVar->SetTitle("");
        hDefault->SetTitle("");
        cRatioPlot->SaveAs("figures/systematics/topo_RatioPlot_" + var + "_" + s + ".pdf");
      }
    }
    
    // Determine significance based on averageDiff/Ratio for both Tight and Loose
    // Set systematic as the average of the absolute values of the two
    double avgBoth;
    if(s == "Sub") {
      avgBoth = (std::abs(varResults["Tight"].avgFromRatio) + std::abs(varResults["Loose"].avgFromRatio)) / 2.0;
    } else {
      avgBoth = (std::abs(varResults["Tight"].avgFromRatio - 1.0) + std::abs(varResults["Loose"].avgFromRatio - 1.0)) / 2.0;
    }
    doSys = varResults["Tight"].significant || varResults["Loose"].significant;

    TH1D* hSys = (TH1D*) hDefault->Clone();
    hSys->SetName((s+"TopoSys").Data());
    if(doSys){
      cout << "avg (Tight): " << varResults["Tight"].avgFromRatio << " +/- " << varResults["Tight"].errFromRatio << endl;
      cout << "avg (Loose): " << varResults["Loose"].avgFromRatio << " +/- " << varResults["Loose"].errFromRatio << endl;
      cout << "average of both: " << avgBoth << endl;
      for(int b = 1; b < hSys->GetNbinsX()+1; b++){
        hSys->SetBinContent(b, avgBoth);
      }
    } else {
      for(int b = 1; b < hSys->GetNbinsX()+1; b++){
        hSys->SetBinContent(b, 0.0);
      }
    }
    res[i] = hSys;
    hSys->Write();
  }
  return res;
}

void SBcheck(TString run){
  cout << "doing SBcheck..." << endl;

  TFile* f = new TFile("plots/"+run+".root", "READ");
  TDirectory* dSBcheck = fOutput->mkdir("SBcheck");
  dSBcheck->cd();

  TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.95);
  for(int i = 0; i < signTypes.size(); i++){
    TString s = signTypes[i];
    TH1D *hNoSBcorrection, *hSBcorrection;
    f->GetObject("XiXi/XiXidphi"+s+"_pT_1.0_8.0", hNoSBcorrection);
    f->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", hSBcorrection);

    if(!hNoSBcorrection || !hSBcorrection){
      std::cout << "Could not retrieve all histograms for " << s << std::endl;
      continue;
    }
    hNoSBcorrection->SetName("NoSBcorrection"+s); hSBcorrection->SetName("SBcorrection"+s);
    applyHistStyle(hNoSBcorrection, 0);  // First color
    applyHistStyle(hSBcorrection, 1);   // Second color
    hNoSBcorrection->SetTitle((";"+xytitle).Data());
    leg->AddEntry(hNoSBcorrection, "No SB correction", "l"); leg->AddEntry(hSBcorrection, "SB correction", "l");
    
    // Draw using drawAndSavePdf
    std::vector<TH1*> hists = {hNoSBcorrection, hSBcorrection};
    TCanvas *cSB = drawAndSavePdf(hists, makePDFs ? "figures/systematics/SBcheck_" + s + ".pdf" : "", false, 800, 800, "E", leg);
    cSB->SetName("SBcheck_"+s);
    dSBcheck->cd();
    cSB->Write();
    if(!makePDFs) delete cSB;

  TH1D *hRatio = DivideAndComputeRogerBarlow(hSBcorrection, hNoSBcorrection);
  hRatio->Write();
  double error;
  double avg = hRatio->IntegralAndError(1, hRatio->GetNbinsX(), error);
    avg /= hSBcorrection->GetNbinsX(); error /= hSBcorrection->GetNbinsX();
    std::cout << "Average ratio for " << s << ": " << avg << " +/- " << error << std::endl;
    std::cout << "Significance: " << std::abs(1.-avg)/error << " sigma" << std::endl;
  }
  fOutput->cd(); // go back to root of output file
}

std::array<TH1D*, 3> SBvariations(TString run4_10, TString run4_8, TString run6_10){
  cout << "doing SB variations..." << endl;
  TFile* f4_10 = new TFile("plots/"+run4_10+".root", "READ");
  TFile* f4_8 = new TFile("plots/"+run4_8+".root", "READ");
  TFile* f6_10 = new TFile("plots/"+run6_10+".root", "READ");
  TDirectory* dSBvariations = fOutput->mkdir("SBvariations");
  dSBvariations->cd();

  std::array<TH1D*, 3> res;

  TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.95);
  for(int i = 0; i < signTypes.size(); i++){
    TString s = signTypes[i];
    TH1D *h4_10, *h4_8, *h6_10;
    f4_10->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", h4_10); f4_8->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", h4_8); f6_10->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", h6_10);

    if(!h4_10 || !h4_8 || !h6_10){
      std::cout << "Could not retrieve all histograms for " << s << std::endl;
      continue;
    }
    h4_10->SetName("SB4_10_"+s); h4_8->SetName("SB4_8_"+s); h6_10->SetName("SB6_10_"+s);
    applyHistStyle(h4_10, 0);
    applyHistStyle(h4_8, 1);
    applyHistStyle(h6_10, 2);
    h4_10->SetTitle((";"+xytitle).Data());
    leg->AddEntry(h4_10, "SB 4-10 sigma", "l"); leg->AddEntry(h4_8, "SB 4-8 sigma", "l"); leg->AddEntry(h6_10, "SB 6-10 sigma", "l");
    
    // Draw using drawAndSavePdf
    std::vector<TH1*> histsSB = {h4_10, h4_8, h6_10};
    TCanvas *cSB = drawAndSavePdf(histsSB, makePDFs ? "figures/systematics/SBvariations_" + s + ".pdf" : "", false, 800, 800, "E", leg);
    cSB->SetName("SBvariation_"+s);
    dSBvariations->cd();
    cSB->Write();
    if(!makePDFs) delete cSB;
    leg->Clear();

    bool doSys = false;
    double avg, error;
    TH1D *hRB1, *hRB2;
    if (s == "Sub") { // in case of subtracted histogram, do subtraction of variations instead of ratio (ratio --> 0/0)
      hRB1 = SubtractAndComputeRogerBarlow(h4_8, h4_10);
      hRB2 = SubtractAndComputeRogerBarlow(h6_10, h4_10);
      avg = computeArithmeticMean(hRB1, error);
      doSys = checkSignificance(avg, error, "difference");
    } else {
      hRB1 = DivideAndComputeRogerBarlow(h4_8, h4_10);;
      hRB2 = DivideAndComputeRogerBarlow(h6_10, h4_10);
      avg = computeArithmeticMean(hRB1, error);
      doSys = checkSignificance(avg, error, "ratio");
    }
    hRB1->Write(); hRB2->Write();

    TH1D* hSys = (TH1D*) hRB1->Clone();
    if(doSys){
      cout << "avg: " << avg << endl;
      if(s != "Sub") // subtract 1 for ratio case
        avg -= 1.0;
      for(int b = 1; b < hSys->GetNbinsX()+1; b++){
        hSys->SetBinContent(b, std::abs(avg));
      }
    } else {
      for(int b = 1; b < hSys->GetNbinsX()+1; b++){
        hSys->SetBinContent(b, 0.);
      }
    }
    hSys->SetName("SBvariationSys_" + s);
    res[i] = hSys;
  }
  return res;
}

std::array<TH1D*, 3> MCClosure(TString runnr){
  cout << "doing MC closure..." << endl;
  TFile* fData = new TFile("plots/"+runnr+".root", "READ");
  TFile* fMC = new TFile("plots/Closure"+runnr+".root", "READ");
  
  TDirectory* dMCClosure = fOutput->mkdir("MCClosure");
  dMCClosure->cd();

  std::array<TH1D*, 3> res;

  TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.95);
  for(int i = 0; i < signTypes.size(); i++){
    TString s = signTypes[i];
    TH1D *hData, *hMC;
    fData->GetObject("XiXi/XiXidphiCorrected"+s+"_pT_1.0_8.0", hData);
    fMC->GetObject("h"+s, hMC);

    if(!hData || !hMC){
      std::cout << "Could not retrieve all histograms for " << s << std::endl;
      continue;
    }

    // plot data vs MC
    hData->SetName("Data"+s); hMC->SetName("MC"+s);
    applyHistStyle(hData, 0);
    applyHistStyle(hMC, 1);
    hData->SetTitle((";"+xytitle).Data());
    leg->AddEntry(hData, "Reco", "l"); leg->AddEntry(hMC, "Gen", "l");
    
    // Draw using drawAndSavePdf
    std::vector<TH1*> histsMC = {hData, hMC};
    TCanvas *cMC = drawAndSavePdf(histsMC, makePDFs ? "figures/systematics/MCClosure_" + s + ".pdf" : "", false, 800, 800, "E", leg);
    cMC->SetName("MCClosure_"+s);
    dMCClosure->cd();
    cMC->Write();
    if(!makePDFs) delete cMC;
    leg->Clear();

    bool doSys = false;
    double avg, error;
    TH1D *hRB;
    // if(s == "Sub"){ // in case of subtracted histogram, do subtraction of variations instead of ratio (ratio --> 0/0)
    if(s == "Sub" || s == "SS"){ // in case of subtracted histogram, do subtraction of variations instead of ratio (ratio --> 0/0)
      hRB = SubtractAndComputeRogerBarlow(hMC, hData); hRB->Write();
      avg = computeArithmeticMean(hRB, error);
      doSys = checkSignificance(avg, error, "difference");
    } else {
      hRB = DivideAndComputeRogerBarlow(hMC, hData); hRB->Write();
      avg = computeArithmeticMean(hRB, error);
      doSys = checkSignificance(avg, error, "ratio");
    }
    TH1D* hSys = (TH1D*) hRB->Clone();
    if(doSys){
      cout << "avg: " << avg << endl;
      if(s != "Sub") // subtract 1 for ratio case
        avg -= 1.0;
      for(int b = 1; b < hSys->GetNbinsX()+1; b++){
        hSys->SetBinContent(b, std::abs(avg));
      }
    } else {
      for(int b = 1; b < hSys->GetNbinsX()+1; b++){
        hSys->SetBinContent(b, 0.);
      }
    }
    hSys->SetName("MCClosureSys_" + s);
    res[i] = hSys;
  }
  return res;
}

int plotSpectra(std::vector<TString> runnumbers, std::vector<TString> labels, TString name){
  spectraDir->cd();
  
  TCanvas *c = new TCanvas("cSpectra"+name, "cSpectra"+name, 800, 800);
  c->cd();
  TLegend *l = new TLegend(0.55, 0.70, 0.88, 0.85);

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
    hSpectra->SetTitle("Xi- Spectra;p_{T} (GeV/c);1/N_{evt} dN/dp_{T} (GeV/c)^{-1}");
    applyHistStyle(hSpectra, i % nColours);
    if (i == 0){
      hSpectra->Draw("E");
    } else {
      hSpectra->Draw("Esame");
    }
    l->AddEntry(hSpectra, labels[i]);
  }
  spectraDir->cd();
  l->Draw();
  c->Write();
  c->Clear(); l->Clear();

  TH1D* def;
  TFile *f = new TFile("spectra"+runnumbers[1]+".root", "READ");
  f->GetObject("hSpectraRebin", def);
  
  // Compute averages for Tight (index 2) and Loose (index 0) variations BEFORE fitting
  TFile *fTight = new TFile("spectra"+runnumbers[2]+".root", "READ");
  TFile *fLoose = new TFile("spectra"+runnumbers[0]+".root", "READ");
  TH1D *hTight, *hLoose;
  fTight->GetObject("hSpectraRebin", hTight); fLoose->GetObject("hSpectraRebin", hLoose);
  
  // Create ratio histograms for averaging using Roger Barlow errors
  TH1D *hRatioTight = DivideAndComputeRogerBarlow(hTight, def); TH1D *hRatioLoose = DivideAndComputeRogerBarlow(hLoose, def);
  
  // Compute averages using the helper functions
  double avgTight = 0.0, errTight = 0.0;
  double avgLoose = 0.0, errLoose = 0.0;
  
  avgTight = computeWeightedAverage(hRatioTight, errTight, 2); // Skip first bin
  avgLoose = computeWeightedAverage(hRatioLoose, errLoose, 2); // Skip first bin
  
  // Average of the two variations
  double avgBoth = (std::abs(avgTight - 1.0) + std::abs(avgLoose - 1.0)) / 2.0;
  cout << "plotSpectra: avgTight = " << avgTight << " +/- " << errTight << endl;
  cout << "plotSpectra: avgLoose = " << avgLoose << " +/- " << errLoose << endl;
  cout << "plotSpectra: average of both variations = " << avgBoth << endl;
  
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
    hSpectra->SetTitle("Xi- Spectra Ratio to Run "+runnumbers[1]+";p_{T} (GeV/c);Ratio");
    hSpectra->SetAxisRange(0.9, 1.1, "Y");
    applyHistStyle(hSpectra, i % nColours);
    if (i == 0){
      hSpectra->Draw("E");
    } else {
      hSpectra->Draw("Esame");
    }
    
    // Add straight line fit for Tight (i == 2) and Loose (i == 0) variations
    if(i == 0 || i == 2) {
      TString fitName = (i == 2) ? "fitTight" : "fitLoose";
      double fitValue = (i == 2) ? avgTight : avgLoose;
      double fitError = (i == 2) ? errTight : errLoose;
      
      // Draw error band around the fit line using TGraphErrors
      double xMin = hSpectra->GetXaxis()->GetXmin();
      double xMax = hSpectra->GetXaxis()->GetXmax();
      TGraphErrors *gFit = new TGraphErrors(2);
      gFit->SetPoint(0, xMin, fitValue); gFit->SetPoint(1, xMax, fitValue);
      gFit->SetPointError(0, 0, fitError); gFit->SetPointError(1, 0, fitError);
      gFit->SetLineColor(kStyleColors[i % nColours]); gFit->SetLineStyle(2); gFit->SetLineWidth(2);
      gFit->SetFillColorAlpha(kStyleColors[i % nColours], 0.2);
      gFit->Draw("E3 L SAME");
    }
    
    // Add legend entry with systematic error information
    if(i == 0) {
      l->AddEntry(hSpectra, Form("%s: %.2f%% #pm %.2f%%", labels[i].Data(), std::abs((avgLoose - 1.0) * 100.0), errLoose * 100.0), "lep");
    } else if(i == 2) {
      l->AddEntry(hSpectra, Form("%s: %.2f%% #pm %.2f%%", labels[i].Data(), std::abs((avgTight - 1.0) * 100.0), errTight * 100.0), "lep");
    } else {
      l->AddEntry(hSpectra, labels[i], "lep");
    }
  }
  spectraDir->cd();
  l->Draw();
  c->Write("ratio"+name);
  if(makePDFs) {
    c->SetTitle("");
    c->Print("figures/systematics/spectra_ratio_" + name + ".pdf");
  }
  c->Clear(); l->Clear();
  fOutput->cd();
  return 0;
}

void plotSystematics(std::array<TH1D*, 3> hSys, TString runnr){
  std::array<TString, 3> types = {"OS", "SS", "Sub"};
  
  for(int i = 0; i < 3; i++){
    if(!hSys[i]){
      cout << "plotSystematics: hSys[" << i << "] is null." << endl;
      continue;
    }

    TFile* f = new TFile("plots/" + runnr + ".root", "READ");
    if(!f || f->IsZombie()){
      cout << "plotSystematics: Could not open file plots/" << runnr << ".root" << endl;
      continue;
    }

    TString hname = "XiXi/XiXidphiCorrected" + types[i] + "_pT_1.0_8.0";
    TH1D* hResult = nullptr;
    f->GetObject(hname, hResult);
    if(!hResult){
      cout << "plotSystematics: Could not retrieve " << hname << endl;
      f->Close();
      continue;
    }

    if(hResult->GetNbinsX() != hSys[i]->GetNbinsX()){
      cout << "plotSystematics: Bin mismatch between result and systematic for " << types[i] << endl;
      f->Close();
      continue;
    }

    TCanvas *c = new TCanvas(("cSystematics_" + types[i]).Data(), ("Systematics " + types[i]).Data(), 800, 800);
    c->cd();

    applyHistStyle(hResult, 0);  // First color (black)
    hResult->SetTitle((";"+xytitle).Data());
    hResult->Draw("E");

    // Draw systematic errors as filled boxes using symmetric errors
    bool isSub = (types[i] == "Sub");
    TGraphErrors *gSys = new TGraphErrors(hResult->GetNbinsX());
    
    for(int b = 1; b <= hResult->GetNbinsX(); b++){
      double yVal = hResult->GetBinContent(b);
      double sys = hSys[i]->GetBinContent(b);
      
      if(!isSub)
        sys *= std::abs(yVal);
      
      gSys->SetPoint(b-1, hResult->GetXaxis()->GetBinCenter(b), yVal);
      gSys->SetPointError(b-1, hResult->GetXaxis()->GetBinWidth(b) / 2.0, sys);
    }
    
    gSys->SetFillColorAlpha(kBlue, 0.3); gSys->SetLineColor(kBlue); gSys->SetLineWidth(0);
    gSys->Draw("2 SAME");
    hResult->Draw("E SAME");

    bool isSS = (types[i] == "SS");
    TLegend *leg = isSS ? new TLegend(0.22, 0.80, 0.55, 0.88)
                        : new TLegend(0.55, 0.80, 0.88, 0.88);
    leg->AddEntry(hResult, "This thesis", "lep");
    leg->Draw();

    if(fOutput){
      fOutput->cd();
      c->Write();
    }
    if(makePDFs) {
      c->SetTitle("");
      c->Print("figures/systematics/systematics_" + types[i] + ".pdf");
    }
    f->Close();
  }
}

int systematics(){
  SetStyle();
  // gEnv->Print(); // debug
  cout << gStyle->GetLabelFont() << endl;
  
  fOutput = new TFile("systematics.root", "RECREATE");
  fOutput->cd();

  TCanvas *c = new TCanvas();
  c->cd();

  TH1D *h = new TH1D("h", "h", 10, 0, 10);
  h->FillRandom("gaus", 1000);
  h->Draw();
  TLatex *t = new TLatex(0.5, 0.5, "The quick brown fox jumps over the lazy dog");
  t->SetNDC();
  t->SetTextAlign(22);
  t->Draw();
  c->Print("mock_data.pdf");
  c->Clear();

  // histograms that contain all the individual systematic uncertainties
  std::array<TH1D*, 3> topoSys, SBSys, MCClosureSys, corrSys;

  // this is the run with the default settings. In case of variations in postprocessing, simply use this one. 
  TString standard = "577711";

  // // topological variations
  // TString topoLoose = "541739";
  // TString topoDefault = standard;
  // TString topoTight = "541740";

  // TString topoLoose = "577712";
  // TString topoDefault = "577711";
  // TString topoTight = "578138";  

  // The ones above give a ~10% sys variation, so we try again with better cuts
  TString topoLoose = "607799";
  TString topoDefault = "607798";
  TString topoTight = "607800";

  topoSys = topo(topoDefault, topoLoose, topoTight);

  // to investigate some of the effect of weighting on the Roger-Barlow uncertainties, here we compare the variations without ME weighting
  // note that efficiency corrections are still applied, so the datapoints are still weighted
  // topoSys = topo("607798NoME", "607799NoME", "607800NoME");

  // check whether sideband correction is significant
  // as the comparison between with and without correction uses the same dataset, the statistical uncertainty is fully correlated.
  // Hence, we use the Roger Barlow method to determine whether the difference is significant
  SBcheck(standard);

  // vary the mass window for the sideband to estimate the systematic uncertainty from the sideband subtraction
  TString SB_4_10 = "607798_SB_4_10";
  TString SB_4_8 = "607798_SB_4_8";
  TString SB_6_10 = "607798_SB_6_10";
  SBSys = SBvariations(SB_4_10, SB_4_8, SB_6_10);

  TString closureRun = "577062";
  // MCClosureSys = MCClosure(closureRun);
  MCClosureSys = MCClosure("617914");

  // Plot spectra with different variations
  spectraDir = fOutput->mkdir("spectra");
  plotSpectra({"603231", "603230", "603229"}, CascRadLabels, "CascRad");
  plotSpectra({"603228", "603227", "603226"}, DCAPosNegToPVLabels, "DCAPosNegToPV");

  // plotSpectra({"608045", "608044", "608046"}, varsLabels, "sysvars"); // small dataset
  plotSpectra({"617288", "617289", "617287"}, varsLabels, "sysvars"); // full 24 dataset

  // now combine all systematics
  fOutput->cd();
  TCanvas *c1 = new TCanvas();
  c1->cd();
  TLegend *leg = new TLegend(0.55, 0.75, 0.88, 0.95);

  // grab the results from the default run
  TFile* f = new TFile("plots/" + topoDefault + ".root", "READ");
  std::array<TH1D*, 3> hCorrFunc;
  for(int i = 0; i < signTypes.size(); i++){
    TString hname = "XiXi/XiXidphiCorrected" + signTypes[i] + "_pT_1.0_8.0";
    hCorrFunc[i] = nullptr;
    f->GetObject(hname, hCorrFunc[i]);
  }

  // add all the different sources of systematic uncertainties in quadrature
  double materialBudgetSys = 0.03;
  fOutput->cd();
  for(int i = 0; i < 3; i++) {
    TH1D* h = topoSys[i];
    for(int n = 1; n < h->GetNbinsX()+1; n++) {
      double sys = 0;
      if(i == 2) { // absolute value in case of subtracted
        sys += std::pow(0.03*hCorrFunc[2]->GetBinContent(n), 2);
      } else {
        sys += std::pow(materialBudgetSys, 2);
      }
      sys += std::pow(topoSys[i]->GetBinContent(n), 2);
      sys += std::pow(SBSys[i]->GetBinContent(n), 2);
      sys += std::pow(MCClosureSys[i]->GetBinContent(n), 2);
      h->SetBinContent(n, std::sqrt(sys));
    }

    h->SetName("CorrSystematic_" + signTypes[i]);
    corrSys[i] = h;
    h->Write();
  }

  // let's do it another way - from the spectra
  double topoSpectraSys = 0.015;
  double sysvalue = std::sqrt(std::pow(materialBudgetSys, 2) + std::pow(topoSpectraSys, 2));
  std::array<TH1D*, 3> spectraSys;
  for(int i = 0; i < 2; i++){
    TH1D* h = new TH1D(*topoSys[i]);
    h->SetName("SpectraSystematic_" + signTypes[i]);
    for(int bin = 1; bin <= h->GetNbinsX(); bin++){
      h->SetBinContent(bin, std::sqrt(sysvalue*sysvalue + std::pow(SBSys[i]->GetBinContent(bin), 2)));
    }
    spectraSys[i] = h;
  }
  spectraSys[2] = new TH1D(*topoSys[2]);
  spectraSys[2]->SetName("SpectraSystematic_Sub");

  for(int bin = 1; bin <= spectraSys[2]->GetNbinsX(); bin++){
    double val = std::abs(hCorrFunc[0]->GetBinContent(bin)*spectraSys[0]->GetBinContent(bin) - hCorrFunc[1]->GetBinContent(bin)*spectraSys[1]->GetBinContent(bin));
    spectraSys[2]->SetBinContent(bin, std::sqrt(val*val + std::pow(SBSys[2]->GetBinContent(bin), 2)));
  }

  std::array<TH1D*, 3> totalSys;
  for(int i = 0; i < 2; i++){
    totalSys[i] = new TH1D(*spectraSys[i]);
    totalSys[i]->SetName("TotalSystematic_" + signTypes[i]);
  }

  totalSys[2] = new TH1D(*corrSys[2]);
  totalSys[2]->SetName("TotalSystematic_Sub");
  for(int bin = 1; bin <= totalSys[2]->GetNbinsX(); bin++){
    double corrVal = corrSys[2]->GetBinContent(bin);
    double spectraVal = spectraSys[2]->GetBinContent(bin);
    totalSys[2]->SetBinContent(bin, TMath::Max(corrVal, spectraVal));
  }

  // plotSystematics(corrSys, topoDefault);
  // plotSystematics(spectraSys, topoDefault);
  plotSystematics(totalSys, topoDefault);

  fOutput->Write();
  fOutput->Close();
  return 0;
} 