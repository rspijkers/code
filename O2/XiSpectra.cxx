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
#include "TLine.h"
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

// using json = nlohmann::json;
using std::cout; using std::endl;

struct mass{
  enum{
    invMass, 
    sign,
    pT, 
    eta, 
    V_z,
    multiplicity
  }; 
};

double rapidityFromEta(double eta, double pt, double m){
    /*Converts pseudorapidity to rapidity.*/
    double cosheta = cosh(eta);
    double m2 = m*m;
    double pt2 = pt*pt;

    return log((sqrt(m2 + pt2*cosheta*cosheta) + pt*sinh(eta)) / sqrt(m2 + pt2));
}

double pTmax = 6.5; //fixme

struct XiFitResults {
  std::vector<double> purity;
  std::vector<double> mu;
  std::vector<double> sigma;
};

XiFitResults doXiInvMassFitsXiBinning(TH3F *hXiMass, const double *xiBinning, int nXiBins,
                                      TFile *outputFile, double signalWindow = 3.0, bool makePDF = false) {
  XiFitResults results;
  results.purity.resize(nXiBins - 1);
  results.mu.resize(nXiBins - 1);
  results.sigma.resize(nXiBins - 1);
  
  TDirectory *XiInvMass = outputFile->mkdir("XiInvMass");
  XiInvMass->cd();

  TF1 *f1 = new TF1("f1", "pol2(0) + gaus(3) + gaus(6)", 1.29, 1.36);
  f1->SetParameters(0, 0, 10, 0, 1.321, 0.005, 0, 1.321, 0.001);
  f1->SetParLimits(4, 1.31, 1.33);
  f1->SetParLimits(7, 1.31, 1.33);
  f1->SetParLimits(5, 0, 0.01);
  f1->SetParLimits(8, 0, 0.01);

  TF1 *fBKGXi = new TF1("fBKGXi", pol2bkgXi, 1.29, 1.36, 3);

  hXiMass->SetDirectory(outputFile);
  cout << "Start Xi mass fitting (Xi binning)..." << endl;

  for (int pTbin = 0; pTbin < nXiBins - 1; pTbin++) {
    const double ptLow = xiBinning[pTbin];
    const double ptHigh = xiBinning[pTbin + 1];
    hXiMass->GetYaxis()->SetRangeUser(ptLow, ptHigh);

    TH1D *h = (TH1D *)hXiMass->Project3D("xe");
    TString lowLabel = TString::Format("%.1f", ptLow);
    TString highLabel = TString::Format("%.1f", ptHigh);
    TString canvasName = "cXiInvMass_" + lowLabel + "_" + highLabel;
    TCanvas *c = new TCanvas(canvasName, canvasName);
    h->SetName("hMassXi_" + lowLabel + "_" + highLabel);
    h->SetTitle(";#Lambda#pi inv. mass (GeV/#it{c^{2}}); counts");

    f1->SetParameter(3, .4 * h->GetMaximum());
    f1->SetParameter(6, .4 * h->GetMaximum());
    h->Fit(fBKGXi, "SLBQRO");
    f1->SetParameter(0, fBKGXi->GetParameter(0));
    f1->SetParameter(1, fBKGXi->GetParameter(1));
    f1->SetParameter(2, fBKGXi->GetParameter(2));

    TFitResultPtr r1 = h->Fit("f1", "SLBQR", "", 1.29, 1.36);
    double chi2Reduced = r1->Chi2() / r1->Ndf();
    cout << "Reduced Chi2 = " << chi2Reduced << " in pT bin " << pTbin << endl;

    double mu = (f1->GetParameter(4) + f1->GetParameter(7)) / 2.;
    double sigma = (f1->GetParameter(5) + f1->GetParameter(8)) / 2.;

    h->GetXaxis()->SetRangeUser(1.28, 1.38);
    h->SetStats(kFALSE);

    // compute purity
    TF1 *fPol2Xi = new TF1("fPol2Xi", "pol2", 1.29, 1.36);
    fPol2Xi->SetParameters(fBKGXi->GetParameter(0), fBKGXi->GetParameter(1), fBKGXi->GetParameter(2));
    double total = h->Integral(h->GetXaxis()->FindBin(mu - signalWindow * sigma),
                               h->GetXaxis()->FindBin(mu + signalWindow * sigma));
    double bkg = fPol2Xi->Integral(mu - signalWindow * sigma, mu + signalWindow * sigma) /
                 h->GetXaxis()->GetBinWidth(1);
    double purity = (total - bkg) / total;
    
    // Store fit results for later use
    results.purity[pTbin] = purity;
    results.mu[pTbin] = mu;
    results.sigma[pTbin] = sigma;

    if (makePDF) {
      h->Draw();
      // vertical lines for visualization of signal region only
      TLine *sig3low = new TLine(mu - signalWindow * sigma, 0, mu - signalWindow * sigma, h->GetMaximum());
      TLine *sig3high = new TLine(mu + signalWindow * sigma, 0, mu + signalWindow * sigma, h->GetMaximum());
      for (auto line : {sig3low, sig3high}) {
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->Draw("same");
      }
      // plot bkg function, also under signal peak
      fPol2Xi->SetLineColor(kGreen + 2);
      fPol2Xi->Draw("same");

      TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
      leg->AddEntry(h, "#Xi invariant mass", "l");
      leg->AddEntry(fPol2Xi, "Background fit function", "l");
      leg->AddEntry(f1, "Total fit function", "l");
      leg->Draw();

      TLatex *purityText = new TLatex();
      purityText->SetNDC();
      purityText->SetTextSize(0.03);
      purityText->DrawLatex(0.55, 0.60, Form("Purity: %.2f%%", purity * 100));

      c->Write();
      c->Clear();
    }
  }

  outputFile->cd();
  
  // Reset histogram ranges to full range
  hXiMass->GetXaxis()->SetRange(0, 0);
  hXiMass->GetYaxis()->SetRange(0, 0);
  hXiMass->GetZaxis()->SetRange(0, 0);
  
  return results;
}

int XiSpectra(TString trainnr, TString mctrainnr, TString filename = "AnalysisResults.root", bool makePDF = false, double signalWindow = 3.0){
  // make single Xi spectra (charge independent)

  // Enable OpenGL for proper alpha transparency support
  // gEnv->SetValue("OpenGL.CanvasPreferGL", 1);
  gStyle->SetCanvasPreferGL(kTRUE);

  TFile* inputFile = new TFile("results/" + trainnr + "/" + filename, "READ");
  THnSparse *hEffCorrXiMass;
  inputFile->GetObject("cascade-correlations/hMassXiEffCorrected", hEffCorrXiMass);
  TFile *mcFile = new TFile("results/" + mctrainnr + "/" + filename, "READ");

  double nGenEvents = mcFile->Get<TH1F>("cascade-selector/gen/hNevents")->GetEntries(); //get this from cascade-selector/gen/hNEvents (get entries? or get bincontent?)
  double nRecEvents = mcFile->Get<TH1F>("cascade-selector/hEventSel")->GetBinContent(6); //get this from cascade-selector/hEventSel->GetBinContent(6) (= selected events)
  double nEvents = inputFile->Get<TH1F>("cascade-selector/hEventSel")->GetBinContent(6);
  
  // these are the values for runs 522962 (data) and 521456 (MC)
  // double nRecEvents = 1209904670;
  // double nGenEvents = 2024331520;
  // double nEvents = 1002102822;

  cout << "NRecEvents from MC: " << nRecEvents << endl;
  cout << "NGenEvents from MC: " << nGenEvents << endl;
  cout << "NEvents from Data: " << nEvents << endl;

  TFile *outputFile = new TFile("spectra" + trainnr + ".root", "RECREATE");
  
  // Create a new THnSparse with your selection criteria
  THnSparse *hFiltered = (THnSparse*)hEffCorrXiMass->Clone("hFiltered");
  hFiltered->Reset();
  // Iterate through all bins in the sparse histogram
  for(Long64_t i = 0; i < hEffCorrXiMass->GetNbins(); i++){
    Int_t binIndices[6];
    Double_t value = hEffCorrXiMass->GetBinContent(i, binIndices);
    Double_t error = hEffCorrXiMass->GetBinError(i);
    if(value == 0) continue; // Skip empty bins
    
    // Get axis values for this bin
    // Double_t m = hEffCorrXiMass->GetAxis(mass::invMass)->GetBinCenter(binIndices[mass::invMass]);
    Double_t pt = hEffCorrXiMass->GetAxis(mass::pT)->GetBinCenter(binIndices[mass::pT]);
    Double_t eta = hEffCorrXiMass->GetAxis(mass::eta)->GetBinCenter(binIndices[mass::eta]);

    if(rapidityFromEta(eta, pt, 1.32171) > -0.5 && rapidityFromEta(eta, pt, 1.32171) < 0.5){
      hFiltered->SetBinContent(binIndices, value);
      hFiltered->SetBinError(binIndices, error);
    }
  }

  // 2D projections of hFiltered on pT and eta (with and without basic inv mass selection)
  TH2D *hFilteredPtEta = (TH2D*)hFiltered->Projection(mass::pT, mass::eta);
  hFilteredPtEta->SetName("hFilteredPtEta");

  const int invMassBinLowBasic = hFiltered->GetAxis(mass::invMass)->FindBin(1.31);
  const int invMassBinHighBasic = hFiltered->GetAxis(mass::invMass)->FindBin(1.33);
  hFiltered->GetAxis(mass::invMass)->SetRange(invMassBinLowBasic, invMassBinHighBasic);
  TH2D *hFilteredPtEtaSignal = (TH2D*)hFiltered->Projection(mass::pT, mass::eta);
  hFilteredPtEtaSignal->SetName("hFilteredPtEtaSignal");
  hFiltered->GetAxis(mass::invMass)->SetRange(0, 0);

  TH1D* hXiSpectra = project(hFiltered, mass::pT, {{mass::pT, {0., pTmax}}, {mass::invMass, {1.31, 1.33}}});
  hXiSpectra->SetName("hXiSpectra");
  hXiSpectra->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  Double_t Xibinning[14] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};
  
  TH3F *hMassXiMinus = inputFile->Get<TH3F>("cascade-correlations/hMassXiMinus");
  outputFile->cd();
  TH2D *hMassXiMinusPtEta = (TH2D*)hMassXiMinus->Project3D("yz");
  hMassXiMinusPtEta->SetName("hMassXiMinusPtEta");

  const int invMassBinLow = hMassXiMinus->GetXaxis()->FindBin(1.31);
  const int invMassBinHigh = hMassXiMinus->GetXaxis()->FindBin(1.33);
  hMassXiMinus->GetXaxis()->SetRange(invMassBinLow, invMassBinHigh);
  TH2D *hMassXiMinusPtEtaSignal = (TH2D*)hMassXiMinus->Project3D("yz");
  hMassXiMinusPtEtaSignal->SetName("hMassXiMinusPtEtaSignal");
  hMassXiMinus->GetXaxis()->SetRange(0, 0);
  XiFitResults fitResults = doXiInvMassFitsXiBinning(hMassXiMinus, Xibinning, 14, outputFile, signalWindow, makePDF);
  
  // Build spectra bin-by-bin using fitted signal windows
  TH1D* hSpectraRebin = new TH1D("hSpectraRebin", "", 13, Xibinning);
  for (int pTbin = 0; pTbin < 13; pTbin++) {
    const double ptLow = Xibinning[pTbin];
    const double ptHigh = Xibinning[pTbin + 1];
    const double invMassLow = fitResults.mu[pTbin] - signalWindow * fitResults.sigma[pTbin];
    const double invMassHigh = fitResults.mu[pTbin] + signalWindow * fitResults.sigma[pTbin];
    
    // Use project() helper function to get spectrum in this pT and mass window
    TH1D* hTemp = project(hFiltered, mass::pT, {{mass::pT, {ptLow, ptHigh}}, {mass::invMass, {invMassLow, invMassHigh}}});
    double error = 0;
    double integral = hTemp->IntegralAndError(1, hTemp->GetNbinsX(), error);
    
    // Correct for purity (multiply by purity)
    integral *= fitResults.purity[pTbin];
    error *= fitResults.purity[pTbin];
    
    hSpectraRebin->SetBinContent(pTbin + 1, integral);
    hSpectraRebin->SetBinError(pTbin + 1, error);
    
    delete hTemp;
  }
  
  hSpectraRebin->Scale(1./(nEvents), "width"); // divide by an extra factor of 2 because of Xi charges ---- 65536 = factor of 2^n used in run 3 AN
  // hSpectraRebin->Scale(1209904670./2024331520.); // Nev(rec) = 1209904670 / Nev(gen) = 2024331520 (run 521456)
  // hSpectraRebin->Scale(1199241785./2019676955.); // Nev(rec) = 1199241785 / Nev(gen) = 2019676955 (run 550154 24f3c_fix_medium)
  hSpectraRebin->Scale(nRecEvents/nGenEvents); // correct for the different number of generated and reconstructed events in MC
  cout << "Event reconstruction efficiency is " << nRecEvents/nGenEvents << endl;
  TH1D *hXiMinusMass1D = hMassXiMinus->ProjectionX();
  hXiMinusMass1D->SetName("hXiMinusMass1D");
  TH1D *hRawXiMass = hMassXiMinus->ProjectionY("hRawXiMass");
  TH1D* hRawSpectraRebin = (TH1D*) hRawXiMass->Rebin(13, "hRawSpectraRebin", Xibinning);
  hRawSpectraRebin->Scale(1./nEvents, "width"); // ---- 65536 = factor of 2^n used in run 3 AN

  TH1D* hXiEta = project(hEffCorrXiMass, mass::eta, {{mass::pT, {0.6, pTmax}}, {mass::invMass, {1.31, 1.33}}});
  hXiEta->SetName("hXiEta");

  TH1D* hXiSpectraMult = project(hEffCorrXiMass, mass::multiplicity, {{mass::pT, {0.6, pTmax}}, {mass::invMass, {1.31, 1.33}}, {mass::multiplicity, {0,99}}});
  hXiSpectraMult->SetName("hXiSpectraMult");

  TFile *run2file = new TFile("run2Xispectra.root", "READ");
  outputFile->cd();
  TH1D *run2Spectra = (TH1D*) run2file->Get<TH1F>("Table 3/Hist1D_y11");
  run2Spectra->SetLineColor(kRed);
  // get the error bars here as well:
  // stat = e1, sys,total = e2
  // set the stat errors as bin errors, add the systematic errors as boxes later on
  TH1D *run2StatError = (TH1D*) run2file->Get<TH1F>("Table 3/Hist1D_y11_e1");
  TH1D *run2SysError = (TH1D*) run2file->Get<TH1F>("Table 3/Hist1D_y11_e2");
  for(Int_t i = 1; i <= run2Spectra->GetNbinsX(); i++){
    run2Spectra->SetBinError(i, run2StatError->GetBinContent(i));
  }
  hSpectraRebin->GetYaxis()->SetRangeUser(1e-5, 0.045); // don't put minimum to zero, so we can use logscale later
  TCanvas *cSpectra = new TCanvas("cSpectra");
  hSpectraRebin->SetTitle("Xi spectra;#it{p}_{T} (GeV/#it{c});1/N_{ev} d^{2}N/(d#it{y} d#it{p}_{T}) (GeV/#it{c})^{-1}");
  hSpectraRebin->SetStats(kFALSE);
  hSpectraRebin->Draw();
  run2Spectra->Draw("SAME");
  TLegend *legend = new TLegend(0.65, 0.75, 0.89, 0.89);
  legend->AddEntry(hSpectraRebin, "This analysis");
  legend->AddEntry(run2Spectra, "EPJC 2020 7673");
  legend->Draw();
  cSpectra->Write();
  TH1D* hSpectraRatio = (TH1D*) hSpectraRebin->Clone();
  hSpectraRatio->SetName("hSpectraRatio");
  hSpectraRatio->GetXaxis()->SetLimits(0.6, 6.5);
  hSpectraRatio->Divide(run2Spectra);
  hSpectraRatio->GetYaxis()->SetTitle("ratio to run 2");

  hSpectraRebin->SetLineWidth(2);
  run2Spectra->SetLineWidth(2);
  hSpectraRatio->SetLineWidth(2);

  // Material budget systematics from run 2:
  std::vector<double> materialBudgetSys = {0.051, 0.036, 0.029, 0.024, 0.021, 0.019, 0.018, 0.016, 0.015, 0.012, 0.008, 0.005, 0.005};

  // Create histogram for material budget systematic uncertainty
  TH1D *hMaterialBudgetSys = new TH1D("hMaterialBudgetSys", "Material budget systematic uncertainty;#it{p}_{T} (GeV/#it{c});Relative systematic uncertainty", 13, Xibinning);
  for(int i = 0; i < 13; i++){
    hMaterialBudgetSys->SetBinContent(i+1, materialBudgetSys[i]);
  }
  hMaterialBudgetSys->SetLineColor(kBlue);
  hMaterialBudgetSys->SetLineWidth(2);

  // Create TGraphErrors for hSpectraRebin with 5% systematic errors
  TGraphErrors *hSpectraSysGraph = new TGraphErrors();
  for(int i = 0; i < hSpectraRebin->GetNbinsX(); i++){
    hSpectraSysGraph->SetPoint(i, hSpectraRebin->GetBinCenter(i+1), hSpectraRebin->GetBinContent(i+1));
    double sysErr = std::sqrt(0.015*0.015 + materialBudgetSys[i]*materialBudgetSys[i]) * hSpectraRebin->GetBinContent(i+1);  // 1.5% systematic error
    hSpectraSysGraph->SetPointError(i, hSpectraRebin->GetBinWidth(i+1)/2, sysErr);
  }

  // Create TGraphErrors for run 2 systematic errors
  TGraphErrors *run2SysGraph = new TGraphErrors();
  for(int i = 0; i < run2SysError->GetNbinsX(); i++){
    run2SysGraph->SetPoint(i, run2Spectra->GetBinCenter(i+1), run2Spectra->GetBinContent(i+1));
    run2SysGraph->SetPointError(i, run2Spectra->GetBinWidth(i+1)/2, run2SysError->GetBinContent(i+1));
  }

  // Okay, TRatioPlot sucks, I'll do it myself  
  TCanvas *cRatio = plotWithRatio("cRatio", hSpectraRebin, run2Spectra, hSpectraRatio, 
                                   "2", "2", "E", legend, 800, 800, hSpectraSysGraph, run2SysGraph);
  TPad *padUpper = (TPad*)cRatio->GetPrimitive("padUpper");
  TPad *padLower = (TPad*)cRatio->GetPrimitive("padLower");
  padUpper->SetLogy();

  cRatio->Write();
  outputFile->Write();
  
  // Save spectra pdfs to figures/spectra
  // save ratio plot first, before messing with y-axis range and log scale of the upper pad
  hSpectraRebin->SetTitle("");
  run2Spectra->SetTitle("");
  cRatio->SaveAs("figures/spectra/spectra_ratio.pdf");

  cSpectra->cd();
  cSpectra->SetLeftMargin(0.15);
  hSpectraRebin->SetMarkerStyle(0);
  run2Spectra->SetMarkerStyle(0);
  hSpectraRebin->GetYaxis()->SetRangeUser(0, 0.025);
  cSpectra->SaveAs("figures/spectra/spectra.pdf");
  
  // Create logarithmic version - need to redraw after changing to log scale
  cSpectra->SetLogy();
  hSpectraRebin->GetYaxis()->SetRangeUser(1e-5, 0.045);
  cSpectra->SaveAs("figures/spectra/spectra_log.pdf");

  outputFile->Close();
  inputFile->Close();

  return 0;
}