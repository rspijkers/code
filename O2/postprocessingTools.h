// std
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <cmath>
// // json parsing
// #include <nlohmann/json.hpp>
// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TString.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TGrid.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TRatioPlot.h"

// using json = nlohmann::json;
using std::cout; using std::endl;
using axranges = std::map<int, std::vector<double>>;

// Preferred colors and markers for consistent plotting style
static const Int_t kStyleFillColors[] = {kGray+1, kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9, kCyan-8, kYellow-7};
static const Int_t kStyleColors[] = {kBlack, kRed+1, kBlue+1, kGreen+3, kMagenta+1, kOrange-1, kCyan+2, kYellow+2};
static const Int_t kStyleMarkers[] = {kFullCircle, kFullSquare, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross, kFullCross, kFullDiamond, kFullStar, kOpenStar};

inline void myPadSetUp(TPad *currentPad, float currentLeft=0.15, float currentTop=0.04, float currentRight=0.04, float currentBottom=0.15){
  if(!currentPad) return;
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);
}

inline void SetStyle(Bool_t graypalette=kFALSE) {
  gStyle->Reset("Plain");
  gROOT->ForceStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(2);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.034,"xyz");
  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.038,"xyz");
  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleOffset(1.6,"y");
  gStyle->SetTitleOffset(1.0,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.026);
  gStyle->SetStripDecimals(kFALSE);
  TGaxis::SetMaxDigits(4);
  TGaxis::SetExponentOffset(0.0f, 0.01f, "y");
}

inline void applyHistStyle(TH1 *hist,
                           int colorIdx = 0,
                           int markerIdx = 0,
                           int lineStyle = 1,
                           int fillColorIdx = -1) {
  if(!hist) return;

  const int nColors = sizeof(kStyleColors) / sizeof(kStyleColors[0]);
  const int nMarkers = sizeof(kStyleMarkers) / sizeof(kStyleMarkers[0]);
  const int nFillColors = sizeof(kStyleFillColors) / sizeof(kStyleFillColors[0]);

  int cIdx = ((colorIdx % nColors) + nColors) % nColors;
  int mIdx = ((markerIdx % nMarkers) + nMarkers) % nMarkers;

  hist->SetLineColor(kStyleColors[cIdx]);
  hist->SetMarkerColor(kStyleColors[cIdx]);
  hist->SetMarkerStyle(kStyleMarkers[mIdx]);
  hist->SetLineStyle(lineStyle);

  if(fillColorIdx >= 0) {
    int fIdx = ((fillColorIdx % nFillColors) + nFillColors) % nFillColors;
    hist->SetFillColor(kStyleFillColors[fIdx]);
  }
}

inline TCanvas* drawAndSavePdf(const std::vector<TH1*>& hists,
                                         const TString& pdfPath,
                                         bool logy = false,
                                         int canvasWidth = 800,
                                         int canvasHeight = 800,
                                         const char* drawOpt = "E",
                                         TLegend* legend = nullptr) {
  if (hists.empty()) {
    cout << "drawAndSavePdf: empty histogram list" << endl;
    return nullptr;
  }

  for (size_t i = 0; i < hists.size(); ++i) {
    if (!hists[i]) {
      cout << "drawAndSavePdf: histogram at index " << i << " is null" << endl;
      return nullptr;
    }
  }
  TH1* firstHist = hists[0];

  static int canvasCounter = 0;
  TString canvasName = Form("cDrawHists_%d", canvasCounter++);
  TCanvas* c = new TCanvas(canvasName.Data(), canvasName.Data(), canvasWidth, canvasHeight);
  myPadSetUp((TPad*)c);
  c->SetLogy(logy);
  c->cd();

  // Determine global y-range from all provided histograms
  double globalMin = logy ? 1e-6 : 1e30;  // Initialize to +inf so first data point sets it
  double globalMax = logy ? 1e-6 : -1e30;  // Initialize to -inf for linear, so first value will update it
  for (auto* h : hists) {
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      const double y = h->GetBinContent(b);
      const double e = h->GetBinError(b);
      double low = y - e;
      double high = y + e;
      if (logy && high <= 0.) continue;
      if (logy && low <= 0.) low = 1e-12;
      if (high < low) continue;
      if (low < globalMin) globalMin = low;
      if (high > globalMax) globalMax = high;
    }
  }

  // Set sensible defaults if no content was found
  if (globalMax < globalMin) {
    globalMin = logy ? 1e-6 : 0.;
    globalMax = 1.;
  }

  if (logy) {
    const double minSafe = std::max(globalMin, 1e-12);
    firstHist->GetYaxis()->SetRangeUser(minSafe * 0.8, globalMax * 3.0);
  } else {
    const double span = globalMax - globalMin;
    const double baseline = (span > 0.) ? span : 0.2 * std::max(std::abs(globalMax), 1.0);
    double padLow  = 0.15 * baseline;
    double padHigh = legend ? 0.50 * baseline : 0.15 * baseline;
    firstHist->GetYaxis()->SetRangeUser(globalMin - padLow, globalMax + padHigh);
  }

  firstHist->Draw(drawOpt);
  for (size_t i = 0; i < hists.size(); ++i) {
    if (i == 0) continue;
    TString opt = TString(drawOpt) + " SAME";
    hists[i]->Draw(opt.Data());
  }

  if (legend) {
    legend->Draw();
  }

  c->SaveAs(pdfPath.Data());
  return c;
}

// Function for bkg fitting
double pol2bkgom(double *x, double *par){
  if (x[0] > 1.655 && x[0] < 1.689) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Function for bkg fitting
double pol2bkgXi(double *x, double *par){
  if (x[0] > 1.305 && x[0] < 1.335) {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Function that creates a projection to TargetAxis, with ranges in other dimensions
TH1D *project(THnSparse *THn,             // input THn
              int targetAxis,             // axis nr on which to project
              axranges map,               // map with [axisnr, {lower, upper}] bounds
              Option_t *options = "") {   // any options to pass to THnSparse->Projection()
  // TH1::SetDefaultSumw2(); // Make sure we propagate the errors
  int ndim = THn->GetNdimensions();
  for(auto const& [axis, bounds] : map){
    TAxis *ax = THn->GetAxis(axis);
    if (!ax) cout << "Error: GetAxis() with number " << axis << " gives a null pointer! The program will crash!" << endl; // null pointer check
    if(bounds[1] == ax->GetBinLowEdge(ax->FindFixBin(bounds[1]))){
      // if the upper boundary is at the lower edge of a bin, we actually want to integrate up to the bin before
      ax->SetRange(ax->FindFixBin(bounds[0]), ax->FindFixBin(bounds[1]) - 1);
    } else {
      ax->SetRangeUser(bounds[0], bounds[1]);
    }
  }
  TH1D *hp = THn->Projection(targetAxis, "E");
  TString axisname = THn->GetAxis(targetAxis)->GetTitle();
  hp->SetTitle("Projection on " + axisname);

  // reset axis ranges
  for(int i=0; i<ndim; i++) {
    TAxis *ax = THn->GetAxis(i);
    ax->SetRange(0,0);
  }
  return hp;
}

TH2D *project2D(THnSparse *THn,             // input THn
                int targetAxis1,            // axis 1 on which to project
                int targetAxis2,            // axis 2 on which to project
                axranges map,               // map with [axisnr, {lower, upper}] bounds
                Option_t *options = "") {   // any options to pass to THnSparse->Projection()
  // TH1::SetDefaultSumw2(); // Make sure we propagate the errors
  int ndim = THn->GetNdimensions();
  for(auto const& [axis, bounds] : map){
    TAxis *ax = THn->GetAxis(axis);
    if (!ax) cout << "Error: GetAxis() with number " << axis << " gives a null pointer! The program will crash!" << endl; // null pointer check
    if(bounds[1] == ax->GetBinLowEdge(ax->FindFixBin(bounds[1]))){
      // if the upper boundary is at the lower edge of a bin, we actually want to integrate up to the bin before
      ax->SetRange(ax->FindFixBin(bounds[0]), ax->FindFixBin(bounds[1]) - 1);
    } else {
      ax->SetRangeUser(bounds[0], bounds[1]);
    }
  }
  TH2D *hp = THn->Projection(targetAxis1, targetAxis2, "E");
  TString axis1 = THn->GetAxis(targetAxis1)->GetTitle();
  TString axis2 = THn->GetAxis(targetAxis2)->GetTitle();
  hp->SetTitle("Projection on " + axis1 + ", " + axis2);

  // reset axis ranges
  for(int i=0; i<ndim; i++) {
    TAxis *ax = THn->GetAxis(i);
    ax->SetRange(0,0);
  }
  return hp;
}

TH1D* DivideAndComputeRogerBarlow(TH1D* h1, TH1D *h2){ 
  // Use Roger Barlow "sigma_{delta}" as errors for ratios
  // Returns a new histogram with the ratio and Roger Barlow errors
  Double_t lh1NBins = h1->GetNbinsX(); 
  Double_t lh2NBins = h2->GetNbinsX(); 

  if( lh1NBins != lh2NBins ){ 
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return nullptr;
  }

  // Create a new histogram as a clone of h1
  TH1D* hRatio = (TH1D*)h1->Clone();
  hRatio->SetDirectory(nullptr);

  Double_t lSigmaDelta[100]; 
  for(Int_t i = 1; i < h1->GetNbinsX() + 1; i++){ 
    // Computation of roger barlow sigma_{delta} 
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    // Computation of relationship to h2 for plotting in ratio plot 
    if (h2->GetBinContent(i) > 1e-12) { 
      lSigmaDelta[i] /= h2->GetBinContent(i); 
    } else { 
      lSigmaDelta[i] = 0; 
    }
  }
  // Regular Division 
  hRatio->Divide(h2);
  // Replace Errors
  for(Int_t i = 1; i < hRatio->GetNbinsX() + 1; i++){
    hRatio->SetBinError(i, lSigmaDelta[i]);
  }
  return hRatio;
}

TH1D* SubtractAndComputeRogerBarlow(TH1D* h1, TH1D *h2){ 
  // Use Roger Barlow "sigma_{delta}" as errors for differences
  // Returns a new histogram with the difference and Roger Barlow errors
  Double_t lh1NBins = h1->GetNbinsX(); 
  Double_t lh2NBins = h2->GetNbinsX(); 

  if( lh1NBins != lh2NBins ){ 
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return nullptr;
  }

  // Create a new histogram as a clone of h1
  TH1D* hDiff = (TH1D*)h1->Clone();
  hDiff->SetDirectory(nullptr);

  Double_t lSigmaDelta[100]; 
  for(Int_t i = 1; i < h1->GetNbinsX() + 1; i++){ 
    // Computation of roger barlow sigma_{delta} 
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
  }
  // Regular Subtraction 
  hDiff->Add(h2, -1);
  // Replace Errors
  for(Int_t i = 1; i < hDiff->GetNbinsX() + 1; i++){
    hDiff->SetBinError(i, lSigmaDelta[i]);
  }
  return hDiff;
}

double computeArithmeticMean(TH1D* h, double &error, int minBin = 1, int maxBin = -1) {
  // Compute the arithmetic mean of histogram bins in range [minBin, maxBin] and its error
  // Default: use all bins (maxBin = -1 means use GetNbinsX())
  if (minBin < 1) minBin = 1;
  if (maxBin < 0) maxBin = h->GetNbinsX();
  if (maxBin > h->GetNbinsX()) maxBin = h->GetNbinsX();
  
  int nBins = maxBin - minBin + 1;
  if (nBins <= 0) {
    error = 0.0;
    return 0.0;
  }
  
  double avg = h->IntegralAndError(minBin, maxBin, error);
  avg /= nBins;
  error /= nBins;
  return avg;
}

bool checkSignificance(double avg, double error, TString method = "ratio"){
  // Check if a value is significantly different from the expected value (1 for ratio, 0 for difference)
  // Returns true if the deviation is >= 1 sigma
  assert(method == "ratio" || method == "difference" && "method must be 'ratio' or 'difference'");

  cout << avg << " +/- " << error << endl;
  double sig;
  if(method == "difference")
    sig = std::abs(avg)/error;
  else if(method == "ratio")
    sig = std::abs(1.-avg)/error;
  cout << "Average " << method << " is " << sig << " sigma away from 1." << endl;
  if (sig < 1){
    cout << "=> Variation is not significant." << endl;
    return false;
  } else {
    cout << "=> Variation is significant!" << endl;
    return true;
  }
  // maybe do some systematics here based on significance?
}
double computeWeightedAverage(TH1D* hist, double &error, int minBin = 1, int maxBin = -1) {
  // Compute the weighted average of the bin values in range [minBin, maxBin] weighted by their errors
  // Uses inverse variance weighting: weight_i = 1 / error_i^2
  // Default: use all bins (maxBin = -1 means use GetNbinsX())
  if (minBin < 1) minBin = 1;
  if (maxBin < 0) maxBin = hist->GetNbinsX();
  if (maxBin > hist->GetNbinsX()) maxBin = hist->GetNbinsX();
  
  double sumWeights = 0.0;
  double sumWeightedValues = 0.0;

  for (Int_t i = minBin; i <= maxBin; i++) {
    double binContent = hist->GetBinContent(i);
    double binError = hist->GetBinError(i);

    if (binError > 0) {
      double weight = 1.0 / (binError * binError);
      sumWeights += weight;
      sumWeightedValues += weight * binContent;
    }
  }

  if (sumWeights <= 0) {
    error = 0.0;
    return 0.0;
  }

  double weightedAverage = sumWeightedValues / sumWeights;
  error = 1.0 / TMath::Sqrt(sumWeights);

  return weightedAverage;
}

double computeAverageRatio(TH1D* hist1, TH1D* hist2, double &error, int minBin = 1, int maxBin = -1) {
  // Compute the ratio of sums of two histograms using IntegralAndError
  // Default: use all bins (maxBin = -1 means use GetNbinsX())
  if (minBin < 1) minBin = 1;
  if (maxBin < 0) maxBin = hist1->GetNbinsX();
  if (maxBin > hist1->GetNbinsX()) maxBin = hist1->GetNbinsX();
  
  double avg1, avg2, error1, error2;

  avg1 = hist1->IntegralAndError(minBin, maxBin, error1);
  avg2 = hist2->IntegralAndError(minBin, maxBin, error2);

  if (TMath::Abs(avg2) < 1e-12) {
    error = 0.0;
    return 0.0;
  }

  double ratio = avg1 / avg2;
  // Roger barlow error on subset
  error = ratio * TMath::Sqrt(TMath::Abs(error1*error1 - error2*error2));
  error /= avg2;

  return ratio;
}

double computeAverageDifference(TH1D* hist1, TH1D* hist2, double &error, int minBin = 1, int maxBin = -1) {
  // Compute the difference of the averages of two histograms using IntegralAndError
  // Default: use all bins (maxBin = -1 means use GetNbinsX())
  if (minBin < 1) minBin = 1;
  if (maxBin < 0) maxBin = hist1->GetNbinsX();
  if (maxBin > hist1->GetNbinsX()) maxBin = hist1->GetNbinsX();
  
  int nBins = maxBin - minBin + 1;
  if (nBins <= 0) {
    error = 0.0;
    return 0.0;
  }
  
  double avg1, avg2, error1, error2;

  avg1 = hist1->IntegralAndError(minBin, maxBin, error1);
  avg1 /= nBins;
  error1 /= nBins;

  avg2 = hist2->IntegralAndError(minBin, maxBin, error2);
  avg2 /= nBins;
  error2 /= nBins;

  double diff = avg1 - avg2;
  // Error propagation for subtraction: dD^2 = dA1^2 + dA2^2
  error = TMath::Sqrt(TMath::Abs(error1 * error1 - error2 * error2));

  return diff;
}


/**
 * Creates a canvas with two vertically stacked pads for ratio plots.
 * 
 * The canvas is divided into:
 * - Upper pad (70% of canvas height): for main plot
 * - Lower pad (30% of canvas height): for ratio plot
 * 
 * @param name The name identifier for the canvas
 * @param title The title displayed in the canvas window
 * 
 * @return Pointer to the created TCanvas object
 * 
 * @note To access the pads after creation:
 *       - Use canvas->GetPrimitive("padUpper") to get the upper pad
 *       - Use canvas->GetPrimitive("padLower") to get the lower pad
 *       - Alternatively, call canvas->cd(1) for upper pad and canvas->cd(2) for lower pad
 *       - Draw histograms on each pad using pad->cd() followed by histogram->Draw()
 * 
 * @warning The caller is responsible for deleting the returned canvas pointer
 */
TCanvas *plotWithRatio(const char *name, TH1 *hist1, TH1 *hist2, TH1 *ratio, 
                       const char *sysDrawOpt1 = "5", const char *sysDrawOpt2 = "5", 
                       const char *ratioDrawOpt = "E", TLegend *legend = nullptr,
                       int canvasWidth = 800, int canvasHeight = 800,
                       TGraphErrors *sysErrors1 = nullptr, TGraphErrors *sysErrors2 = nullptr,
                       bool isDifference = false) {
  // Pad split: upper takes upperFrac, lower takes (1 - upperFrac)
  const Double_t upperFrac = 0.7;
  const Double_t lowerFrac = 1.0 - upperFrac;
  // Sizes in ROOT are relative to pad height, so to get the same visual size
  // in the lower pad its sizes need to be scaled up by this factor:
  const Double_t padScale = upperFrac / lowerFrac;

  // Clone histograms to avoid modifying originals
  TH1 *h1 = (TH1*)hist1->Clone();
  h1->SetDirectory(nullptr);
  h1->SetBit(kCanDelete);
  TH1 *h2 = nullptr;
  if (hist2) {
    h2 = (TH1*)hist2->Clone();
    h2->SetDirectory(nullptr);
    h2->SetBit(kCanDelete);
  }
  
  // Create canvas and pads
  TCanvas *canvas = new TCanvas(name, name, canvasWidth, canvasHeight);
  
  TPad *padUpper = new TPad("padUpper", "Upper pad", 0, lowerFrac, 1, 1);
  myPadSetUp(padUpper, 0.14, 0.08, 0.04, 0.02);
  padUpper->SetTickx(1); padUpper->SetTicky(1);
  padUpper->Draw();
  
  TPad *padLower = new TPad("padLower", "Lower pad", 0, 0, 1, lowerFrac);
  myPadSetUp(padLower, 0.14, 0.02, 0.04, 0.25);
  padLower->SetTickx(1); padLower->SetTicky(1);
  padLower->Draw();
  
  // --- Upper pad ---
  padUpper->cd();
  h1->SetTitle("");
  h1->GetXaxis()->SetLabelSize(0);
  h1->GetXaxis()->SetTitleSize(0);
  h1->GetXaxis()->SetTickLength(0.03);
  h1->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("y") / upperFrac);
  h1->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("y") / upperFrac);
  h1->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("y") * upperFrac);
  h1->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("y"));
  // h1->SetMarkerStyle(0);
  h1->Draw("E");
  if (h2) {
    h2->SetTitle("");
    // h2->SetMarkerStyle(0);
    h2->Draw("E SAME");
  }
  
  // Helper lambda to style a TGraphErrors consistent with its associated histogram
  auto styleGraph = [](TGraphErrors *g, TH1 *h, bool fill, const char *opt) {
    if (fill) {
      g->SetFillStyle(1001);
      g->SetFillColorAlpha(h->GetLineColor(), 0.3);
      g->SetLineColor(h->GetLineColor());
      g->SetLineWidth(1);
    } else {
      g->SetFillStyle(0);
      g->SetLineColor(h->GetLineColor());
      g->SetLineWidth(2);
    }
    g->Draw(opt);
  };

  bool fillSysErrors1 = TString(sysDrawOpt1).Contains("2") || TString(sysDrawOpt1).Contains("3");
  bool fillSysErrors2 = TString(sysDrawOpt2).Contains("2") || TString(sysDrawOpt2).Contains("3");

  if (sysErrors1) {
    sysErrors1->SetBit(kCanDelete);
    styleGraph(sysErrors1, h1, fillSysErrors1, sysDrawOpt1);
  }
  if (sysErrors2 && h2) {
    sysErrors2->SetBit(kCanDelete);
    styleGraph(sysErrors2, h2, fillSysErrors2, sysDrawOpt2);
  }

  if (legend) {
    legend->SetTextSize(gStyle->GetLegendTextSize() / upperFrac);
    legend->Draw();
  }
  
  // --- Lower pad ---
  padLower->cd();
  TH1 *ratioCopy = (TH1*)ratio->Clone();
  ratioCopy->SetDirectory(nullptr);
  ratioCopy->SetBit(kCanDelete);
  
  ratioCopy->SetTitle("");
  ratioCopy->SetLineColor(kBlack);
  ratioCopy->SetLineWidth(h1->GetLineWidth());
  ratioCopy->SetLineStyle(h1->GetLineStyle());
  ratioCopy->SetMarkerColor(kBlack);
  ratioCopy->SetMarkerStyle(h1->GetMarkerStyle());
  ratioCopy->SetMarkerSize(h1->GetMarkerSize());  // marker size is in absolute units, no scaling needed
  ratioCopy->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());

  // Scale all text sizes by padScale so they match the upper pad visually
  ratioCopy->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("x") / lowerFrac);
  ratioCopy->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("x") / lowerFrac);
  ratioCopy->GetXaxis()->SetTitleOffset(0.9);
  ratioCopy->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("x"));
  ratioCopy->GetXaxis()->SetTickLength(0.03 * padScale);
  ratioCopy->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("y") / lowerFrac);
  ratioCopy->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("y") / lowerFrac);
  ratioCopy->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("y") * lowerFrac);
  ratioCopy->GetYaxis()->CenterTitle();
  ratioCopy->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("y"));
  ratioCopy->GetYaxis()->SetNdivisions(505);
  ratioCopy->GetYaxis()->SetTickLength(0.03);
  // ratioCopy->GetYaxis()->SetNoExponent(true);
  
  // Auto y-axis range: symmetric around reference value (1 for ratio, 0 for difference)
  const Double_t ref = isDifference ? 0.0 : 1.0;

  Double_t minVal = 1e10, maxVal = -1e10;
  for (Int_t i = 1; i <= ratioCopy->GetNbinsX(); i++) {
    Double_t c = ratioCopy->GetBinContent(i);
    Double_t e = ratioCopy->GetBinError(i);
    if (c != 0 || e != 0) { minVal = TMath::Min(minVal, c - e); maxVal = TMath::Max(maxVal, c + e); }
  }
  if (minVal > maxVal) { minVal = ref - 0.1; maxVal = ref + 0.1; }  // fallback if all zero
  Double_t maxDev = std::max({std::abs(maxVal - ref), std::abs(minVal - ref), isDifference ? 1e-8 : 0.05});
  maxDev *= 1.3;  // 30% padding

  // For differences, rescale to avoid tiny labels; encode scale in axis title
  if (isDifference && maxDev > 0) {
    int exp = (int)std::floor(std::log10(maxDev));
    if (exp < -1) {
      double scale = std::pow(10.0, -exp);
      for (int i = 1; i <= ratioCopy->GetNbinsX(); i++) {
        ratioCopy->SetBinContent(i, ratioCopy->GetBinContent(i) * scale);
        ratioCopy->SetBinError(i, ratioCopy->GetBinError(i) * scale);
      }
      maxDev *= scale;
      ratioCopy->GetYaxis()->SetTitle(Form("%s (#times10^{%d})", ratioCopy->GetYaxis()->GetTitle(), exp));
    }
  }

  ratioCopy->GetYaxis()->SetRangeUser(ref - maxDev, ref + maxDev);
  
  ratioCopy->Draw(ratioDrawOpt);
  
  TLine *line = new TLine(ratioCopy->GetXaxis()->GetXmin(), ref, ratioCopy->GetXaxis()->GetXmax(), ref);
  line->SetLineStyle(2); line->SetLineColor(kGray+2); line->SetLineWidth(2);
  line->Draw("SAME");
  
  // Systematic error bands in ratio pad
  if (sysErrors1 && hist2) {
    TGraphErrors *sysRatio1 = new TGraphErrors();
    sysRatio1->SetBit(kCanDelete);
    for (int i = 0; i < sysErrors1->GetN(); i++) {
      double x, y, ex, ey;
      sysErrors1->GetPoint(i, x, y); ex = sysErrors1->GetErrorX(i); ey = sysErrors1->GetErrorY(i);
      int bin = hist2->FindBin(x);
      double denom = hist2->GetBinContent(bin);
      if (denom > 0) sysRatio1->SetPoint(i, x, y / denom), sysRatio1->SetPointError(i, ex, ey / denom);
    }
    sysRatio1->SetLineColor(kBlack); sysRatio1->SetFillStyle(0); sysRatio1->SetLineWidth(2);
    sysRatio1->Draw("SAME 5");
  }
  if (sysErrors2 && hist2) {
    TGraphErrors *sysRatio2 = new TGraphErrors();
    sysRatio2->SetBit(kCanDelete);
    for (int i = 0; i < sysErrors2->GetN(); i++) {
      double x, y, ex, ey;
      sysErrors2->GetPoint(i, x, y); ex = sysErrors2->GetErrorX(i); ey = sysErrors2->GetErrorY(i);
      if (y > 0) sysRatio2->SetPoint(i, x, 1.0), sysRatio2->SetPointError(i, ex, ey / y);
    }
    sysRatio2->SetFillColorAlpha(kRed, 0.3); sysRatio2->SetFillStyle(1001);
    sysRatio2->Draw("SAME 3");
  }
  
  return canvas;
}

