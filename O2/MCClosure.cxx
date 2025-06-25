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
#include "TGrid.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRatioPlot.h"

#include "postprocessingTools.h"

// using json = nlohmann::json;
using std::cout; using std::endl;

// another enum for MC closure test:
struct mc{
  enum{
    dPhi, 
    dY, 
    ptTrigg, 
    ptAssoc,
    V_z,
    multiplicity
  }; 
};
// for more info on the THnSparse see O2 task: PWGLF/Tasks/cascadecorrelations.cxx

// Wrapper for lower and upper bounds of an axis
using axranges = std::map<int, std::vector<double>>;

// // pT bins corresponding to efficiency 
// const std::vector<double> pTbins = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0};
// const std::vector<TString> pTlabels = {"0.4", "0.6", "0.8", "1.0", "1.2", "1.4", "1.8", "2.0", "2.2", "2.4", "2.6", "2.8", "3.0", "3.5", "4.0", "4.5", "5.0", "6.0", "8.0", "10.0"};
// assert(pTbins.size() == pTlabels.size() && "pTbins and pTlabels have different sizes, something is wrong!");
// const int maxPtBins = pTbins.size();
// const double pTmin = pTbins[0];
// const double pTmax = pTbins[maxPtBins - 1];

// broader bins:
const std::vector<double> pTbins = {0.6, 1.0, 2.0, 3.0, 5.0, 12.0};
const std::vector<TString> pTlabels = {"0.6", "1.0", "2.0", "3.0", "5.0", "12.0"};
assert(pTbins.size() == pTlabels.size() && "pTbins and pTlabels have different sizes, something is wrong!");
const int maxPtBins = pTbins.size();
const double pTmin = pTbins[0];
const double pTmax = pTbins[maxPtBins - 1];

// vectors/arrays with boundaries of signal/bkg regions of inv mass plots
double sigXi[100][2], bkgXi[100][2]; // can't be variable length (i.e. maxPtBins - 1) so just put it to 100.
double sigOm[100][2], bkgOm[100][2]; // can't be variable length (i.e. maxPtBins - 1) so just put it to 100.

// initialize some variables, to be set in main
TFile *inputFile;
TFile *outputFile;
TDirectory *inputDir;

int MCClosure(TString trainnr, TString filename = "AnalysisResults.root", bool makePDF = false) {
  TH1::SetDefaultSumw2(); // Make sure we propagate the errors
  gStyle->SetHistLineWidth(3);
  gROOT->ForceStyle();

  if(trainnr == "test"){ // if not a trainnr but just a test, look for AnalysisResults.root in current dir
    inputFile = new TFile(filename, "READ");
  } else {
    inputFile = new TFile("results/" + trainnr + "/" + filename, "READ");
  }
  inputFile->GetObject("cascade-correlations", inputDir);
  outputFile = new TFile("plots/Closure" + trainnr + ".root", "RECREATE");
  
  // make sure inputFiles and stuff are not null
  assert((!inputFile && !inputDir && !outputFile) && "input or output file is null pointer!");

  THnSparse *hMCPlusMinus, *hMCPlusPlus, *hMCMinusMinus, *hMCMinusPlus;
  inputDir->GetObject("MC/hMCPlusMinus", hMCPlusMinus);
  inputDir->GetObject("MC/hMCPlusPlus", hMCPlusPlus);
  inputDir->GetObject("MC/hMCMinusMinus", hMCMinusMinus);
  inputDir->GetObject("MC/hMCMinusPlus", hMCMinusPlus);
  cout << hMCMinusMinus->GetAxis(1)->GetXmin() << endl;

  // get normalization from efficiency
  TH2F *hXiMinGen, *hXiPlusGen;
  inputFile->GetObject("cascade-selector/gen/hXiMinus", hXiMinGen);
  inputFile->GetObject("cascade-selector/gen/hXiPlus", hXiPlusGen);
  TH1D *hXiMinGen1D = hXiMinGen->ProjectionX();
  TH1D *hXiPlusGen1D = hXiPlusGen->ProjectionX();
  double nXiMin = hXiMinGen1D->Integral(hXiMinGen1D->FindBin(1.0), hXiMinGen1D->FindBin(10.));
  double nXiPlus = hXiPlusGen1D->Integral(hXiPlusGen1D->FindBin(1.0), hXiPlusGen1D->FindBin(10.));
  cout << "nXiMin = " << nXiMin << ", nXiPlus = " << nXiPlus << endl;
  
  // 1D to start
  axranges aMC{{mc::ptTrigg, {1., 10}}, {mc::ptAssoc, {1., 10}},
               {mc::dY, {-1., 1.}}
              };
  TH1D* hPlusMinus = project(hMCPlusMinus, mc::dPhi, aMC);
  hPlusMinus->SetName("hPlusMinus");
  hPlusMinus->SetTitle("#Xi^{+} - #Xi^{-} correlations");
  hPlusMinus->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  hPlusMinus->Scale(1./nXiPlus);
  hPlusMinus->Rebin(9);
  TH1D* hPlusPlus = project(hMCPlusPlus, mc::dPhi, aMC);
  hPlusPlus->SetName("hPlusPlus");
  hPlusPlus->SetTitle("#Xi^{+} - #Xi^{+} correlations");
  hPlusPlus->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  hPlusPlus->Scale(1./nXiPlus);
  hPlusPlus->Rebin(9);
  TH1D* hMinusPlus = project(hMCMinusPlus, mc::dPhi, aMC);
  hMinusPlus->SetName("hMinusPlus");
  hMinusPlus->SetTitle("#Xi^{-} - #Xi^{+} correlations");
  hMinusPlus->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  hMinusPlus->Scale(1./nXiMin);
  hMinusPlus->Rebin(9);
  TH1D* hMinusMinus = project(hMCMinusMinus, mc::dPhi, aMC);
  hMinusMinus->SetName("hMinusMinus");
  hMinusMinus->SetTitle("#Xi^{-} - #Xi^{-} correlations");
  hMinusMinus->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  hMinusMinus->Scale(1./nXiMin);
  hMinusMinus->Rebin(9);

  // DeltaRapidity
  TH1D* hRapidityMinusPlus = project(hMCMinusPlus, mc::dY, aMC);
  hRapidityMinusPlus->SetName("hRapidityMinusPlus");
  hRapidityMinusPlus->SetTitle("#Xi^{-} - #Xi^{+} correlations #Delta y");
  TH1D* hRapidityMinusMinus = project(hMCMinusMinus, mc::dY, aMC);
  hRapidityMinusMinus->SetName("hRapidityMinusMinus");
  hRapidityMinusMinus->SetTitle("#Xi^{-} - #Xi^{-} correlations #Delta y");
  TH1D* hRapidityPlusMinus = project(hMCPlusMinus, mc::dY, aMC);
  hRapidityPlusMinus->SetName("hRapidityPlusMinus");
  hRapidityPlusMinus->SetTitle("#Xi^{+} - #Xi^{-} correlations #Delta y");
  TH1D* hRapidityPlusPlus = project(hMCPlusPlus, mc::dY, aMC);
  hRapidityPlusPlus->SetName("hRapidityPlusPlus");
  hRapidityPlusPlus->SetTitle("#Xi^{+} - #Xi^{+} correlations #Delta y");

  TH1D* hPlusSubtracted = new TH1D(*hPlusMinus);
  hPlusSubtracted->SetName("hPlusSubtracted");
  hPlusSubtracted->Add(hPlusMinus, hPlusPlus, 1, -1);
  TH1D* hMinusSubtracted = new TH1D(*hPlusMinus);
  hMinusSubtracted->SetName("hMinusSubtracted");
  hMinusSubtracted->Add(hMinusPlus, hMinusMinus, 1, -1);
  TH1D* hAveraged = new TH1D(*hPlusMinus);
  hAveraged->SetName("hAveraged");
  hAveraged->Add(hPlusSubtracted, hMinusSubtracted, .5, .5);

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  return 0;
}