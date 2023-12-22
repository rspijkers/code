// This script makes the relevant plots/projections from a given AnalysisResults.root.
// Run me like `root 'postprocessingResults.cxx("trainnr", "outfilename")' -q`, where the second 
//  variable is optional in case the filename diverges from "AnalysisResults.root".

// std
#include <iostream>
#include <vector>
#include <set>
#include <map>
// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TString.h"

using std::cout; using std::endl;

// enum for axisnumbers/names, in- or outside main()?
enum {
  dPhi, 
  dEta, 
  ptTrigg, 
  ptAssoc, 
  invMassTrigg, 
  invMassAssoc, 
  selflagTrigg, 
  selflagAssoc, 
  V_z, 
  multiplicity
}; // for more info on the THnSparse see O2 task: PWGLF/Tasks/cascadecorrelations.cxx

// Wrapper for lower and upper bounds of an axis
using axranges = std::map<int, std::vector<double>>;

// Function that creates a projection to TargetAxis, with ranges in other dimensions
TH1D *project(THnSparse *THn,             // input THn
              int targetAxis,             // axis nr on which to project
              axranges map,               // map with [axisnr, {lower, upper}] bounds
              Option_t *options = "") {   // any options to pass to THnSparse->Projection()
  int ndim = THn->GetNdimensions();
  for(auto const& [axis, bounds] : map){
    TAxis *ax = THn->GetAxis(axis);
    if (!ax) cout << "Error: GetAxis() with number " << axis << " gives a null pointer! The program will crash!" << endl; // null pointer check
    ax->SetRangeUser(bounds[0], bounds[1]);
  }
  TH1D *hp = THn->Projection(targetAxis, options);
  TString axisname = THn->GetAxis(targetAxis)->GetTitle();
  hp->SetTitle("Projection on " + axisname);

  // reset axis ranges
  for(int i=0; i<ndim; i++) {
    TAxis *ax = THn->GetAxis(i);
    ax->SetRange(0,0);
  }
  return hp;
}

// define function that projects all the QA histograms?

int postprocessingResults(TString trainnr, TString filename = "AnalysisResults.root") {
  TH1::SetDefaultSumw2(); // Make sure we propagate the errors

  TFile *infile = new TFile("results/" + trainnr + "/" + filename, "READ");
  TDirectory *dir;
  infile->GetObject("cascade-correlations", dir);

  // ROOT is being a bitch about TBrowser in batch mode, so make an outfile:
  TFile *outfile = new TFile("plots/" + trainnr + ".root", "RECREATE");

  // Let's define some global variables:
  double pTmin = 0.15;
  double pTmax = 15.0;
  const int maxPtBins = 10;
  double pTbins[maxPtBins] = {pTmin, 1.0, 1.9, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, pTmax}; // edges of the different (trigger) pT ranges
  // const TString pTlabels[3] = {"Low", "Med", "Hig"};
  const TString pTlabels[maxPtBins - 1] = {"1", "1.9", "2", "3", "4", "5", "6", "8", "max"};

  // QA plots
  TH1F *hPhi = dir->Get<TH1F>("hPhi");
  hPhi->SetDirectory(outfile);
  hPhi->Draw();

  TF1 *f1 = new TF1("f1", "pol2(0) + gaus(3) + gaus(6)", 1.29, 1.42); // pol2 doesn't work, maybe normalize histo first? --> discrepancy between y-axis O(10^6) x-axis O(0.1)
  f1->SetParameters(0, 0, 0, 0, 1.32, 0.02, 0, 1.32, 0.005); // set the pol parameters to 0 (less important than gaussian), determine the constant of the gaussian in the loop
  // do inv mass fit in pT bins
  TH2F *hMassXiMinus = dir->Get<TH2F>("hMassXiMinus");
  hMassXiMinus->SetDirectory(outfile);
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hMassXiMinus->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin+1]);
    TH1D *h = hMassXiMinus->ProjectionX("hMassXiMinus_"+ pTlabels[pTbin]);
    // Double_t nentries = h->Integral();
    // h->Scale(1. / nentries, "width");
    f1->SetParameter(3, h->GetMaximum()/2);
    f1->SetParameter(6, h->GetMaximum()/2);
    // f1->SetParameter(0, h->GetBinContent(h->FindBin(1.4)));
    TFitResultPtr r = h->Fit("f1", "LQ", "", 1.29, 1.42);
    std::cout << f1->GetParameter(0) << " and " << f1->GetParameter(3) << std::endl;
  }

  // TH1D* h1DXiMass = hMassXiMinus->ProjectionX();
  TH2F *hMassOmegaMinus = dir->Get<TH2F>("hMassOmegaMinus");
  hMassOmegaMinus->SetDirectory(outfile);
  TH1D* h1DOmegaMass = hMassOmegaMinus->ProjectionX();

  // Load the THnSparses, 8 in total (4 particle combinations * 2 sign combinations)
  THnSparse *hXiXiOS, *hXiXiSS;
  dir->GetObject("hXiXiOS", hXiXiOS);
  dir->GetObject("hXiXiSS", hXiXiSS);

  // Let's do inv-mass projections, only Xi mass hypothesis for now
  axranges massXi{{selflagTrigg, {0.5, 1.5}}};
  axranges massOm{{selflagTrigg, {1.5, 2.5}}};
  axranges massBoth{{selflagTrigg, {2.5, 3.5}}};
  TH1D *hInvMassXi = project(hXiXiOS, invMassTrigg, massXi);
  hInvMassXi->SetName("hInvMassXi");
  hInvMassXi->SetTitle("inv Mass of Xi (bachelor == pion and != kaon)");
  TH1D *hInvMassOm = project(hXiXiOS, invMassTrigg, massOm);
  hInvMassOm->SetName("hInvMassOm");
  hInvMassOm->SetTitle("inv Mass of Xi (bachelor == kaon and != pion)");
  TH1D *hInvMassBoth = project(hXiXiOS, invMassTrigg, massBoth);
  hInvMassBoth->SetName("hInvMassBoth");
  hInvMassBoth->SetTitle("inv Mass of Xi (bachelor consistent with both pion, kaon)");

  // TODO pT spectra inv mass afhankelijk

  // let's see
  // loop over different correlations (XiXi, XiOm, etc. = different THnSparses)
  // loop over pT --> 3 ranges
  // do sidebands and signal correlations (9 combinations?)
  for (int i = 0; i < maxPtBins-1; i++){
    axranges a{{ptTrigg, {pTbins[i], pTbins[i + 1]}}, {ptAssoc, {pTmin, pTmax}}};
    TH1D *h = project(hXiXiOS, invMassTrigg, a);
    h->SetName("hInvMassXi"+ pTlabels[i]);
    f1->SetParameter(3, h->GetMaximum());
    TFitResultPtr r = h->Fit("f1", "LQ", "", 1.29, 1.42);
    
    // Okay, now sidebands/signal loop? How to do this? 
    // so for signal define a mass window I guess, and for sidebands as well. 
    // instead of nsigma (sigma = width of peak), maybe use bare nrs?
  }

  // projection configurations
  axranges xiLow{{ptTrigg, {0., 4.0}}, {ptAssoc, {0., 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiMed{{ptTrigg, {4.0, 8.0}},  {ptAssoc, {0.15, 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiHig{{ptTrigg, {8.0, 15.0}}, {ptAssoc, {0.15, 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  // axranges xiOmega{{selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {1.5, 2.5}}};
  // axranges omegaXi{{selflagTrigg, {1.5, 2.5}}, {selflagAssoc, {0.5, 1.5}}};
  // axranges omegaOmega{{selflagTrigg, {1.5, 2.5}}, {selflagAssoc, {1.5, 2.5}}};

  // test
  axranges at{{selflagAssoc, {2.5, 3.5}}};
  TH1D *test = project(hXiXiOS, dPhi, at);
  test->SetName("test");

  // Xi invMass for different pT:
  axranges massLowpT{{ptTrigg, {0., 4.0}}, {selflagTrigg, {0.5, 1.5}}};
  axranges massMedpT{{ptTrigg, {4.0, 8.0}}, {selflagTrigg, {0.5, 1.5}}};
  axranges massHigpT{{ptTrigg, {8.0, 15.0}}, {selflagTrigg, {0.5, 1.5}}};
  TH1D *hInvMassXiLow = project(hXiXiOS, invMassTrigg, massLowpT);
  hInvMassXiLow->SetName("hInvMassXiLow");
  hInvMassXiLow->SetTitle("inv Mass of Xi trigger (w PID response, 0 < pT < 4)");
  TH1D *hInvMassXiMed = project(hXiXiOS, invMassTrigg, massMedpT);
  hInvMassXiMed->SetName("hInvMassXiMed");
  hInvMassXiMed->SetTitle("inv Mass of Xi trigger (w PID response, 4 < pT < 8)");
  TH1D *hInvMassXiHig = project(hXiXiOS, invMassTrigg, massHigpT);
  hInvMassXiHig->SetName("hInvMassXiHig");
  hInvMassXiHig->SetTitle("inv Mass of Xi trigger (w PID response, 8 < pT < 15)");

  // OS
  TH1D *hdphiLowOS = project(hXiXiOS, dPhi, xiLow);
  hdphiLowOS->SetName("hdphiLowOS");
  TH1D *hdphiMedOS = project(hXiXiOS, dPhi, xiMed);
  hdphiMedOS->SetName("hdphiMedOS");
  TH1D *hdphiHigOS = project(hXiXiOS, dPhi, xiHig);
  hdphiHigOS->SetName("hdphiHigOS");

  // TH1D *hdphiXiOmOS = project(hXiXiOS, dPhi, xiOmega);
  // hdphiXiOmOS->SetName("hdphiXiOmOS");
  // TH1D *hdphiOmXiOS = project(hXiXiOS, dPhi, omegaXi);
  // hdphiOmXiOS->SetName("hdphiOmXiOS");
  // TH1D *hdphiOmOmOS = project(hXiXiOS, dPhi, omegaOmega);
  // hdphiOmOmOS->SetName("hdphiOmOmOS");

  // SS
  TH1D *hdphiLowSS = project(hXiXiSS, dPhi, xiLow);
  hdphiLowSS->SetName("hdphiLowSS");
  TH1D *hdphiMedSS = project(hXiXiSS, dPhi, xiMed);
  hdphiMedSS->SetName("hdphiMedSS");
  TH1D *hdphiHigSS = project(hXiXiSS, dPhi, xiHig);
  hdphiHigSS->SetName("hdphiHigSS");

  // TH1D *hdphiXiOmSS = project(hXiXiSS, dPhi, xiOmega);
  // hdphiXiOmSS->SetName("hdphiXiOmSS");
  // TH1D *hdphiOmXiSS = project(hXiXiSS, dPhi, omegaXi);
  // hdphiOmXiSS->SetName("hdphiOmXiSS");
  // TH1D *hdphiOmOmSS = project(hXiXiSS, dPhi, omegaOmega);
  // hdphiOmOmSS->SetName("hdphiOmOmSS");

  // OS - SS
  TH1D *hdphiLow = new TH1D(*hdphiLowOS);
  hdphiLow->Add(hdphiLowOS, hdphiLowSS, 1, -1);
  hdphiLow->SetName("hdphiLow");
  TH1D *hdphiMed = new TH1D(*hdphiMedOS);
  hdphiMed->Add(hdphiMedOS, hdphiMedSS, 1, -1);
  hdphiMed->SetName("hdphiMed");
  TH1D *hdphiHig = new TH1D(*hdphiHigOS);
  hdphiHig->Add(hdphiHigOS, hdphiHigSS, 1, -1);
  hdphiHig->SetName("hdphiHig");

  // TH1D *hdphiXiOm = new TH1D(*hdphiXiOmOS);
  // hdphiXiOm->Add(hdphiXiOmOS, hdphiXiOmSS, 1, -1);
  // hdphiXiOm->SetName("hdphiXiOm");
  // TH1D *hdphiOmXi = new TH1D(*hdphiOmXiOS);
  // hdphiOmXi->Add(hdphiOmXiOS, hdphiOmXiSS, 1, -1);
  // hdphiOmXi->SetName("hdphiOmXi");
  // TH1D *hdphiOmOm = new TH1D(*hdphiOmOmOS);
  // hdphiOmOm->Add(hdphiOmOmOS, hdphiOmOmSS, 1, -1);
  // hdphiOmOm->SetName("hdphiOmOm");

  // Now do the same with inv mass selection for signal:
  axranges xiLowInvMass{{ptTrigg, {0., 4.0}}, {ptAssoc, {0., 15.0}}, {invMassTrigg, {1.31, 1.335}}, {invMassAssoc, {1.31, 1.335}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiMedInvMass{{ptTrigg, {4.0, 8.0}},  {ptAssoc, {0.15, 15.0}}, {invMassTrigg, {1.31, 1.335}}, {invMassAssoc, {1.31, 1.335}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiHigInvMass{{ptTrigg, {8.0, 15.0}}, {ptAssoc, {0.15, 15.0}}, {invMassTrigg, {1.31, 1.335}}, {invMassAssoc, {1.31, 1.335}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};

  // OS
  TH1D *hdphiLowInvOS = project(hXiXiOS, dPhi, xiLowInvMass);
  hdphiLowInvOS->SetName("hdphiLowInvOS");
  TH1D *hdphiMedInvOS = project(hXiXiOS, dPhi, xiMedInvMass);
  hdphiMedInvOS->SetName("hdphiMedInvOS");
  TH1D *hdphiHigInvOS = project(hXiXiOS, dPhi, xiHigInvMass);
  hdphiHigInvOS->SetName("hdphiHigInvOS");

  // SS
  TH1D *hdphiLowInvSS = project(hXiXiSS, dPhi, xiLowInvMass);
  hdphiLowInvSS->SetName("hdphiLowInvSS");
  TH1D *hdphiMedInvSS = project(hXiXiSS, dPhi, xiMedInvMass);
  hdphiMedInvSS->SetName("hdphiMedInvSS");
  TH1D *hdphiHigInvSS = project(hXiXiSS, dPhi, xiHigInvMass);
  hdphiHigInvSS->SetName("hdphiHigInvSS");

  // OS - SS
  TH1D *hdphiLowInv = new TH1D(*hdphiLowOS);
  hdphiLowInv->Add(hdphiLowInvOS, hdphiLowInvSS, 1, -1);
  hdphiLowInv->SetName("hdphiLowInv");
  TH1D *hdphiMedInv = new TH1D(*hdphiMedOS);
  hdphiMedInv->Add(hdphiMedInvOS, hdphiMedInvSS, 1, -1);
  hdphiMedInv->SetName("hdphiMedInv");
  TH1D *hdphiHigInv = new TH1D(*hdphiHigOS);
  hdphiHigInv->Add(hdphiHigInvOS, hdphiHigInvSS, 1, -1);
  hdphiHigInv->SetName("hdphiHigInv");

  outfile->Write();
  outfile->Close();
  return 0;
}