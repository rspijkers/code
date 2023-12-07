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

int postprocessingResults(TString trainnr, TString filename = "AnalysisResults.root") {
  TH1::SetDefaultSumw2(); // Make sure we propagate the errors

  TFile *infile = new TFile("results/" + trainnr + "/" + filename, "READ");
  TDirectory *dir;
  infile->GetObject("cascade-correlations", dir);

  // TODO: get info like Nevents from other dirs/plots in inputfile

  // ROOT is being a bitch about TBrowser in batch mode, so make an outfile:
  TFile *outfile = new TFile("plots/" + trainnr + ".root", "RECREATE");

  // plot a QA plot just to check
  // TH1F *hPhi1;
  // dir->GetObject("hPhi", hPhi1);
  // hPhi1->Draw();
  TH1F *hPhi = dir->Get<TH1F>("hPhi");
  hPhi->SetDirectory(outfile);
  hPhi->Draw();

  // QA plots
  // TODO: inv mass of everything, in the pT bins
  TH2F *hMassXiMinus = dir->Get<TH2F>("hMassXiMinus");
  hMassXiMinus->SetDirectory(outfile);
  TH1D* h1DXiMass = hMassXiMinus->ProjectionX();
  h1DXiMass->Draw();
  TH2F *hMassOmegaMinus = dir->Get<TH2F>("hMassOmegaMinus");
  hMassOmegaMinus->SetDirectory(outfile);
  TH1D* h1DOmegaMass = hMassOmegaMinus->ProjectionX();
  h1DOmegaMass->Draw();

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
  hInvMassXi->Draw();
  TH1D *hInvMassOm = project(hXiXiOS, invMassTrigg, massOm);
  hInvMassOm->SetName("hInvMassOm");
  hInvMassOm->SetTitle("inv Mass of Xi (bachelor == kaon and != pion)");
  hInvMassOm->Draw();
  TH1D *hInvMassBoth = project(hXiXiOS, invMassTrigg, massBoth);
  hInvMassBoth->SetName("hInvMassBoth");
  hInvMassBoth->SetTitle("inv Mass of Xi (bachelor consistent with both pion, kaon)");
  hInvMassBoth->Draw();

  // TODO pT spectra inv mass afhankelijk

  // let's see
  // loop over different correlations (XiXi, XiOm, etc. = different THnSparses)
  // loop over pT --> 3 ranges
  // do sidebands and signal correlations (9 combinations?)
  
  // Let's define some global variables:
  double pTmin = 0.15;
  double pTmax = 15.0;
  double pTbins[4] = {pTmin, 4.0, 8.0, pTmax}; // edges of the different (trigger) pT ranges
  for (int i = 0; i < 4 - 1; i++){
    axranges a{{ptTrigg, {pTbins[i], pTbins[i + 1]}}, {ptAssoc, {pTmin, pTmax}}};
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
  test->Draw();

  // Xi invMass for different pT:
  axranges massLowpT{{ptTrigg, {0., 4.0}}, {selflagTrigg, {0.5, 1.5}}};
  axranges massMedpT{{ptTrigg, {4.0, 8.0}}, {selflagTrigg, {0.5, 1.5}}};
  axranges massHigpT{{ptTrigg, {8.0, 15.0}}, {selflagTrigg, {0.5, 1.5}}};
  TH1D *hInvMassXiLow = project(hXiXiOS, invMassTrigg, massLowpT);
  hInvMassXiLow->SetName("hInvMassXiLow");
  hInvMassXiLow->SetTitle("inv Mass of Xi trigger (w PID response, 0 < pT < 4)");
  hInvMassXiLow->Draw();
  TH1D *hInvMassXiMed = project(hXiXiOS, invMassTrigg, massMedpT);
  hInvMassXiMed->SetName("hInvMassXiMed");
  hInvMassXiMed->SetTitle("inv Mass of Xi trigger (w PID response, 4 < pT < 8)");
  hInvMassXiMed->Draw();
  TH1D *hInvMassXiHig = project(hXiXiOS, invMassTrigg, massHigpT);
  hInvMassXiHig->SetName("hInvMassXiHig");
  hInvMassXiHig->SetTitle("inv Mass of Xi trigger (w PID response, 8 < pT < 15)");
  hInvMassXiHig->Draw();

  // OS
  TH1D *hdphiLowOS = project(hXiXiOS, dPhi, xiLow);
  hdphiLowOS->SetName("hdphiLowOS");
  hdphiLowOS->Draw();
  TH1D *hdphiMedOS = project(hXiXiOS, dPhi, xiMed);
  hdphiMedOS->SetName("hdphiMedOS");
  hdphiMedOS->Draw();
  TH1D *hdphiHigOS = project(hXiXiOS, dPhi, xiHig);
  hdphiHigOS->SetName("hdphiHigOS");
  hdphiHigOS->Draw();

  // TH1D *hdphiXiOmOS = project(hXiXiOS, dPhi, xiOmega);
  // hdphiXiOmOS->SetName("hdphiXiOmOS");
  // hdphiXiOmOS->Draw();
  // TH1D *hdphiOmXiOS = project(hXiXiOS, dPhi, omegaXi);
  // hdphiOmXiOS->SetName("hdphiOmXiOS");
  // hdphiOmXiOS->Draw();
  // TH1D *hdphiOmOmOS = project(hXiXiOS, dPhi, omegaOmega);
  // hdphiOmOmOS->SetName("hdphiOmOmOS");
  // hdphiOmOmOS->Draw();

  // SS
  TH1D *hdphiLowSS = project(hXiXiSS, dPhi, xiLow);
  hdphiLowSS->SetName("hdphiLowSS");
  hdphiLowSS->Draw();
  TH1D *hdphiMedSS = project(hXiXiSS, dPhi, xiMed);
  hdphiMedSS->SetName("hdphiMedSS");
  hdphiMedSS->Draw();
  TH1D *hdphiHigSS = project(hXiXiSS, dPhi, xiHig);
  hdphiHigSS->SetName("hdphiHigSS");
  hdphiHigSS->Draw();

  // TH1D *hdphiXiOmSS = project(hXiXiSS, dPhi, xiOmega);
  // hdphiXiOmSS->SetName("hdphiXiOmSS");
  // hdphiXiOmSS->Draw();
  // TH1D *hdphiOmXiSS = project(hXiXiSS, dPhi, omegaXi);
  // hdphiOmXiSS->SetName("hdphiOmXiSS");
  // hdphiOmXiSS->Draw();
  // TH1D *hdphiOmOmSS = project(hXiXiSS, dPhi, omegaOmega);
  // hdphiOmOmSS->SetName("hdphiOmOmSS");
  // hdphiOmOmSS->Draw();

  // OS - SS
  TH1D *hdphiLow = new TH1D(*hdphiLowOS);
  hdphiLow->Add(hdphiLowOS, hdphiLowSS, 1, -1);
  hdphiLow->SetName("hdphiLow");
  hdphiLow->Draw();
  TH1D *hdphiMed = new TH1D(*hdphiMedOS);
  hdphiMed->Add(hdphiMedOS, hdphiMedSS, 1, -1);
  hdphiMed->SetName("hdphiMed");
  hdphiMed->Draw();
  TH1D *hdphiHig = new TH1D(*hdphiHigOS);
  hdphiHig->Add(hdphiHigOS, hdphiHigSS, 1, -1);
  hdphiHig->SetName("hdphiHig");
  hdphiHig->Draw();

  // TH1D *hdphiXiOm = new TH1D(*hdphiXiOmOS);
  // hdphiXiOm->Add(hdphiXiOmOS, hdphiXiOmSS, 1, -1);
  // hdphiXiOm->SetName("hdphiXiOm");
  // hdphiXiOm->Draw();
  // TH1D *hdphiOmXi = new TH1D(*hdphiOmXiOS);
  // hdphiOmXi->Add(hdphiOmXiOS, hdphiOmXiSS, 1, -1);
  // hdphiOmXi->SetName("hdphiOmXi");
  // hdphiOmXi->Draw();
  // TH1D *hdphiOmOm = new TH1D(*hdphiOmOmOS);
  // hdphiOmOm->Add(hdphiOmOmOS, hdphiOmOmSS, 1, -1);
  // hdphiOmOm->SetName("hdphiOmOm");
  // hdphiOmOm->Draw();

  // Now do the same with inv mass selection for signal:
  axranges xiLowInvMass{{ptTrigg, {0., 4.0}}, {ptAssoc, {0., 15.0}}, {invMassTrigg, {1.31, 1.335}}, {invMassAssoc, {1.31, 1.335}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiMedInvMass{{ptTrigg, {4.0, 8.0}},  {ptAssoc, {0.15, 15.0}}, {invMassTrigg, {1.31, 1.335}}, {invMassAssoc, {1.31, 1.335}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiHigInvMass{{ptTrigg, {8.0, 15.0}}, {ptAssoc, {0.15, 15.0}}, {invMassTrigg, {1.31, 1.335}}, {invMassAssoc, {1.31, 1.335}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};

  // OS
  TH1D *hdphiLowInvOS = project(hXiXiOS, dPhi, xiLowInvMass);
  hdphiLowInvOS->SetName("hdphiLowInvOS");
  hdphiLowInvOS->Draw();
  TH1D *hdphiMedInvOS = project(hXiXiOS, dPhi, xiMedInvMass);
  hdphiMedInvOS->SetName("hdphiMedInvOS");
  hdphiMedInvOS->Draw();
  TH1D *hdphiHigInvOS = project(hXiXiOS, dPhi, xiHigInvMass);
  hdphiHigInvOS->SetName("hdphiHigInvOS");
  hdphiHigInvOS->Draw();

  // SS
  TH1D *hdphiLowInvSS = project(hXiXiSS, dPhi, xiLowInvMass);
  hdphiLowInvSS->SetName("hdphiLowInvSS");
  hdphiLowInvSS->Draw();
  TH1D *hdphiMedInvSS = project(hXiXiSS, dPhi, xiMedInvMass);
  hdphiMedInvSS->SetName("hdphiMedInvSS");
  hdphiMedInvSS->Draw();
  TH1D *hdphiHigInvSS = project(hXiXiSS, dPhi, xiHigInvMass);
  hdphiHigInvSS->SetName("hdphiHigInvSS");
  hdphiHigInvSS->Draw();

  // OS - SS
  TH1D *hdphiLowInv = new TH1D(*hdphiLowOS);
  hdphiLowInv->Add(hdphiLowInvOS, hdphiLowInvSS, 1, -1);
  hdphiLowInv->SetName("hdphiLowInv");
  hdphiLowInv->Draw();
  TH1D *hdphiMedInv = new TH1D(*hdphiMedOS);
  hdphiMedInv->Add(hdphiMedInvOS, hdphiMedInvSS, 1, -1);
  hdphiMedInv->SetName("hdphiMedInv");
  hdphiMedInv->Draw();
  TH1D *hdphiHigInv = new TH1D(*hdphiHigOS);
  hdphiHigInv->Add(hdphiHigInvOS, hdphiHigInvSS, 1, -1);
  hdphiHigInv->SetName("hdphiHigInv");
  hdphiHigInv->Draw();

  outfile->Write();
  outfile->Close();
  return 0;
}