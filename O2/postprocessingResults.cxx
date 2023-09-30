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
  V_z
}; // for more info on the THnSparse see O2 task: PWGLF/Tasks/cascadecorrelations.cxx

// Wrapper for lower and upper bounds
using axranges = std::map<int, std::vector<double>>;

// Function that creates a projection to TargetAxis, with ranges in other dimensions
TH1D *project(THnSparse *THn,             // input THn
              int targetAxis,             // axis nr on which to project
              axranges map,               // map with [axisnr, {lower, upper}] bounds
              Option_t *options = "") {   // any options to pass to THnSparse->Projection()
  int ndim = THn->GetNdimensions();
  for(auto const& [axis, bounds] : map){
    TAxis *ax = THn->GetAxis(axis);
    ax->SetRangeUser(bounds[0], bounds[1]);
  }
  TH1D *hp = THn->Projection(targetAxis, options);

  // reset axis ranges
  for(int i=0; i<ndim; i++) {
    TAxis *ax = THn->GetAxis(i);
    ax->SetRange(0,0);
  }
  return hp;
}

int postprocessingResults(const char *filename) {
  TH1::SetDefaultSumw2(); // Make sure we propagate the errors

  TFile *infile = new TFile(filename);
  TDirectory *dir;
  infile->GetObject("cascade-correlations", dir);

  // ROOT is being a bitch about TBrowser in batch mode, so make an outfile:
  TFile *outfile = new TFile("plots.root", "RECREATE");

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

  // Load the THnSparses
  THnSparse *hSparseSS, *hSparseOS;
  dir->GetObject("hSparseSS", hSparseSS);
  dir->GetObject("hSparseOS", hSparseOS);

  // projection configurations
  axranges xiLow{{ptTrigg, {0.15, 4.0}}, {ptAssoc, {0.15, 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiMed{{ptTrigg, {4.0, 8.0}},  {ptAssoc, {0.15, 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiHig{{ptTrigg, {8.0, 15.0}}, {ptAssoc, {0.15, 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiOmega{{selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {1.5, 2.5}}};
  axranges omegaXi{{selflagTrigg, {1.5, 2.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges omegaOmega{{selflagTrigg, {1.5, 2.5}}, {selflagAssoc, {1.5, 2.5}}};

  // OS
  TH1D *hdphiLowOS = project(hSparseOS, dPhi, xiLow);
  hdphiLowOS->SetName("hdphiLowOS");
  hdphiLowOS->Draw();
  TH1D *hdphiMedOS = project(hSparseOS, dPhi, xiMed);
  hdphiMedOS->SetName("hdphiMedOS");
  hdphiMedOS->Draw();
  TH1D *hdphiHigOS = project(hSparseOS, dPhi, xiHig);
  hdphiHigOS->SetName("hdphiHigOS");
  hdphiHigOS->Draw();
  TH1D *hdphiXiOmOS = project(hSparseOS, dPhi, xiOmega);
  hdphiXiOmOS->SetName("hdphiXiOmOS");
  hdphiXiOmOS->Draw();
  TH1D *hdphiOmXiOS = project(hSparseOS, dPhi, omegaXi);
  hdphiOmXiOS->SetName("hdphiOmXiOS");
  hdphiOmXiOS->Draw();
  TH1D *hdphiOmOmOS = project(hSparseOS, dPhi, omegaOmega);
  hdphiOmOmOS->SetName("hdphiOmOmOS");
  hdphiOmOmOS->Draw();

  // SS
  TH1D *hdphiLowSS = project(hSparseSS, dPhi, xiLow);
  hdphiLowSS->SetName("hdphiLowSS");
  hdphiLowSS->Draw();
  TH1D *hdphiMedSS = project(hSparseSS, dPhi, xiMed);
  hdphiMedSS->SetName("hdphiMedSS");
  hdphiMedSS->Draw();
  TH1D *hdphiHigSS = project(hSparseSS, dPhi, xiHig);
  hdphiHigSS->SetName("hdphiHigSS");
  hdphiHigSS->Draw();
  TH1D *hdphiXiOmSS = project(hSparseSS, dPhi, xiOmega);
  hdphiXiOmSS->SetName("hdphiXiOmSS");
  hdphiXiOmSS->Draw();
  TH1D *hdphiOmXiSS = project(hSparseSS, dPhi, omegaXi);
  hdphiOmXiSS->SetName("hdphiOmXiSS");
  hdphiOmXiSS->Draw();
  TH1D *hdphiOmOmSS = project(hSparseSS, dPhi, omegaOmega);
  hdphiOmOmSS->SetName("hdphiOmOmSS");
  hdphiOmOmSS->Draw();

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
  TH1D *hdphiXiOm = new TH1D(*hdphiXiOmOS);
  hdphiXiOm->Add(hdphiXiOmOS, hdphiXiOmSS, 1, -1);
  hdphiXiOm->SetName("hdphiXiOm");
  hdphiXiOm->Draw();
  TH1D *hdphiOmXi = new TH1D(*hdphiOmXiOS);
  hdphiOmXi->Add(hdphiOmXiOS, hdphiOmXiSS, 1, -1);
  hdphiOmXi->SetName("hdphiOmXi");
  hdphiOmXi->Draw();
  TH1D *hdphiOmOm = new TH1D(*hdphiOmOmOS);
  hdphiOmOm->Add(hdphiOmOmOS, hdphiOmOmSS, 1, -1);
  hdphiOmOm->SetName("hdphiOmOm");
  hdphiOmOm->Draw();

  outfile->Write();
  outfile->Close();
  return 0;
}