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

// Structure of the THnSparse:
// 0 dPhi, 1 dEta, 2 ptTrigger, 3 ptAssoc, 4 selflagTrigger, 5 selflagAssoc, 6 V_z

// more info:
// AxisSpec invMassAxis = {3000, 0.0f, 3.0f, "Inv. Mass (GeV/c^{2})"};
// AxisSpec deltaPhiAxis = {100, -PI / 2, 1.5  *PI, "#Delta#varphi"};
// AxisSpec deltaEtaAxis = {40, -2, 2, "#Delta#eta"};
// AxisSpec ptAxis = {200, 0, 15, "#it{p}_{T}"};
// AxisSpec selectionFlagAxis = {4, -0.05f, 3.5f, "Selection flag of casc candidate"};
// AxisSpec vertexAxis = {1000, -10.0f, 10.0f, "cm"};

// TODO: make project function that takes THnSparse, target axis, and range of any other axis/axes

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

  // TODO FIX THE SELECTIONFLAG MESS
  double bin0 = -0.05;
  double bin1 = -0.05+0.8875;
  double bin2 = -0.05+2*0.8875;
  double bin3 = -0.05+3*0.8875;

  axranges xiLow{{2, {0.15, 4.0}}, {3, {0.15, 15.0}}, {4, {bin1, bin2}}, {5, {bin1, bin2}}};
  axranges xiMed{{2, {4.0, 8.0}},  {3, {0.15, 15.0}}, {4, {bin1, bin2}}, {5, {bin1, bin2}}};
  axranges xiHig{{2, {8.0, 15.0}}, {3, {0.15, 15.0}}, {4, {bin1, bin2}}, {5, {bin1, bin2}}};
  axranges xiOmega{{4, {bin1, bin2}}, {5, {bin2, bin3}}};
  axranges omegaXi{{4, {bin2, bin3}}, {5, {bin1, bin2}}};
  axranges omegaOmega{{4, {bin2, bin3}}, {5, {bin2, bin3}}};

  // OS
  TH1D *hdphiLowOS = project(hSparseOS, 0, xiLow);
  hdphiLowOS->SetName("hdphiLowOS");
  hdphiLowOS->Draw();
  TH1D *hdphiMedOS = project(hSparseOS, 0, xiMed);
  hdphiMedOS->SetName("hdphiMedOS");
  hdphiMedOS->Draw();
  TH1D *hdphiHigOS = project(hSparseOS, 0, xiHig);
  hdphiHigOS->SetName("hdphiHigOS");
  hdphiHigOS->Draw();
  TH1D *hdphiXiOmOS = project(hSparseOS, 0, xiOmega);
  hdphiXiOmOS->SetName("hdphiXiOmOS");
  hdphiXiOmOS->Draw();
  TH1D *hdphiOmXiOS = project(hSparseOS, 0, omegaXi);
  hdphiOmXiOS->SetName("hdphiOmXiOS");
  hdphiOmXiOS->Draw();
  TH1D *hdphiOmOmOS = project(hSparseOS, 0, omegaOmega);
  hdphiOmOmOS->SetName("hdphiOmOmOS");
  hdphiOmOmOS->Draw();

  // SS
  TH1D *hdphiLowSS = project(hSparseSS, 0, xiLow);
  hdphiLowSS->SetName("hdphiLowSS");
  hdphiLowSS->Draw();
  TH1D *hdphiMedSS = project(hSparseSS, 0, xiMed);
  hdphiMedSS->SetName("hdphiMedSS");
  hdphiMedSS->Draw();
  TH1D *hdphiHigSS = project(hSparseSS, 0, xiHig);
  hdphiHigSS->SetName("hdphiHigSS");
  hdphiHigSS->Draw();
  TH1D *hdphiXiOmSS = project(hSparseSS, 0, xiOmega);
  hdphiXiOmSS->SetName("hdphiXiOmSS");
  hdphiXiOmSS->Draw();
  TH1D *hdphiOmXiSS = project(hSparseSS, 0, omegaXi);
  hdphiOmXiSS->SetName("hdphiOmXiSS");
  hdphiOmXiSS->Draw();
  TH1D *hdphiOmOmSS = project(hSparseSS, 0, omegaOmega);
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