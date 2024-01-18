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
void doQAprojections(TFile* infile) {
  std::vector<TString> plotnames = {"hV0Radius", "hCascRadius", "hV0CosPA", "hCascCosPA", "hDCAPosToPV", "hDCANegToPV", "hDCABachToPV", "hDCAV0ToPV", "hDCAV0Dau", "hDCACascDau", "hLambdaMass", 
    "hITSnClustersPos", "hITSnClustersNeg", "hITSnClustersBach", "hTPCnCrossedRowsPos", "hTPCnCrossedRowsNeg", "hTPCnCrossedRowsBach"};
  TDirectory *dir;
  infile->GetObject("cascade-selector", dir);
  for(TString name : plotnames){
    TH3F *h; 
    dir->GetObject(name, h);
    TH1D* projection = h->ProjectionX();
  }
}

int postprocessingResults(TString trainnr, TString filename = "AnalysisResults.root") {
  TH1::SetDefaultSumw2(); // Make sure we propagate the errors

  TFile *infile;
  if(trainnr == "test"){ // if not a trainnr but just a test, look for AnalysisResults.root in current dir
    infile = new TFile(filename, "READ");
  } else {
    infile = new TFile("results/" + trainnr + "/" + filename, "READ");
  }
  TDirectory *dir;
  infile->GetObject("cascade-correlations", dir);
  // ROOT is being a bitch about TBrowser in batch mode, so make an outfile:
  TFile *outfile = new TFile("plots/" + trainnr + ".root", "RECREATE");
  TDirectory* XiXidir = outfile->mkdir("XiXi");
  TDirectory* OmXidir = outfile->mkdir("OmXi");
  TDirectory* OmOmdir = outfile->mkdir("OmOm");

  // Let's define some global variables:
  double pTmin = 0.15;
  double pTmax = 15.0;

  // ad hoc pT bins:
  // const int maxPtBins = 4;
  // double pTbins[maxPtBins] = {pTmin, 4., 8., pTmax}; // edges of the different (trigger) pT ranges
  // const TString pTlabels[3] = {"Low", "Med", "High"};
  const int maxPtBins = 10;
  double pTbins[maxPtBins] = {pTmin, 1.0, 1.9, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, pTmax}; // edges of the different (trigger) pT ranges
  const TString pTlabels[maxPtBins] = {"0.15", "1", "1.9", "2", "3", "4", "5", "6", "8", "15"};

  // QA plots
  TH1F *hPhi = dir->Get<TH1F>("hPhi");
  hPhi->SetDirectory(outfile);
  hPhi->Draw();

  // doQAprojections(infile);

  // TODO: setrangeuser for mass plots...
  
  // Xi mass
  std::vector<double> siglowXi, sighighXi;
  std::vector<double> bkglowXi, bkghighXi;
  TF1 *f1 = new TF1("f1", "pol2(0) + gaus(3) + gaus(6)", 1.29, 1.42); 
  f1->SetParameters(0, 0, 10, 0, 1.321, 0.005, 0, 1.321, 0.001); 
  f1->SetParLimits(0, -1000, 0);
  f1->SetParLimits(4, 1.321, 1.322);
  f1->SetParLimits(7, 1.321, 1.322);
  f1->SetParLimits(5, 0.001, 0.01);
  f1->SetParLimits(8, 0, 0.01);
  // do inv mass fit in pT bins
  TH2F *hMassXiMinus = dir->Get<TH2F>("hMassXiMinus");
  hMassXiMinus->SetDirectory(outfile);
  cout << "Start Xi mass fitting..." << endl;
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hMassXiMinus->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin+1]);
    TH1D *h = hMassXiMinus->ProjectionX("hMassXiMinus_"+ pTlabels[pTbin + 1]);
    h->GetXaxis()->SetRangeUser(1.2, 1.5);
    h->SetTitle("#Xi^{-} inv. mass for " + pTlabels[pTbin] + " GeV < p_{T} < " + pTlabels[pTbin+1] + " GeV");
    Double_t nentries = h->Integral();
    h->Scale(1. / nentries, "width");
    f1->SetParameter(3, .4* h->GetMaximum());
    f1->SetParameter(6, .4* h->GetMaximum());
    TFitResultPtr r1 = h->Fit("f1", "SLQBR", "", 1.29, 1.42);
    // TODO: automated checks on goodness of fit & wider gaussian > narrow gaussian
    if (f1->GetParameter(5) < f1->GetParameter(8)) cout << "Warning: Gaussian widths are flipped in pT bin " << pTbin << endl;
    std::cout << r1->Chi2() << endl;
    // std::cout << "Check " << f1->GetParameter(5) << " > " << f1->GetParameter(8) << std::endl;
    // std::cout << f1->GetParameter(4) << " and " << f1->GetParameter(7) << std::endl;
    // cout << f1->GetParameter(0) << ", " << f1->GetParameter(1) << ", " << f1->GetParameter(2) << endl; // pol pars
    siglowXi.push_back(f1->GetParameter(4) - 3*f1->GetParameter(5));
    sighighXi.push_back(f1->GetParameter(4) + 3*f1->GetParameter(5));
    bkglowXi.push_back(f1->GetParameter(4) + 6*f1->GetParameter(5));
    bkghighXi.push_back(f1->GetParameter(4) + 12*f1->GetParameter(5));
    cout << "bkg upper edge Xi " << f1->GetParameter(4) + 12*f1->GetParameter(5) << endl;
  }
  // for(auto x : siglowXi) cout << x << endl;

  // Omega mass
  std::vector<double> siglowOm, sighighOm;
  std::vector<double> bkglowOm, bkghighOm;
  TF1 *f2 = new TF1("f2", "pol3(0) + gaus(4) + gaus(7)", 1.64, 1.76); 
  f2->SetParameters(0, 0, 0, 0, 0, 1.672, 0.01, 0, 1.672, 0.005); 
  f2->SetParLimits(5, 1.672, 1.673);
  f2->SetParLimits(8, 1.672, 1.673);
  f2->SetParLimits(6, 0.004, 0.02);
  f2->SetParLimits(9, 0, 0.01);
  // do inv mass fit in pT bins
  TH2F *hMassOmegaMinus = dir->Get<TH2F>("hMassOmegaMinus");
  hMassOmegaMinus->SetDirectory(outfile);
  cout << "Start Omega mass fitting..." << endl;
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hMassOmegaMinus->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin+1]);
    TH1D *h = hMassOmegaMinus->ProjectionX("hMassOmegaMinus_"+ pTlabels[pTbin + 1]);
    h->GetXaxis()->SetRangeUser(1.6, 1.9);
    h->SetTitle("#Omega^{-} inv. mass for " + pTlabels[pTbin] + " GeV < p_{T} < " + pTlabels[pTbin+1] + " GeV");
    Double_t nentries = h->Integral();
    h->Scale(1. / nentries, "width");
    f2->SetParameter(4, 0.3* h->GetMaximum());
    f2->SetParameter(7, 0.4* h->GetMaximum());
    TFitResultPtr r2 = h->Fit("f2", "LQBR", "", 1.64, 1.76);
    // std::cout << "Check " << f2->GetParameter(6) << " > " << f2->GetParameter(9) << std::endl;
    if (f2->GetParameter(6) < f2->GetParameter(9)) cout << "Warning: Gaussian widths are flipped in pT bin " << pTbin << endl;
    siglowOm.push_back(f2->GetParameter(5) - 3*f2->GetParameter(6));
    sighighOm.push_back(f2->GetParameter(5) + 3*f2->GetParameter(6));
    bkglowOm.push_back(f2->GetParameter(5) + 6*f2->GetParameter(6));
    double upperlim = 1.76; // make sure we're not too far out
    if(f2->GetParameter(5) + 12*f2->GetParameter(6) < 1.76) {
      upperlim = f2->GetParameter(5) + 12*f2->GetParameter(6);
    } else {
      cout << "we've cut away (sigma)" << (f2->GetParameter(5) + 12*f2->GetParameter(6) - 1.76)/f2->GetParameter(6) << endl;
    }
    bkghighOm.push_back(upperlim);
    cout << f2->GetParameter(5) + 6*f2->GetParameter(6) << " low bkg high Om " << upperlim << endl;
  }
  // TH1D* h1DXiMass = hMassXiMinus->ProjectionX();
  // TH1D* h1DOmegaMass = hMassOmegaMinus->ProjectionX();

  // Load the THnSparses, 8 in total (4 particle combinations * 2 sign combinations)
  THnSparse *hXiXiOS, *hXiXiSS;
  dir->GetObject("hXiXiOS", hXiXiOS);
  dir->GetObject("hXiXiSS", hXiXiSS);
  THnSparse *hOmXiOS, *hOmXiSS;
  dir->GetObject("hOmXiOS", hOmXiOS);
  dir->GetObject("hOmXiSS", hOmXiSS);
  THnSparse *hOmOmOS, *hOmOmSS;
  dir->GetObject("hOmOmOS", hOmOmOS);
  dir->GetObject("hOmOmSS", hOmOmSS);

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

  // Xi invMass for different pT:
  axranges massLowpT{{ptTrigg, {pTmin, 4.0}}, {selflagTrigg, {0.5, 1.5}}};
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

  // TODO: correlations in pT loop (bins of 1 GeV, can be added/rebinned later). Need tight bins for invMassSelection
  for (int i = 0; i < maxPtBins - 1; i++){
    axranges aXiSig{{ptTrigg, {pTbins[i], pTbins[i+1]}}, {ptAssoc, {pTmin, pTmax}}, 
                    {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}, 
                    {invMassTrigg, {siglowXi[i], sighighXi[i]}}, {invMassAssoc, {siglowXi[i], sighighXi[i]}}
    };
    axranges aXiBkg{{ptTrigg, {pTbins[i], pTbins[i+1]}}, {ptAssoc, {pTmin, pTmax}}, 
                    {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}, 
                    {invMassTrigg, {bkglowXi[i], bkghighXi[i]}}, {invMassAssoc, {bkglowXi[i], bkghighXi[i]}}
    };
    axranges aOmSig{{ptTrigg, {pTbins[i], pTbins[i+1]}}, {ptAssoc, {pTmin, pTmax}}, 
                    {selflagTrigg, {1.5, 3.5}}, {selflagAssoc, {0.5, 1.5}},
                    {invMassTrigg, {siglowOm[i], sighighOm[i]}}, {invMassAssoc, {siglowXi[i], sighighXi[i]}}
    };
    axranges aOmBkg{{ptTrigg, {pTbins[i], pTbins[i+1]}}, {ptAssoc, {pTmin, pTmax}}, 
                    {selflagTrigg, {1.5, 3.5}}, {selflagAssoc, {0.5, 1.5}},
                    {invMassTrigg, {bkglowOm[i], bkghighOm[i]}}, {invMassAssoc, {bkglowXi[i], bkghighXi[i]}}
    };
    XiXidir->cd();
    // Xi-Xi sig
    TH1D *hXiOS = project(hXiXiOS, dPhi, aXiSig);
    hXiOS->SetName("hXiOS"+pTlabels[i]);
    TH1D *hXiSS = project(hXiXiSS, dPhi, aXiSig);
    hXiSS->SetName("hXiSS"+pTlabels[i]);
    TH1D *hXi = new TH1D(*hXiOS);
    hXi->Add(hXiOS, hXiSS, 1, -1);
    hXi->SetName("hXi"+pTlabels[i]);
    // Xi-Xi bkg
    TH1D *hXiOSbkg = project(hXiXiOS, dPhi, aXiBkg);
    hXiOSbkg->SetName("hXiOSbkg"+pTlabels[i]);
    TH1D *hXiSSbkg = project(hXiXiSS, dPhi, aXiBkg);
    hXiSSbkg->SetName("hXiSSbkg"+pTlabels[i]);
    TH1D *hXibkg = new TH1D(*hXiOSbkg);
    hXibkg->Add(hXiOSbkg, hXiSSbkg, 1, -1);
    hXibkg->SetName("hXibkg"+pTlabels[i]);
    // Om-Xi sig
    OmXidir->cd();
    TH1D *hOmOS = project(hOmXiOS, dPhi, aOmSig);
    hOmOS->SetName("hOmOS"+pTlabels[i]);
    TH1D *hOmSS = project(hOmXiSS, dPhi, aOmSig);
    hOmSS->SetName("hOmSS"+pTlabels[i]);
    TH1D *hOm = new TH1D(*hOmOS);
    hOm->Add(hOmOS, hOmSS, 1, -1);
    hOm->SetName("hOm"+pTlabels[i]);
    // Om-Xi bkg
    TH1D *hOmOSbkg = project(hOmXiOS, dPhi, aOmBkg);
    hOmOSbkg->SetName("hOmOSbkg"+pTlabels[i]);
    TH1D *hOmSSbkg = project(hOmXiSS, dPhi, aOmBkg);
    hOmSSbkg->SetName("hOmSSbkg"+pTlabels[i]);
    TH1D *hOmbkg = new TH1D(*hOmOSbkg);
    hOmbkg->Add(hOmOSbkg, hOmSSbkg, 1, -1);
    hOmbkg->SetName("hOmbkg"+pTlabels[i]);

    // Let's do Om Om
    axranges aOmOmSig{{ptTrigg, {pTbins[i], pTbins[i+1]}}, {ptAssoc, {pTmin, pTmax}}, 
                    {selflagTrigg, {1.5, 3.5}}, {selflagAssoc, {1.5, 3.5}},
                    {invMassTrigg, {siglowOm[i], sighighOm[i]}}, {invMassAssoc, {siglowOm[i], sighighOm[i]}}
    };
    axranges aOmOmBkg{{ptTrigg, {pTbins[i], pTbins[i+1]}}, {ptAssoc, {pTmin, pTmax}}, 
                    {selflagTrigg, {1.5, 3.5}}, {selflagAssoc, {1.5, 3.5}},
                    {invMassTrigg, {bkglowOm[i], bkghighOm[i]}}, {invMassAssoc, {bkglowOm[i], bkghighOm[i]}}
    };
    // Om-Om sig
    OmOmdir->cd();
    TH1D *hOmOmOSsig = project(hOmOmOS, dPhi, aOmOmSig);
    hOmOmOSsig->SetName("hOmOmOSsig"+pTlabels[i]);
    TH1D *hOmOmSSsig = project(hOmOmSS, dPhi, aOmOmSig);
    hOmOmSSsig->SetName("hOmOmSSsig"+pTlabels[i]);
    TH1D *hOmOmsig = new TH1D(*hOmOmOSsig);
    hOmOmsig->Add(hOmOmOSsig, hOmOmSSsig, 1, -1);
    hOmOmsig->SetName("hOmOmsig"+pTlabels[i]);
    // Om-Om bkg
    TH1D *hOmOmOSbkg = project(hOmOmOS, dPhi, aOmOmBkg);
    hOmOmOSbkg->SetName("hOmOmOSbkg"+pTlabels[i]);
    TH1D *hOmOmSSbkg = project(hOmOmSS, dPhi, aOmOmBkg);
    hOmOmSSbkg->SetName("hOmOmSSbkg"+pTlabels[i]);
    TH1D *hOmOmbkg = new TH1D(*hOmOmOSbkg);
    hOmOmbkg->Add(hOmOmOSbkg, hOmOmSSbkg, 1, -1);
    hOmOmbkg->SetName("hOmOmbkg"+pTlabels[i]);
    // sig-sig
      // do Xi, Om, both (OS, SS, subtracted)
      // then do logic with Xi + both, Om + both
    // repeat exercise for bkg-bkg
    // remember to save histos in relevant dirs
  }
  outfile->cd();

  // projection configurations
  // axranges xiLow{{ptTrigg, {pTmin, 4.0}}, {ptAssoc, {pTmin, 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  // axranges xiMed{{ptTrigg, {4.0, 8.0}},  {ptAssoc, {0.15, 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  // axranges xiHig{{ptTrigg, {8.0, 15.0}}, {ptAssoc, {0.15, 15.0}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiLow{{ptTrigg, {pTmin, 4.0}}, {ptAssoc, {pTmin, 15.0}}};
  axranges xiMed{{ptTrigg, {4.0, 8.0}},  {ptAssoc, {0.15, 15.0}}};
  axranges xiHig{{ptTrigg, {8.0, 15.0}}, {ptAssoc, {0.15, 15.0}}};

  // OS
  TH1D *hdphiLowOS = project(hXiXiOS, dPhi, xiLow);
  hdphiLowOS->SetName("hdphiLowOS");
  TH1D *hdphiMedOS = project(hXiXiOS, dPhi, xiMed);
  hdphiMedOS->SetName("hdphiMedOS");
  TH1D *hdphiHigOS = project(hXiXiOS, dPhi, xiHig);
  hdphiHigOS->SetName("hdphiHigOS");

  // SS
  TH1D *hdphiLowSS = project(hXiXiSS, dPhi, xiLow);
  hdphiLowSS->SetName("hdphiLowSS");
  TH1D *hdphiMedSS = project(hXiXiSS, dPhi, xiMed);
  hdphiMedSS->SetName("hdphiMedSS");
  TH1D *hdphiHigSS = project(hXiXiSS, dPhi, xiHig);
  hdphiHigSS->SetName("hdphiHigSS");

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

  // Now do the same with inv mass selection for signal: // TODO: FIXME!!!
  axranges xiLowInvMass{{ptTrigg, {0., 4.0}}, {ptAssoc, {0., 15.0}}, {invMassTrigg, {siglowXi[0], sighighXi[0]}}, {invMassAssoc, {siglowXi[0], sighighXi[0]}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiMedInvMass{{ptTrigg, {4.0, 8.0}},  {ptAssoc, {0.15, 15.0}}, {invMassTrigg, {siglowXi[1], sighighXi[1]}}, {invMassAssoc, {siglowXi[1], sighighXi[1]}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};
  axranges xiHigInvMass{{ptTrigg, {8.0, 15.0}}, {ptAssoc, {0.15, 15.0}}, {invMassTrigg, {siglowXi[2], sighighXi[2]}}, {invMassAssoc, {siglowXi[2], sighighXi[2]}}, {selflagTrigg, {0.5, 1.5}}, {selflagAssoc, {0.5, 1.5}}};

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

  // double test = hInvMassXiLow->Integral(hInvMassXiLow->FindBin(siglowXi[0]), hInvMassXiLow->FindBin(sighighXi[0]));
  // cout << test << endl;

  // OS - SS
  TH1D *hdphiLowInv = new TH1D(*hdphiLowOS);
  hdphiLowInv->Add(hdphiLowInvOS, hdphiLowInvSS, 1, -1);
  hdphiLowInv->SetName("hdphiLowInv");
  // hdphiLowInv->Scale()
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