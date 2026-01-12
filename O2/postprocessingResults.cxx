// This script makes the relevant plots/projections from a given AnalysisResults.root.
// Run me like `root 'postprocessingResults.cxx("trainnr", "outputFilename")' -q`, where the second 
//  variable is optional in case the filename diverges from "AnalysisResults.root".

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

// enum for axisnumbers/names:
struct corr{
  enum{
    dPhi, 
    dY, 
    signTrigg,
    signAssoc,
    ptTrigg, 
    ptAssoc, 
    invMassTrigg, 
    invMassAssoc, 
    V_z, 
    multiplicity
  }; 
};
// another enum for the effeciency corrected invariant mass plots:
struct mass{
  enum{
    invMass, 
    sign,
    pT, 
    y, 
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

// // pT bins corresponding to spectra analysis
// const std::vector<double> pTbins = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};
// const std::vector<TString> pTlabels = {"0.6", "1.0", "1.2", "1.4", "1.6", "1.8", "2.0", "2.2", "2.5", "2.9", "3.4", "4.0", "5.0", "6.5"};

// one inclusive bin
const std::vector<double> pTbins = {1.0, 8.0};
const std::vector<TString> pTlabels = {"1.0", "8.0"};

// broader bins:
// const std::vector<double> pTbins = {1.0, 2.5, 4.0, 8.0, 10.0};
// const std::vector<TString> pTlabels = {"1.0", "2.5", "4.0", "8.0", "10.0"};

const int maxPtBins = pTbins.size();
const double pTmin = pTbins[0];
const double pTmax = pTbins[maxPtBins - 1];

// global analysis variables:
const double dYmax = 1.0; // delta-rapidity window to consider [-dYmax, dYmax]
const double signalWindow = 3.0; // in units of sigma
const double sidebandWindow[2] = {4.0, 10.0}; // in units of sigma

// vectors/arrays with boundaries of signal/bkg regions of inv mass plots
// TODO: Sign dependent!
double sigXi[100][2], SBlowXi[100][2], SBhighXi[100][2]; // can't be variable length (i.e. maxPtBins - 1) so just put it to 100.
double sigOm[100][2], bkgOm[100][2]; // can't be variable length (i.e. maxPtBins - 1) so just put it to 100.
double purityXi[100];
// similar for ME histograms, but no sign dependence here
TH2D* MEXiXiSig[100]; TH2D* MEXiXiSB[100];

// initialize some variables, to be set in main
TFile *inputFile;
TFile *outputFile;
TDirectory *inputDir;

// define function that projects all the QA histograms?
void doQAprojections(TFile* infile, bool makePDF = false) { // todo fix infile
  std::vector<TString> plotnames = {"hV0Radius", "hCascRadius", "hV0CosPA", "hCascCosPA", "hDCAPosToPV", "hDCANegToPV", "hDCABachToPV", "hDCAV0ToPV", "hDCAV0Dau", "hDCACascDau", "hLambdaMass", 
    "hITSnClustersPos", "hITSnClustersNeg", "hITSnClustersBach", "hTPCnCrossedRowsPos", "hTPCnCrossedRowsNeg", "hTPCnCrossedRowsBach"};
  TDirectory *dir;
  infile->GetObject("cascade-selector", dir);
  TCanvas *c = new TCanvas();
  for(TString name : plotnames){
    // cout << "   doing " << name << endl; // DEBUG
    TH3F *h; 
    dir->GetObject(name, h);
    TH1D* projection = h->ProjectionX();
    projection->Draw();
    if(makePDF) c->Print("figures/QA/" + name + ".pdf"); 
    c->Clear();
  }
}

// Function for bkg fitting
double pol2bkg(double *x, double *par){
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

void doXiInvMassFits(TString charge, bool makePDF = false){
  assert(charge == "Minus" || charge == "Plus" && "Error in Xi inv. mass fits: charge has to be either 'Minus' or 'Plus'");
  TString sign;
  charge == "Minus" ? sign = "-" : sign = "+";

  TDirectory* XiInvMass = outputFile->mkdir("Xi"+charge+"InvMass");
  XiInvMass->cd();

  TF1 *f1 = new TF1("f1", "pol2(0) + gaus(3) + gaus(6)", 1.29, 1.36); 
  // parameters are: pol{0,1,2}, A1, mu1, sigma1, A2, mu2, sigma2
  f1->SetParameters(0, 0, 10, 0, 1.321, 0.005, 0, 1.321, 0.001); 
  f1->SetParLimits(4, 1.31, 1.33);
  f1->SetParLimits(7, 1.31, 1.33);
  f1->SetParLimits(5, 0, 0.01);
  f1->SetParLimits(8, 0, 0.01);

  TF1 *fBKGXi = new TF1("fBKGXi", pol2bkgXi, 1.29, 1.36, 3);
  TH3F *hXiMass = inputDir->Get<TH3F>("hMassXi" + charge);
  hXiMass->SetDirectory(outputFile);
  cout << "Start Xi mass fitting..." << endl;
  TCanvas *c = new TCanvas();
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hXiMass->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin + 1]);
    TH1D *h = (TH1D*) hXiMass->Project3D("xe");
    h->SetName("hMassXi" + charge + "_" + pTlabels[pTbin + 1]);
    h->SetTitle(";#Lambda#pi inv. mass (GeV/#it{c^{2}}); counts");
    f1->SetParameter(3, .4* h->GetMaximum());
    f1->SetParameter(6, .4* h->GetMaximum());
    h->Fit(fBKGXi, "SLBQRO");
    f1->SetParameter(0, fBKGXi->GetParameter(0));
    f1->SetParameter(1, fBKGXi->GetParameter(1));
    f1->SetParameter(2, fBKGXi->GetParameter(2));

    TFitResultPtr r1 = h->Fit("f1", "SLBQR", "", 1.29, 1.36);
    if (r1->Chi2() > 1) cout << "Warning: Fit Chi2(" << r1->Chi2() << ") > 1 in pT bin " << pTbin << endl;

    double mu = (f1->GetParameter(4) + f1->GetParameter(7)) / 2.;
    double sigma = (f1->GetParameter(5) + f1->GetParameter(8)) / 2.;

    h->GetXaxis()->SetRangeUser(1.28, 1.38);
    h->SetStats(kFALSE);

    // save the mass windows for later analysis
    sigXi[pTbin][0] = mu - signalWindow*sigma; sigXi[pTbin][1] = mu + signalWindow*sigma;
    SBhighXi[pTbin][0] = mu + sidebandWindow[0]*sigma; SBhighXi[pTbin][1] = mu + sidebandWindow[1]*sigma;
    SBlowXi[pTbin][0] = mu - sidebandWindow[1]*sigma; SBlowXi[pTbin][1] = mu - sidebandWindow[0]*sigma;
    cout << sigXi[pTbin][0] << " " << sigXi[pTbin][1] << endl;

    // compute purity
    TF1 *fPol2Xi = new TF1("fPol2Xi", "pol2", 1.29, 1.36);
    fPol2Xi->SetParameters(fBKGXi->GetParameter(0), fBKGXi->GetParameter(1), fBKGXi->GetParameter(2));
    double total = h->Integral(h->GetXaxis()->FindBin(mu - signalWindow*sigma), h->GetXaxis()->FindBin(mu + signalWindow*sigma));
    double bkg = fPol2Xi->Integral(mu - signalWindow*sigma, mu + signalWindow*sigma) / h->GetXaxis()->GetBinWidth(1);
    double purity = (total - bkg) / total;
    purityXi[pTbin] = purity;

    // make nice inv mass plots with fits & purity
    if(makePDF) {
      h->Draw();
      // vertical lines for visualization of signal and bkg region  
      TLine *sig3low = new TLine(mu - signalWindow*sigma, 0, mu - signalWindow*sigma, h->GetMaximum());
      TLine *sig3high = new TLine(mu + signalWindow*sigma, 0, mu + signalWindow*sigma, h->GetMaximum());
      TLine *sig4low = new TLine(mu - sidebandWindow[0]*sigma, 0, mu - sidebandWindow[0]*sigma, h->GetMaximum());
      TLine *sig4high = new TLine(mu + sidebandWindow[0]*sigma, 0, mu + sidebandWindow[0]*sigma, h->GetMaximum());
      TLine *sig10low = new TLine(mu - sidebandWindow[1]*sigma, 0, mu - sidebandWindow[1]*sigma, h->GetMaximum());
      TLine *sig10high = new TLine(mu + sidebandWindow[1]*sigma, 0, mu + sidebandWindow[1]*sigma, h->GetMaximum());
      for (auto line : {sig3low, sig3high, sig4low, sig4high, sig10low, sig10high}){
        line->SetLineColor(kBlack);
        line->SetLineStyle(kDashed);
        line->Draw("same");
      }
      // plot bkg function, also under signal peak
      fPol2Xi->SetLineColor(kGreen+2);
      fPol2Xi->Draw("same");

      TLegend *leg = new TLegend(0.55, 0.65, 0.85, 0.85);
      leg->AddEntry(h, "#Xi^{"+sign+"} invariant mass", "l");
      leg->AddEntry(fPol2Xi, "Background fit function", "l");
      leg->AddEntry(f1, "Total fit function", "l");
      leg->Draw();

      TLatex *purityText = new TLatex();
      purityText->SetNDC();
      purityText->SetTextSize(0.03);
      purityText->DrawLatex(0.55, 0.60, Form("Purity: %.2f%%", purity * 100));
      c->Print("figures/mass/Xi" + charge + pTlabels[pTbin] + "_" + pTlabels[pTbin+1] + ".pdf");
      c->Clear();
    }
  }
  outputFile->cd();
}

void doOmInvMassFits(TString charge, bool makePDF = false){
  assert(charge == "Minus" || charge == "Plus" && "Error in Xi inv. mass fits: charge has to be either 'Minus' or 'Plus'");
  TString sign;
  charge == "Minus" ? sign = "-" : sign = "+";

  TDirectory* OmInvMass = outputFile->mkdir("Om"+charge+"InvMass");
  OmInvMass->cd();

  TF1 *f2 = new TF1("f2", "pol2(0) + gaus(3) + gaus(6)", 1.64, 1.74); 
  // parameters are: pol{0,1,2}, A1, mu1, sigma1, A2, mu2, sigma2
  f2->SetParameters(0, 0, 0, 0, 1.672, 0.01, 0, 1.672, 0.005); 
  f2->SetParLimits(4, 1.665, 1.675);
  f2->SetParLimits(7, 1.665, 1.675);
  f2->SetParLimits(5, 0, 0.02);
  f2->SetParLimits(8, 0, 0.02);
  // do inv mass fit in pT bins

  TF1 *fBKG = new TF1("fBKG", pol2bkg, 1.64, 1.74, 3);

  TH3F *hOmegaMass = inputDir->Get<TH3F>("hMassOmega"+charge);
  hOmegaMass->SetDirectory(outputFile);
  cout << "Start Omega mass fitting..." << endl;
  TCanvas *c = new TCanvas();
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hOmegaMass->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin+1]);
    // TH1D *h = hOmegaMass->ProjectionX("hMassOmega"+charge+"_"+ pTlabels[pTbin + 1]);
    TH1D *h = (TH1D*) hOmegaMass->Project3D("x");
    h->SetName("hMassOmega" + charge + "_" + pTlabels[pTbin + 1]);
    // h->GetXaxis()->SetRangeUser(1.6, 1.9);
    // h->SetTitle("#Omega^{"+sign+"} inv. mass for " + pTlabels[pTbin] + " GeV < p_{T} < " + pTlabels[pTbin+1] + " GeV");
    h->SetTitle(";#Lambda K inv. mass (GeV/#it{c^{2}}); counts");
    f2->SetParameter(3, 0.3* h->GetMaximum());
    f2->SetParameter(6, 0.4* h->GetMaximum());
    h->Fit(fBKG, "SLBQRO");
    f2->SetParameter(0, fBKG->GetParameter(0));
    f2->SetParameter(1, fBKG->GetParameter(1));
    f2->SetParameter(2, fBKG->GetParameter(2));

    TFitResultPtr r2 = h->Fit("f2", "SLQBR", "", 1.64, 1.74);
    if (r2->Chi2() > 1) cout << "Warning: Fit Chi2(" << r2->Chi2() << ") > 1 in pT bin " << pTbin << endl;

    double mu = (f2->GetParameter(4) + f2->GetParameter(7)) / 2.;
    double sigma = (f2->GetParameter(5) + f2->GetParameter(8)) / 2.;

    h->GetXaxis()->SetRangeUser(1.63, 1.72);
    h->SetStats(kFALSE);
    // h->SetLineWidth(3);
    h->Draw();
    if(makePDF) c->Print("figures/Omega" + charge + pTlabels[pTbin] + "_" + pTlabels[pTbin+1] + ".pdf");
    c->Clear();

    if(mu + 10*sigma > 1.76) cout << "Warning: Bkg region exceeds limit of 1.76 GeV in pT bin " << pTbin << endl;
    sigOm[pTbin][0] = mu - 3*sigma; sigOm[pTbin][1] = mu + 3*sigma;
    bkgOm[pTbin][0] = mu + 4*sigma; bkgOm[pTbin][1] = mu + 10*sigma;
    // cout << sigOm[pTbin][0] << " " << sigOm[pTbin][1] << " " << bkgOm[pTbin][0] << " " << bkgOm[pTbin][1] << endl;
  }
  outputFile->cd();
}

void doMixedEvents(bool makePDF = false){

  TCanvas *c = new TCanvas();

  THnSparse *hMEXiXi, *hMEXiOm, *hMEOmOm;
  inputDir->GetObject("MixedEvents/hMEXiXi", hMEXiXi);
  inputDir->GetObject("MixedEvents/hMEXiOm", hMEXiOm);
  inputDir->GetObject("MixedEvents/hMEOmOm", hMEOmOm);
  
  // do ME pT dependent
  for(int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    axranges aMEXiXiSig{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}}, 
                     {corr::dY, {-dYmax, dYmax}},
                     {corr::invMassTrigg, {sigXi[pTbin][0], sigXi[pTbin][1]}}, {corr::invMassAssoc, {sigXi[pTbin][0], sigXi[pTbin][1]}} // fixme read invmass from fits (dedicated function?)
    };

    TH2D *hMEXiXi2D = project2D(hMEXiXi, corr::dY, corr::dPhi, aMEXiXiSig);
    hMEXiXi2D->SetName("hMEXiXi2D"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
    hMEXiXi2D->SetTitle("");
    hMEXiXi2D->RebinX(9);
    // normalize by taking the integral -0.1 < dy < 0.1 and -pi/2 < dphi < pi/2
    int nbinsy = hMEXiXi2D->GetNbinsY();
    assert(nbinsy%2 == 0 && "Error: nbinsy of ME histogram is not even, cannot find center bins!");
    int centerbinY1 = nbinsy/2;
    int centerbinY2 = centerbinY1 + 1;
    hMEXiXi2D->Scale(20. / hMEXiXi2D->Integral(1,10,centerbinY1,centerbinY2)); 
    hMEXiXi2D->SetStats(false);
    hMEXiXi2D->Draw("surf1");
    if(makePDF) c->Print("figures/ME/XiXi"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]+".pdf");
    c->Clear();

    TH1D *hMEXiXidPhi = hMEXiXi2D->ProjectionX();
    hMEXiXidPhi->SetName("hMEXiXidPhi"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
    hMEXiXidPhi->SetTitle("");
    hMEXiXidPhi->SetStats(false);
    hMEXiXidPhi->Scale(1./hMEXiXi2D->GetNbinsY()); // projecting on X requires averaging over bins in Y
    hMEXiXidPhi->Draw();
    if(makePDF) c->Print("figures/ME/XiXidPhi"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]+".pdf");
    c->Clear();
    TH1D *hMEXiXidY = hMEXiXi2D->ProjectionY();
    hMEXiXidY->SetName("hMEXiXidY"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
    hMEXiXidY->SetTitle("");
    hMEXiXidY->SetStats(false);
    hMEXiXidY->Scale(1./hMEXiXi2D->GetNbinsX()); // projecting on Y requires averaging over bins in X
    hMEXiXidY->Draw();
    if(makePDF) c->Print("figures/ME/XiXidY"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]+".pdf");
    c->Clear();

    // Sideband ME
    axranges aMEXiXiSBlow{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}}, 
              {corr::dY, {-dYmax, dYmax}},
              {corr::invMassTrigg, {SBlowXi[pTbin][0], SBlowXi[pTbin][1]}}, {corr::invMassAssoc, {sigXi[pTbin][0], sigXi[pTbin][1]}}
    };
    axranges aMEXiXiSBhigh{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}}, 
              {corr::dY, {-dYmax, dYmax}},
              {corr::invMassTrigg, {SBhighXi[pTbin][0], SBhighXi[pTbin][1]}}, {corr::invMassAssoc, {sigXi[pTbin][0], sigXi[pTbin][1]}}
    };
    TH2D *hMEXiXi2DSBlow = project2D(hMEXiXi, corr::dY, corr::dPhi, aMEXiXiSBlow);
    hMEXiXi2DSBlow->SetName("hMEXiXi2DSBlow"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
    hMEXiXi2DSBlow->SetTitle("");
    hMEXiXi2DSBlow->RebinX(9);
    TH2D *hMEXiXi2DSBhigh = project2D(hMEXiXi, corr::dY, corr::dPhi, aMEXiXiSBhigh);
    hMEXiXi2DSBhigh->SetName("hMEXiXi2DSBhigh"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
    hMEXiXi2DSBhigh->SetTitle("");
    hMEXiXi2DSBhigh->RebinX(9);
    // simply add the two sidebands before normalizing
    TH2D *hMEXiXi2DSB = new TH2D(*hMEXiXi2DSBlow);
    hMEXiXi2DSB->SetName("hMEXiXi2DSB"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
    hMEXiXi2DSB->Add(hMEXiXi2DSBhigh, hMEXiXi2DSBlow);
    hMEXiXi2DSB->Scale(20. / hMEXiXi2DSB->Integral(1,10,centerbinY1,centerbinY2)); 
    hMEXiXi2DSB->SetStats(false);
    hMEXiXi2DSB->Draw("surf1");
    if(makePDF) c->Print("figures/ME/XiXiSB"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]+".pdf");
    c->Clear();

    TH1D *hMEXiXi2DSBdPhi = hMEXiXi2DSB->ProjectionX();
    hMEXiXi2DSBdPhi->SetName("hMEXiXi2DSBdPhi"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
    hMEXiXi2DSBdPhi->SetTitle("");
    hMEXiXi2DSBdPhi->SetStats(false);
    hMEXiXi2DSBdPhi->Scale(1./hMEXiXi2DSB->GetNbinsY());

    TH1D *hMEXiXi2DSBdY = hMEXiXi2DSB->ProjectionY();
    hMEXiXi2DSBdY->SetName("hMEXiXi2DSBdY"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
    hMEXiXi2DSBdY->SetTitle("");
    hMEXiXi2DSBdY->SetStats(false);
    hMEXiXi2DSBdY->Scale(1./hMEXiXi2DSB->GetNbinsX());

    MEXiXiSB[pTbin] = hMEXiXi2DSB;
    MEXiXiSig[pTbin] = hMEXiXi2D;

    // // vvv XI-OMEGA ME vvv
    // // vvvvvvvvvvvvvvvvvvv
    // axranges aMEXiOm{{corr::ptTrigg, {1., 5.0}}, {corr::ptAssoc, {1., 5.0}},
    //                    {corr::dY, {-dYmax, dYmax}},
    //                   {corr::invMassTrigg, {1.66, 1.685}}, {corr::invMassAssoc, {1.31, 1.33}}
    // };
    // // Okay let's take Om-Xi OS combination to make ME projection in rapidity and phi (all pT ranges)
    // TH1D *hXiOmdY = project(hXiOm, corr::dY, aMEXiOm);
    // hXiOmdY->SetName("hMEXiOmdY");
    // hXiOmdY->SetTitle("ME #Xi-#Omega (pT integrated)");
    // hXiOmdY->SetYTitle("pairs OS");
    // hXiOmdY->SetLineWidth(3);
    // hXiOmdY->Scale(1. / hXiOmdY->GetMaximum());

    // TH1D *hXiOmdPhi = project(hXiOm, corr::dPhi, aMEXiOm);
    // hXiOmdPhi->SetName("hMEXiOmdPhi");
    // hXiOmdPhi->SetTitle("ME #Xi-#Omega (pT integrated)");
    // hXiOmdPhi->SetYTitle("pairs OS");
    // hXiOmdPhi->SetLineWidth(3);
    // hXiOmdPhi->Rebin(9);
    // hXiOmdPhi->Scale(1. / hXiOmdPhi->GetMaximum());
  }
}

// Helper function to perform OS/SS combinations and subtraction
template<typename T>
void doSubAndAvg(T** array, const TString& basename, 
  int pTbin, const TString& pTlabel1, const TString& pTlabel2) {
  // Takes an array[7] = {-+, +-, --, ++, _, _, _} of histograms
  // computes the OS and SS averages and OS-SS subtraction, saving the results in array[4], array[5], array[6] respectively
  // can deal with both TH1D and TH2D

  T *hOS = new T(*array[0]);
  T *hSS = new T(*array[0]);

  hOS->Add(array[0], array[1], 0.5, 0.5);
  hOS->SetName(basename + "OS_pT_" + pTlabel1 + "_" + pTlabel2);
  hOS->SetTitle("");
  array[4] = hOS;

  hSS->Add(array[2], array[3], 0.5, 0.5);
  hSS->SetName(basename + "SS_pT_" + pTlabel1 + "_" + pTlabel2);
  hSS->SetTitle("");
  array[5] = hSS;

  T *hSub = new T(*hOS);
  hSub->Add(hOS, hSS, 1, -1);
  hSub->SetName(basename + "Sub_pT_" + pTlabel1 + "_" + pTlabel2);
  array[6] = hSub;
}

int postprocessingResults(TString trainnr, TString filename = "AnalysisResults.root", bool makePDF = false) {
  // checks based on the global scope variables:
  assert((pTbins.size() == pTlabels.size()) && "pTbins and pTlabels have different sizes, something is wrong!");

  TH1::SetDefaultSumw2(); // Make sure we propagate the errors
  gStyle->SetHistLineWidth(3);
  gStyle->SetAxisMaxDigits(2); // max digits for axis labels, prevents labels like 0.00001
  gStyle->SetPadLeftMargin(0.1); // increase left margin for y axis titles
  gStyle->SetPadRightMargin(0.05); 
  gROOT->ForceStyle();

  if(trainnr == "test"){ // if not a trainnr but just a test, look for AnalysisResults.root in current dir
    inputFile = new TFile(filename, "READ");
  } else {
    inputFile = new TFile("results/" + trainnr + "/" + filename, "READ");
  }
  inputFile->GetObject("cascade-correlations", inputDir);
  outputFile = new TFile("plots/" + trainnr + ".root", "RECREATE");
  TDirectory* XiXidir = outputFile->mkdir("XiXi");
  TDirectory* OmXidir = outputFile->mkdir("OmXi");
  TDirectory* OmOmdir = outputFile->mkdir("OmOm");
  
  // make sure inputFiles and stuff are not null
  assert((!inputFile && !inputDir && !outputFile) && "input or output file is null pointer!");
  TCanvas *c = new TCanvas("c");
  c->cd();
  gPad->SetTheta(60); // test POV for 2D correlations

  // QA plots
  cout << "-- QA PROJECTIONS --" << endl;
  TH1F *hPhi = inputDir->Get<TH1F>("hPhi");
  hPhi->SetDirectory(outputFile);
  // hPhi->Draw();

  doQAprojections(inputFile);

  // test inv mass fitting
  cout << "-- INV MASS FITS --" << endl;
  doXiInvMassFits("Minus", makePDF);
  doXiInvMassFits("Plus", makePDF);
  doOmInvMassFits("Minus", makePDF);
  doOmInvMassFits("Plus", makePDF);
  // todo: assert some properties of the inv mass region boundaries to ensure everything went ok.

  // test ME
  cout << "-- MIXED EVENTS --" << endl;
  doMixedEvents(makePDF);

  // now that we have the inv mass regions, let's project and integrate the eff corrected inv mass plots to determine the number of triggers
  cout << "-- NORMALIZATION --" << endl;
  TDirectory* XiEffdir = outputFile->mkdir("XiEffdir");
  TDirectory* OmEffdir = outputFile->mkdir("OmEffdir");
  THnSparse *hEffCorrXiMass, *hEffCorrOmegaMass;
  inputDir->GetObject("hMassXiEffCorrected", hEffCorrXiMass);
  inputDir->GetObject("hMassOmegaEffCorrected", hEffCorrOmegaMass);
  // put the projections in a vector so we can access them later
  // array of Ntriggers for Xi, Omega
  double NTrigXi[maxPtBins - 1], NTrigOm[maxPtBins - 1];
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    axranges aXi{{mass::pT, {pTbins[pTbin], pTbins[pTbin + 1]}},
                 {mass::invMass, {sigXi[pTbin][0], sigXi[pTbin][1]}},
    };
    axranges aOm{{mass::pT, {pTbins[pTbin], pTbins[pTbin + 1]}},
                 {mass::invMass, {sigOm[pTbin][0], sigOm[pTbin][1]}},
    };
    TH1D *hXi = project(hEffCorrXiMass, mass::invMass, aXi);
    hXi->SetName("hEffXiMass_" + pTlabels[pTbin + 1]);
    hXi->SetDirectory(XiEffdir);
    TH1D *hXi1 = project(hEffCorrXiMass, mass::pT, aXi);
    hXi1->SetName("hEffXiMasspT_" + pTlabels[pTbin + 1]);
    hXi1->SetDirectory(XiEffdir);
    TH1D *hXi2 = project(hEffCorrXiMass, mass::y, aXi);
    hXi2->SetName("hEffXiMassy_" + pTlabels[pTbin + 1]);
    hXi2->SetDirectory(XiEffdir);
    TH2D *hXi3 = project2D(hEffCorrXiMass, mass::pT, mass::y, aXi);
    hXi3->SetName("hEffXiMasspTy_" + pTlabels[pTbin + 1]);
    hXi3->SetDirectory(XiEffdir);
    
    double xiInt = hXi->Integral();
    if(std::isinf(xiInt)) {
      xiInt = 1.; // set to 1 in case of infinity, warn the user
      cout << "WARNING: Ntrig Xi in pT bin " << pTlabels[pTbin] << " - " << pTlabels[pTbin + 1] << " is inf - Ntrig set to 1." << endl;
    }
    NTrigXi[pTbin] = xiInt;
    // cout << "ntrig = " << ntrig << endl;
    TH1D *hOm = project(hEffCorrOmegaMass, mass::invMass, aOm);
    hOm->SetName("hEffOmMass_" + pTlabels[pTbin + 1]);
    hOm->SetDirectory(OmEffdir);
    double omInt = hOm->Integral();
    if(std::isinf(omInt)) {
      omInt = 1.; // set to 1 in case of infinity, warn the user
      cout << "WARNING: Ntrig Om in pT bin " << pTlabels[pTbin] << " - " << pTlabels[pTbin + 1] << " is inf - Ntrig set to 1." << endl;
    }
    NTrigOm[pTbin] = omInt;
  }
  
  // Load the THnSparses
  THnSparse *hXiXi, *hXiOm, *hOmOm;
  inputDir->GetObject("hXiXi", hXiXi);
  inputDir->GetObject("hXiOm", hXiOm);
  inputDir->GetObject("hOmOm", hOmOm);

  // test projection of Xi-Xi OS onto 2D dphi-dy
  axranges aPtIntMass{{mass::pT, {1., 10.}}, 
                      // {mass::y, {-0.5, 0.5}},
                      {mass::invMass, {sigXi[0][0], sigXi[0][1]}}
  };
  // make single Xi spectra (charge independent)
  // get Nevents:
  // double Nevents = (double) inputFile->Get<TH1I>("cascade-selector/hEventSel")->GetBinContent(6);
  double Nevents = 12999723968.; // REMOVE ME AFTER FIXING MAX INT
  cout << "Nevents: " << Nevents << endl;
  TH1D* hXiSpectra = project(hEffCorrXiMass, mass::pT, {{mass::pT, {0., pTmax}}, {mass::y, {-0.5, 0.5}}, {mass::invMass, {1.31, 1.33}}});
  hXiSpectra->SetName("hXiSpectra");
  Double_t Xibinning[14] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};
  TH1D* hSpectraRebin = (TH1D*) hXiSpectra->Rebin(13, "hSpectraRebin", Xibinning);
  hSpectraRebin->Scale(1./(Nevents), "width"); // divide by an extra factor of 2 because of Xi charges ---- 65536 = factor of 2^n used in run 3 AN
  // hSpectraRebin->Scale(1209904670./2024331520.); // Nev(rec) = 1209904670 / Nev(gen) = 2024331520 (run 521456)
  hSpectraRebin->Scale(1199241785./2019676955.); // Nev(rec) = 1199241785 / Nev(gen) = 2019676955 (run 550154 24f3c_fix_medium)
  TH3F *hMassXiMinus = inputDir->Get<TH3F>("hMassXiMinus");
  TH1D *hXiMinusMass1D = hMassXiMinus->ProjectionX();
  hXiMinusMass1D->SetName("hXiMinusMass1D");
  TH1D *hRawXiMass = hMassXiMinus->ProjectionY("hRawXiMass");
  TH1D* hRawSpectraRebin = (TH1D*) hRawXiMass->Rebin(13, "hRawSpectraRebin", Xibinning);
  hRawSpectraRebin->Scale(1./Nevents, "width"); // ---- 65536 = factor of 2^n used in run 3 AN

  TH1D* hXiSpectraMult = project(hEffCorrXiMass, mass::multiplicity, {{mass::pT, {0.6, pTmax}}, {mass::y, {-0.5, 0.5}}, {mass::invMass, {1.31, 1.33}}, {mass::multiplicity, {0,99}}});
  hXiSpectraMult->SetName("hXiSpectraMult");

  TFile *run2file = new TFile("run2Xispectra.root", "READ");
  outputFile->cd();
  TH1D *run2Spectra = (TH1D*) run2file->Get<TH1F>("Table 3/Hist1D_y11");
  run2Spectra->SetLineColor(kRed);
  hSpectraRebin->GetYaxis()->SetRangeUser(1e-5, 0.045); // don't put minimum to zero, so we can use logscale later
  TCanvas *cSpectra = new TCanvas("cSpectra");
  hSpectraRebin->SetTitle("Xi spectra");
  hSpectraRebin->GetYaxis()->SetTitle("\\frac{1}{N_{ev}}\\frac{d^{2}N}{dy dp_{T}} (GeV/c)^{-1}");
  hSpectraRebin->SetStats(kFALSE);
  hSpectraRebin->Draw();
  run2Spectra->Draw("SAME");
  TLegend *legend = new TLegend();
  legend->AddEntry(hSpectraRebin, "this analysis");
  legend->AddEntry(run2Spectra, "run 2");
  legend->Draw();
  cSpectra->Write();
  TH1D* hSpectraRatio = (TH1D*) hSpectraRebin->Clone();
  hSpectraRatio->SetName("hSpectraRatio");
  hSpectraRatio->GetXaxis()->SetLimits(0.6, 6.5);
  Double_t errors[13] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  run2Spectra->SetError(errors);
  hSpectraRatio->Divide(run2Spectra);

  TCanvas *cRatio = new TCanvas("cRatio");
  TRatioPlot *rp = new TRatioPlot(hSpectraRebin, run2Spectra, "divsym");
  rp->SetH1DrawOpt("E");
  rp->Draw();
  rp->GetUpperPad()->cd();
  rp->GetUpperPad()->SetLogy();
  rp->GetXaxis()->SetLimits(0.6, 6.5);
  rp->GetXaxis()->SetTitle("p_T");
  legend->Draw();
  if(makePDF) cRatio->Print("figures/spectraRatio.pdf");
  cRatio->Write();

  // CORRELATIONS

  // TEMP remove this
  // make quick pT spectrum of eff corrected xi, to make ratio
  TH1D *hXiEffpT = project(hEffCorrXiMass, mass::pT, aPtIntMass);
  hXiEffpT->SetName("TEMPXiEffpT");

  axranges aPtIntegrated{{corr::signTrigg, {-1, -1}}, {corr::signAssoc, {1, 1}},
                        {corr::ptTrigg, {1., 10}}, {corr::ptAssoc, {1., 10}},
                        {corr::dY, {-dYmax, dYmax}},
                        {corr::invMassTrigg, {sigXi[0][0], sigXi[0][1]}}, {corr::invMassAssoc, {sigXi[0][0], sigXi[0][1]}}
  };

  axranges aXiMinXiPlus{{corr::signTrigg, {-1, -1}}, {corr::signAssoc, {1, 1}},
                        {corr::ptTrigg, {1., 10}}, {corr::ptAssoc, {1., 10}},
                        {corr::dY, {-dYmax, dYmax}},
                        {corr::invMassTrigg, {sigXi[0][0], sigXi[0][1]}}, {corr::invMassAssoc, {sigXi[0][0], sigXi[0][1]}}
  };
  axranges aXiMinXiMin{{corr::signTrigg, {-1, -1}}, {corr::signAssoc, {-1, -1}},
                      {corr::ptTrigg, {1., 10}}, {corr::ptAssoc, {1., 10}},
                      {corr::dY, {-dYmax, dYmax}},
                      {corr::invMassTrigg, {sigXi[0][0], sigXi[0][1]}}, {corr::invMassAssoc, {sigXi[0][0], sigXi[0][1]}}
  };
  axranges aXiPlusXiMin{{corr::signTrigg, {1, 1}}, {corr::signAssoc, {-1, -1}},
                        {corr::ptTrigg, {1., 10}}, {corr::ptAssoc, {1., 10}},
                        {corr::dY, {-dYmax, dYmax}},
                        {corr::invMassTrigg, {sigXi[0][0], sigXi[0][1]}}, {corr::invMassAssoc, {sigXi[0][0], sigXi[0][1]}}
  };
  axranges aXiPlusXiPlus{{corr::signTrigg, {1, 1}}, {corr::signAssoc, {1, 1}},
                        {corr::ptTrigg, {1., 10}}, {corr::ptAssoc, {1., 10}},
                        {corr::dY, {-dYmax, dYmax}},
                        {corr::invMassTrigg, {sigXi[0][0], sigXi[0][1]}}, {corr::invMassAssoc, {sigXi[0][0], sigXi[0][1]}}
  };

  // inv mass of Xi trig & assoc
  TH1D *hXiInvMassTrig = project(hXiXi, corr::invMassTrigg, aPtIntegrated);
  hXiInvMassTrig->SetName("hXiInvMassTrig");
  TH1D *hXiInvMassAssoc = project(hXiXi, corr::invMassAssoc, aPtIntegrated);
  hXiInvMassAssoc->SetName("hXiInvMassAssoc");

  // TRIGGER NORMALIZATION
  axranges aMassXiMin{{mass::pT, {1., 10.}}, 
                      {mass::sign, {-1, -1}},
                      {mass::invMass, {sigXi[0][0], sigXi[0][1]}}
  };
  axranges aMassXiPlus{{mass::pT, {1., 10.}}, 
                      {mass::sign, {1, 1}},
                      {mass::invMass, {sigXi[0][0], sigXi[0][1]}}
  };
  
  TH1D *hMassTest = hEffCorrXiMass->Projection(0);
  hMassTest->SetName("hMassTest");
  TH1D *hXiMinMassPtInt = project(hEffCorrXiMass, mass::invMass, aMassXiMin);
  hXiMinMassPtInt->SetName("hEffXiMinMass_pTint");
  double nXiMin = hXiMinMassPtInt->Integral();
  TH1D *hXiPlusMassPtInt = project(hEffCorrXiMass, mass::invMass, aMassXiPlus);
  hXiPlusMassPtInt->SetName("hEffXiPlusMass_pTint");
  double nXiPlus = hXiPlusMassPtInt->Integral();
  cout << "NXiMin: " << nXiMin << ", NXiPlus: " << nXiPlus << endl;

  // test ME with perfect triangle (dY) and uniform (dPhi) correction:
  TH2D *hTestME = new TH2D("hTestME", "hTestME", 20, -0.5*M_PI, 1.5*M_PI, 20, -1, 1);
  for (int dphibin = 0; dphibin < 20; dphibin++){
    double x = (dphibin+0.5)*2*M_PI/20. - 0.5*M_PI;
    for (int dybin = 0; dybin < 20; dybin++){
      double y = (dybin+0.5)*0.1 - 1;
      hTestME->Fill(x, y, 1-std::abs(y));
      hTestME->SetBinError(hTestME->FindBin(x,y), 0);
    }
  }
  
  c->cd(); //not sure why, but need this here
  // actual 2D ME, just take Xi-Xi all sign combinations
  // this is for the run 2 comparison
  THnSparse *hMEXiXi;
  inputDir->GetObject("MixedEvents/hMEXiXi", hMEXiXi);
  axranges aMEXiXi{{corr::ptTrigg, {1.0, 10.}}, {corr::ptAssoc, {1.0, 10.}}, 
                  {corr::dY, {-dYmax, dYmax}},
                  {corr::invMassTrigg, {sigXi[0][0], sigXi[0][1]}}, {corr::invMassAssoc, {sigXi[0][0], sigXi[0][1]}} // fixme read invmass from fits (dedicated function?)
  };
  TH2D *hMEXiXi2D = project2D(hMEXiXi, corr::dY, corr::dPhi, aMEXiXi);
  hMEXiXi2D->SetName("hMEXiXi2D");
  hMEXiXi2D->SetTitle("");
  hMEXiXi2D->RebinX(9);
  // normalize by taking the integral -0.1 < dy < 0.1 and -pi/2 < dphi < pi/2
  // hardcoded! integrate over 20 bins, so divide (multiply) by 20 to take average
  int nbinsy = hMEXiXi2D->GetNbinsY();
  cout << "nbinsy ME " << nbinsy << endl;
  assert(nbinsy%2 == 0 && "Error: nbinsy of ME histogram is not even, cannot find center bins!");
  int centerbinY1 = nbinsy/2;
  int centerbinY2 = centerbinY1 + 1;
  hMEXiXi2D->Scale(20. / hMEXiXi2D->Integral(1,10,centerbinY1,centerbinY2)); 
  hMEXiXi2D->SetStats(false);
  hMEXiXi2D->Draw("surf1");
  // c->Write("cMEXiXi2D");
  if(makePDF) c->Print("figures/ME/XiXi.pdf");
  c->Clear();

  TH1D *hMEXiXidPhi = hMEXiXi2D->ProjectionX();
  hMEXiXidPhi->SetName("hMEXiXidPhi");
  hMEXiXidPhi->SetTitle("");
  hMEXiXidPhi->SetStats(false);
  hMEXiXidPhi->Scale(1./hMEXiXi2D->GetNbinsY()); // projecting on X requires averaging over bins in Y
  hMEXiXidPhi->Draw();
  if(makePDF) c->Print("figures/ME/XiXidPhi.pdf");
  c->Clear();
  TH1D *hMEXiXidY = hMEXiXi2D->ProjectionY();
  hMEXiXidY->SetName("hMEXiXidY");
  hMEXiXidY->SetTitle("");
  hMEXiXidY->SetStats(false);
  hMEXiXidY->Scale(1./hMEXiXi2D->GetNbinsX()); // projecting on Y requires averaging over bins in X
  hMEXiXidY->Draw();
  if(makePDF) c->Print("figures/ME/XiXidY.pdf");
  c->Clear();

  //2D

  TH2D* hXiMinXiPlus = project2D(hXiXi, corr::dY, corr::dPhi, aXiMinXiPlus);
  hXiMinXiPlus->SetName("hXiMinXiPlus");
  hXiMinXiPlus->RebinX(9);
  hXiMinXiPlus->Scale(1. / nXiMin);
  hXiMinXiPlus->Divide(hMEXiXi2D);
  hXiMinXiPlus->SetTitle("#Xi^{-} - #Xi^{+} (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hXiMinXiPlus->SetStats(false);
  TH1D* hXiMinXiPlus1D = hXiMinXiPlus->ProjectionX("hXiMinXiPlus1D");
  hXiMinXiPlus1D->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  TH2D* hXiMinXiMin = project2D(hXiXi, corr::dY, corr::dPhi, aXiMinXiMin);
  hXiMinXiMin->SetName("hXiMinXiMin");
  hXiMinXiMin->RebinX(9);
  hXiMinXiMin->Scale(1. / nXiMin);
  hXiMinXiMin->Divide(hMEXiXi2D);
  hXiMinXiMin->SetTitle("#Xi^{-} - #Xi^{-} (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hXiMinXiMin->SetStats(false);
  TH1D* hXiMinXiMin1D = hXiMinXiMin->ProjectionX("hXiMinXiMin1D");
  hXiMinXiMin1D->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  TH2D* hXiPlusXiMin = project2D(hXiXi, corr::dY, corr::dPhi, aXiPlusXiMin);
  hXiPlusXiMin->SetName("hXiPlusXiMin");
  hXiPlusXiMin->RebinX(9);
  hXiPlusXiMin->Scale(1. / nXiPlus);
  hXiPlusXiMin->Divide(hMEXiXi2D);
  hXiPlusXiMin->SetTitle("#Xi^{+} - #Xi^{-} (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hXiPlusXiMin->SetStats(false);
  TH1D* hXiPlusXiMin1D = hXiPlusXiMin->ProjectionX("hXiPlusXiMin1D");
  hXiPlusXiMin1D->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  TH2D* hXiPlusXiPlus = project2D(hXiXi, corr::dY, corr::dPhi, aXiPlusXiPlus);
  hXiPlusXiPlus->SetName("hXiPlusXiPlus");
  hXiPlusXiPlus->RebinX(9);
  hXiPlusXiPlus->Scale(1. / nXiPlus);
  hXiPlusXiPlus->SetTitle("#Xi^{+} - #Xi^{+} (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hXiPlusXiPlus->SetStats(false);
  TH1D* hXiPlusXiPlus1D = hXiPlusXiPlus->ProjectionX("hXiPlusXiPlus1D");
  hXiPlusXiPlus1D->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)

  TH2D *hXiXiOS2D = new TH2D(*hXiMinXiPlus);
  hXiXiOS2D->Add(hXiMinXiPlus, hXiPlusXiMin, 0.5, 0.5);
  hXiXiOS2D->SetName("hXiXiOS2D");
  hXiXiOS2D->SetTitle("OS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c})");
  TH2D *hXiXiSS2D = new TH2D(*hXiMinXiMin);
  hXiXiSS2D->Add(hXiMinXiMin, hXiPlusXiPlus, 0.5, 0.5);
  hXiXiSS2D->SetName("hXiXiSS2D");
  hXiXiSS2D->SetTitle("SS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c})");
  TH2D *hXiXi2D = new TH2D(*hXiXiOS2D);
  hXiXi2D->Add(hXiXiOS2D, hXiXiSS2D, 1, -1);
  hXiXi2D->SetName("hXiXi2D");
  hXiXi2D->SetTitle("OS - SS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c})");

  TH1D *hXiXiOS = hXiXiOS2D->ProjectionX("hXiXiOS");
  hXiXiOS->SetTitle("OS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hXiXiOS->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  TH1D *hXiXiSS = hXiXiSS2D->ProjectionX("hXiXiSS");
  hXiXiSS->SetTitle("SS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hXiXiSS->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  TH1D *hXiXisubtracted = hXiXi2D->ProjectionX("hXiXisubtracted");
  hXiXisubtracted->SetTitle("OS - SS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hXiXisubtracted->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)

  ///// RUN 2 COMPARISONS

  const double datapointsOS[21] = {0, 0.00772392, 0.00763424, 0.0102429, 0.00970434, 0.0147818, 0.0131757, 0.0114212, 0.00877517, 0.00667961, 0.00695319, 0.00729105, 0.00602286, 0.00449109, 0.00573827, 0.00629415, 0.00505972, 0.00429942, 0.00724581, 0.00570471, 0.00750361};
  const double staterrorsOS[21] = {0, 0.00140528, 0.00105404, 0.00141499, 0.000748423, 0.00128475, 0.00101308, 0.00136926, 0.00085206, 0.000793002, 0.000902522, 0.000955686, 0.000843041, 0.000948397, 0.000851573, 0.000820212, 0.00101491, 0.00104776, 0.000878481, 0.000767664, 0.00100554};
  TH1D* hRun2OS = new TH1D(*hXiXiOS);
  hRun2OS->SetName("hRun2OS");
  hRun2OS->SetTitle("Run 2 OS #Xi - #Xi");
  hRun2OS->SetContent(datapointsOS);
  hRun2OS->SetError(staterrorsOS);
  hRun2OS->SetStats(false);
  hRun2OS->SetLineColor(kRed);
  TRatioPlot *hRun2OSRatio = new TRatioPlot(hXiXiOS, hRun2OS, "divsym");
  hRun2OSRatio->SetH1DrawOpt("E");
  hRun2OSRatio->Draw("nogrid");
  hRun2OSRatio->GetUpperRefYaxis()->SetRangeUser(0., 0.02);
  hRun2OSRatio->GetLowerRefYaxis()->SetRangeUser(0., 1.);
  hRun2OSRatio->GetLowYaxis()->SetNdivisions(505);
  TLegend *legendOS = new TLegend(0.55, 0.75, 0.75, 0.85);
  hRun2OSRatio->GetUpperPad()->cd(); // draw legend in upper pad
  legendOS->AddEntry(hXiXiOS, "this analysis");
  legendOS->AddEntry(hRun2OS, "run 2");
  legendOS->Draw();
  c->Write("cOSRatio");
  if(makePDF) c->Print("figures/OSRatio.pdf");
  c->Clear();
  
  const double datapointsSS[21] = {0, 0.00314547, 0.00574206, 0.0049273, 0.00321404, 0.00306671, 0.00141205, 0.00276105, 0.00273416, 0.00406806, 0.00376302, 0.00390338, 0.0048735, 0.00529755, 0.00359733, 0.0080713, 0.00588222, 0.00452432, 0.00337606, 0.00428579, 0.00485876};
  const double staterrorsSS[21] = {0, 0.000688171, 0.00119578, 0.000846638, 0.00051905, 0.000517229, 0.00148369, 0.000717417, 0.000644928, 0.000621341, 0.000561713, 0.000813817, 0.00107704, 0.000920041, 0.00073849, 0.00245958, 0.000746747, 0.000804554, 0.000785611, 0.000791499, 0.000860004};
  TH1D* hRun2SS = new TH1D(*hXiXiSS);
  hRun2SS->SetName("hRun2SS");
  hRun2SS->SetTitle("Run 2 SS #Xi - #Xi");
  hRun2SS->SetContent(datapointsSS);
  hRun2SS->SetError(staterrorsSS);
  hRun2SS->SetStats(false);
  hRun2SS->SetLineColor(kRed);
  TRatioPlot *hRun2SSRatio = new TRatioPlot(hXiXiSS, hRun2SS, "divsym");
  hRun2SSRatio->SetH1DrawOpt("E");
  hRun2SSRatio->Draw("nogrid");
  hRun2SSRatio->GetUpperRefYaxis()->SetRangeUser(0., 0.01);
  hRun2SSRatio->GetLowerRefYaxis()->SetRangeUser(0., 1.);
  hRun2SSRatio->GetLowYaxis()->SetNdivisions(505);
  TLegend *legendSS = new TLegend(0.2, 0.75, 0.4, 0.85);
  hRun2SSRatio->GetUpperPad()->cd(); // draw legend in upper pad
  legendSS->AddEntry(hXiXiSS, "this analysis");
  legendSS->AddEntry(hRun2SS, "run 2");
  legendSS->Draw();
  c->Write("cSSRatio");
  if(makePDF) c->Print("figures/SSRatio.pdf");
  c->Clear();
  
  const double datapoints[21] = {0, 0.00457845, 0.00189219, 0.00531563, 0.00649031, 0.0117151, 0.0117636, 0.00866012, 0.00604101, 0.00261155, 0.00319017, 0.00338767, 0.00114936, -0.000806458, 0.00214094, -0.00177715, -0.000822495, -0.000224905, 0.00386975, 0.00141892, 0.00264486};
  const double staterrors[21] = {0, 0.00156473, 0.00159401, 0.00164894, 0.000910796, 0.00138496, 0.00179657, 0.00154582, 0.00106861, 0.00100743, 0.00106305, 0.00125524, 0.00136775, 0.00132134, 0.00112718, 0.00259274, 0.00126003, 0.00132103, 0.00117852, 0.00110262, 0.00132315};
  TH1D* hRun2 = new TH1D(*hXiXisubtracted);
  hRun2->SetName("hRun2");
  hRun2->SetTitle("Run 2 #Xi - #Xi");
  hRun2->SetContent(datapoints);
  hRun2->SetError(staterrors);
  hRun2->SetStats(false);
  hRun2->SetLineColor(kRed);
  TRatioPlot *hRun2Ratio = new TRatioPlot(hXiXisubtracted, hRun2, "divsym");
  // TH1D* hRun2Ratio = new TH1D(*hDPhicorrected);
  hRun2Ratio->SetH1DrawOpt("E");
  hRun2Ratio->Draw("nogrid");
  hRun2Ratio->GetUpperRefYaxis()->SetRangeUser(0., 0.015);
  hRun2Ratio->GetLowerRefYaxis()->SetRangeUser(0., 1.);
  hRun2Ratio->GetLowYaxis()->SetNdivisions(505);
  TLegend *legend2 = new TLegend(0.55, 0.75, 0.75, 0.85);
  hRun2Ratio->GetUpperPad()->cd(); // draw legend in upper pad
  legend2->AddEntry(hXiXisubtracted, "this analysis");
  legend2->AddEntry(hRun2, "run 2");
  legend2->Draw();
  c->Write("cRun2Ratio");
  if(makePDF) c->Print("figures/Run2Ratio.pdf");
  c->Clear();

  ///// END RUN 2 COMPARISONS

  cout << "Xi-Xi yield " << hXiXisubtracted->Integral() << endl;

  // arrays to store sign combinations in: 
  TH2D *arrayXi[maxPtBins - 1][7]; // {-+, +-, --, ++, OS, SS, OS-SS}
  TH1D *arrayXi1D[maxPtBins - 1][7]; // {-+, +-, --, ++, OS, SS, OS-SS}
  TH1D *arrayXiRaw[maxPtBins - 1][7]; // {-+, +-, --, ++, OS, SS, OS-SS}
  TH1D *arrayXi1DSBlow[maxPtBins - 1][7]; // {-+, +-, --, ++, OS, SS, OS-SS}
  TH1D *arrayXi1DSBhigh[maxPtBins - 1][7]; // {-+, +-, --, ++, OS, SS, OS-SS}
  TH1D *arrayXi1Dcorrected[maxPtBins - 1][7]; // {-+, +-, --, ++, OS, SS, OS-SS}
  TH2D *arrayOm[maxPtBins - 1][7]; // {-+, +-, --, ++, OS, SS, OS-SS}
  double signsarray[4][2] = {{-1, 1}, {1, -1}, {-1, -1}, {1, 1}}; // sign combinations

  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    XiXidir->cd();
    int count = 0;
    for(auto signs : signsarray){
      TString triggerSign, assocSign;
      signs[0] < 0 ? triggerSign = "-" : triggerSign = "+";
      signs[1] < 0 ? assocSign = "-" : assocSign = "+";
      axranges aXiSig{{corr::signTrigg, {signs[0], signs[0]}}, {corr::signAssoc, {signs[1], signs[1]}},
                      {corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}},
                      {corr::dY, {-dYmax, dYmax}},
                      {corr::invMassTrigg, {sigXi[pTbin][0], sigXi[pTbin][1]}}, {corr::invMassAssoc, {sigXi[pTbin][0], sigXi[pTbin][1]}}
      };
      axranges aXiOmSig{{corr::signTrigg, {signs[0], signs[0]}}, {corr::signAssoc, {signs[1], signs[1]}},
                      {corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}},
                      {corr::dY, {-dYmax, dYmax}},
                      {corr::invMassTrigg, {sigXi[pTbin][0], sigXi[pTbin][1]}}, {corr::invMassAssoc, {sigOm[pTbin][0], sigOm[pTbin][1]}}
      };
      axranges aXiSBlow{{corr::signTrigg, {signs[0], signs[0]}}, {corr::signAssoc, {signs[1], signs[1]}},
                      {corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}},
                      {corr::dY, {-dYmax, dYmax}},
                      {corr::invMassTrigg, {SBlowXi[pTbin][0], SBlowXi[pTbin][1]}}, {corr::invMassAssoc, {sigXi[pTbin][0], sigXi[pTbin][1]}}
      };
      axranges aXiSBhigh{{corr::signTrigg, {signs[0], signs[0]}}, {corr::signAssoc, {signs[1], signs[1]}},
                      {corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}},
                      {corr::dY, {-dYmax, dYmax}},
                      {corr::invMassTrigg, {SBhighXi[pTbin][0], SBhighXi[pTbin][1]}}, {corr::invMassAssoc, {sigXi[pTbin][0], sigXi[pTbin][1]}}
      };

      // vvv trigger normalization vvv
      axranges aXiMass{{mass::pT, {pTbins[pTbin], pTbins[pTbin+1]}}, 
                      {mass::sign, {signs[0], signs[0]}},
                      {mass::invMass, {sigXi[0][0], sigXi[0][1]}}
      };
      TH1D *hMass = project(hEffCorrXiMass, mass::invMass, aXiMass);
      double nTriggers = hMass->Integral();
      cout << "n Xi triggers pTbin " << pTbin << " sign " << triggerSign << ": " << nTriggers << endl;
      hMass->Delete();
      // ^^^ trigger normalization ^^^

      // vvv XI-XI CORRELATIONS vvv
      // vvvvvvvvvvvvvvvvvvvvvvvvvv
      TH2D *h = new TH2D();
      TH1D *hdphi = new TH1D();
      h = project2D(hXiXi, corr::dY, corr::dPhi, aXiSig);
      h->SetName("hXi"+triggerSign+"Xi"+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
      h->RebinX(9);
      h->Divide(MEXiXiSig[pTbin]);
      h->SetStats(false);
      hdphi = h->ProjectionX("hXi"+triggerSign+"Xi"+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]+"_dphi");
      TH1D *hdphiRaw = (TH1D*) hdphi->Clone("hdphiRaw_"+triggerSign+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]); // make a raw copy before normalisation, to compare with sideband correlations
      hdphiRaw->Scale(1., "width");
      h->Scale(1. / nTriggers, "width"); // scale after projecting to get correct values (multiplying by width before projecting leads to wrong results)
      h->SetTitle("#Xi^{"+triggerSign+"} - #Xi^{"+assocSign+"} ("+pTlabels[pTbin]+" GeV/#it{c} < #it{p}_{T,trigger} < "+pTlabels[pTbin+1]+" GeV/#it{c});#Delta#varphi;#Delta y;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
      hdphi->Scale(1. / nTriggers, "width"); // scale after projecting to get correct values (multiplying by width before projecting leads to wrong results)
      hdphi->SetTitle("#Xi^{"+triggerSign+"} - #Xi^{"+assocSign+"} ("+pTlabels[pTbin]+" GeV/#it{c} < #it{p}_{T,trigger} < "+pTlabels[pTbin+1]+" GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
      hdphi->SetStats(false);
      arrayXi[pTbin][count] = h;
      arrayXi1D[pTbin][count] = hdphi;
      arrayXiRaw[pTbin][count] = hdphiRaw;
      
      // sidebands
      TH2D *hSBlow = new TH2D();
      TH1D *hdphiSBlow = new TH1D();
      hSBlow = project2D(hXiXi, corr::dY, corr::dPhi, aXiSBlow);;
      hSBlow->SetName("hSBlow_Xi"+triggerSign+"Xi"+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
      hSBlow->RebinX(9);
      hSBlow->Divide(MEXiXiSB[pTbin]);
      hSBlow->SetStats(false);
      hdphiSBlow = hSBlow->ProjectionX("hSBlow_Xi"+triggerSign+"Xi"+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]+"_dphi");
      hdphiSBlow->Scale(1., "width");
      hdphiSBlow->SetTitle("#Xi^{"+triggerSign+"} - #Xi^{"+assocSign+"} SB ("+pTlabels[pTbin]+" GeV/#it{c} < #it{p}_{T,trigger} < "+pTlabels[pTbin+1]+" GeV/#it{c});#Delta#varphi;d#it{N}/d(#Delta#varphi) (rad^{-1})");
      hdphiSBlow->SetStats(false);
      arrayXi1DSBlow[pTbin][count] = hdphiSBlow;

      TH2D *hSBhigh = new TH2D();
      TH1D *hdphiSBhigh = new TH1D();
      hSBhigh = project2D(hXiXi, corr::dY, corr::dPhi, aXiSBhigh);
      hSBhigh->SetName("hSBhigh_Xi"+triggerSign+"Xi"+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
      hSBhigh->RebinX(9);
      hSBhigh->Divide(MEXiXiSB[pTbin]);
      hSBhigh->SetStats(false);
      hdphiSBhigh = hSBhigh->ProjectionX("hSBhigh_Xi"+triggerSign+"Xi"+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]+"_dphi");
      hdphiSBhigh->Scale(1., "width");
      hdphiSBhigh->SetTitle("#Xi^{"+triggerSign+"} - #Xi^{"+assocSign+"} SB ("+pTlabels[pTbin]+" GeV/#it{c} < #it{p}_{T,trigger} < "+pTlabels[pTbin+1]+" GeV/#it{c});#Delta#varphi;d#it{N}/d(#Delta#varphi) (rad^{-1})");
      hdphiSBhigh->SetStats(false);
      arrayXi1DSBhigh[pTbin][count] = hdphiSBhigh;

      // sideband corrected correlations.
      TH1D *hCorrected = new TH1D(*hdphiRaw);
      double scale = signalWindow / (sidebandWindow[1] - sidebandWindow[0]); // scale factor to account for different widths of sideband and signal windows
      hCorrected->Add(hdphiSBlow, -scale);
      hCorrected->Add(hdphiSBhigh, -scale);
      hCorrected->SetName("hCorrected_Xi"+triggerSign+"Xi"+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
      hCorrected->SetTitle("#Xi^{"+triggerSign+"} - #Xi^{"+assocSign+"} SB corrected ("+pTlabels[pTbin]+" GeV/#it{c} < #it{p}_{T,trigger} < "+pTlabels[pTbin+1]+" GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
      hCorrected->Scale(1./(purityXi[pTbin] * nTriggers)); // already scaled by width
      hCorrected->SetStats(false);
      arrayXi1Dcorrected[pTbin][count] = hCorrected;

      // ^^^^^^^^^^^^^^^^^^^^^^^^^^
      // ^^^ XI-XI CORRELATIONS ^^^

      TH2D *hOm = new TH2D();
      hOm = project2D(hXiOm, corr::dY, corr::dPhi, aXiOmSig);
      hOm->SetName("hXi"+triggerSign+"Om"+assocSign+"pT_"+pTlabels[pTbin]+"_"+pTlabels[pTbin+1]);
      hOm->RebinX(9);
      hOm->Scale(1. / nTriggers); 
      hOm->Divide(hMEXiXi2D); // TODO: XI-OM ME
      hOm->SetTitle("#Xi^{"+triggerSign+"} - #Om^{"+assocSign+"} ("+pTlabels[pTbin]+" GeV/#it{c} < #it{p}_{T,trigger} < "+pTlabels[pTbin+1]+" GeV/#it{c});#Delta#varphi;#Delta y;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
      hOm->SetStats(false);
      arrayOm[pTbin][count] = hOm;
      count++;
    }
    // vvv DO SUBTRACTIONS AND AVERAGES vvv
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    doSubAndAvg(arrayXi[pTbin], "XiXi2D", pTbin, pTlabels[pTbin], pTlabels[pTbin+1]);
    doSubAndAvg(arrayXi1D[pTbin], "XiXidphi", pTbin, pTlabels[pTbin], pTlabels[pTbin+1]);
    doSubAndAvg(arrayXi1DSBlow[pTbin], "XiXidphiSBlow", pTbin, pTlabels[pTbin], pTlabels[pTbin+1]);
    doSubAndAvg(arrayXi1DSBhigh[pTbin], "XiXidphiSBhigh", pTbin, pTlabels[pTbin], pTlabels[pTbin+1]);
    doSubAndAvg(arrayXi1Dcorrected[pTbin], "XiXidphiCorrected", pTbin, pTlabels[pTbin], pTlabels[pTbin+1]);
    doSubAndAvg(arrayXiRaw[pTbin], "XiXidphiRaw", pTbin, pTlabels[pTbin], pTlabels[pTbin+1]);
    doSubAndAvg(arrayOm[pTbin], "XiOm2D", pTbin, pTlabels[pTbin], pTlabels[pTbin+1]);
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    // ^^^ DO SUBTRACTIONS AND AVERAGES ^^^

    if(makePDF){
      for (TH2D* h : arrayXi[pTbin]){
        if(!h) {
          cout << "found null histogram in array while drawing, skipping it..." << endl;
          continue;
        }
        h->GetZaxis()->SetTitleOffset(1.5);
        h->Draw("surf1");
        TString hname = h->GetName();

        float x = 0.10;
        float y = 0.90;
        float spacing = 0.05;
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.DrawLatex(x, y, h->GetTitle());
        latex.DrawLatex(x, y - spacing, "1.0 GeV/#it{c} < #it{p}_{T,assoc} < #it{p}_{T,trigger}");
        latex.DrawLatex(x, y - 2*spacing, "pp, #sqrt{s} = 13 TeV");

        h->SetTitle("");
        c->Print("figures/corr/"+hname+".pdf");
        c->Clear();
      }
      for (TH1D* h : arrayXi1D[pTbin]){
        if(!h) {
          cout << "found null histogram in array while drawing, skipping it..." << endl;
          continue;
        }
        h->Draw();
        TString hname = h->GetName();

        float x = 0.50;
        float y = 0.90;
        float spacing = 0.05;
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.DrawLatex(x, y, h->GetTitle());
        latex.DrawLatex(x, y - spacing, "1.0 GeV/#it{c} < #it{p}_{T,assoc} < #it{p}_{T,trigger}");
        latex.DrawLatex(x, y - 2*spacing, "pp, #sqrt{s} = 13 TeV");

        h->SetTitle("");
        // c->GetPad(0)->SetLeftMargin(0.15);
        c->Print("figures/corr/"+hname+".pdf");
        c->Clear();
      }
    }
  }

  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  return 0;
}