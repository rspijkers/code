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

// enum for axisnumbers/names:
struct corr{
  enum{
    dPhi, 
    dY, 
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


// compare with run 2
void plotRun2(){
  std::vector<double> datapoints = {0.021814, 0.0223102, 0.0368555, 0.0459319, 0.0563678, 0.0604721, 0.0484733, 0.0334666, 0.0292751, 0.019142, 0.0122185, 0.0153656, 0.010596, 0.0109022, 0.0129868, 0.0104814, 0.00605552, 0.0126087, 0.0124405, 0.0161898};
  std::vector<double> staterrors = {0.00363465, 0.0041518, 0.00355898, 0.00356545, 0.003877, 0.0037702, 0.0035463, 0.00385414, 0.00379416, 0.00342748, 0.00335034, 0.00348166, 0.00373101, 0.00344667, 0.00506973, 0.00411657, 0.00352012, 0.00346011, 0.00344019, 0.00342756};
}

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

// TODO: write function that gets the correct efficiency map.
//    access the dpl-config.json in the folder of the corresponding train number
//    check the path of the efficiency map
//    connect to the ccdb
//    download/access the efficiency map(s).
//
// OKAY NEVER MIND just do efficiency correction on the fly in O2Physics and don't look at it again
void getEfficiencyMaps(TString trainnr){
//   TGrid *grid = TGrid::Connect("alien://");
//   TString jsonpath("results/" + trainnr + "/dpl-config.json"); 
//   std::ifstream f(jsonpath);
//   json test = json::parse(f);

//   // this piece of magic gets the configuration for the cascade-correlations task, 
//   // then the value of the corresponding configurable (efficiencyCCDBPath) as a string
//   TString effPath(test.at("cascade-correlations")["efficiencyCCDBPath"].get<std::string>());
//   cout << effPath << endl;

//   // use TFile::Open for remote files
//   // TFile *effFile = TFile::Open("alien:///alice/data/CCDB/" + effPath + "/CascadeEfficienciesPID.root", "READ");

//   cout << grid->Ls("alien:///alice/data/CCDB/" + effPath) << endl;

//   UInt_t size = grid->Ls("alien:///alice/data/CCDB/" + effPath)->GetFileInfoList()->GetSize();
//   for (UInt_t i = 0; i < size; i++){
//     TString subdir(grid->Ls("alien:///alice/data/CCDB/" + effPath)->GetPath(i));
//     cout << subdir << endl;
//     grid->Ls("alien:///alice/data/CCDB/" + effPath)
//   }
//   for (auto const& test : *(grid->Ls("alien:///alice/data/CCDB/" + effPath)->GetFileInfoList())){
//     cout << test->ClassName() << endl;
//   }

//   // for (auto const& test : *(grid->Ls("alien:///alice/data/CCDB/" + effPath))){
//   //   cout << test->GetPath() << endl;
//   // }

//   // "Users\/r\/rspijker\/test\/Efftest"
}

void doXiInvMassFits(bool makePDF = false){
  TDirectory* XiInvMass = outputFile->mkdir("XiInvMass");
  XiInvMass->cd();

  TF1 *f1 = new TF1("f1", "pol2(0) + gaus(3) + gaus(6)", 1.29, 1.42); 
  // parameters are: pol{0,1,2}, A1, mu1, sigma1, A2, mu2, sigma2
  f1->SetParameters(0, 0, 10, 0, 1.321, 0.005, 0, 1.321, 0.001); 
  f1->SetParLimits(4, 1.31, 1.33);
  f1->SetParLimits(7, 1.31, 1.33);
  f1->SetParLimits(5, 0, 0.01);
  f1->SetParLimits(8, 0, 0.01);

  TF1 *fBKGXi = new TF1("fBKGXi", pol2bkgXi, 1.29, 1.42, 3);
  TH2F *hMassXiMinus = inputDir->Get<TH2F>("hMassXiMinus");
  hMassXiMinus->SetDirectory(outputFile);
  cout << "Start Xi mass fitting..." << endl;
  TCanvas *c = new TCanvas();
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hMassXiMinus->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin + 1]);
    TH1D *h = hMassXiMinus->ProjectionX("hMassXiMinus_"+ pTlabels[pTbin + 1]);
    h->GetXaxis()->SetRangeUser(1.2, 1.5);
    h->SetTitle("#Xi^{-} inv. mass for " + pTlabels[pTbin] + " GeV < p_{T} < " + pTlabels[pTbin+1] + " GeV");
    f1->SetParameter(3, .4* h->GetMaximum());
    f1->SetParameter(6, .4* h->GetMaximum());
    h->Fit(fBKGXi, "SLBQRO");
    f1->SetParameter(0, fBKGXi->GetParameter(0));
    f1->SetParameter(1, fBKGXi->GetParameter(1));
    f1->SetParameter(2, fBKGXi->GetParameter(2));

    TFitResultPtr r1 = h->Fit("f1", "SLBQR", "", 1.29, 1.42);
    if (r1->Chi2() > 1) cout << "Warning: Fit Chi2(" << r1->Chi2() << ") > 1 in pT bin " << pTbin << endl;

    double mu = (f1->GetParameter(4) + f1->GetParameter(7)) / 2.;
    double sigma = (f1->GetParameter(5) + f1->GetParameter(8)) / 2.;

    h->GetXaxis()->SetRangeUser(1.28, 1.38);
    h->SetStats(kFALSE);
    // h->SetLineWidth(3);
    h->Draw();
    if(makePDF) c->Print("figures/Ximass" + pTlabels[pTbin] + "_" + pTlabels[pTbin+1] + ".pdf");
    c->Clear();

    sigXi[pTbin][0] = mu - 3*sigma; sigXi[pTbin][1] = mu + 3*sigma;
    bkgXi[pTbin][0] = mu + 4*sigma; bkgXi[pTbin][1] = mu + 10*sigma;
    // cout << sigXi[pTbin][0] << " " << sigXi[pTbin][1] << " " << bkgXi[pTbin][0] << " " << bkgXi[pTbin][1] << endl;
  }
  outputFile->cd();
}

void doOmInvMassFits(){
  TDirectory* OmInvMass = outputFile->mkdir("OmInvMass");
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

  TH2F *hMassOmegaMinus = inputDir->Get<TH2F>("hMassOmegaMinus");
  hMassOmegaMinus->SetDirectory(outputFile);
  cout << "Start Omega mass fitting..." << endl;
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hMassOmegaMinus->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin+1]);
    TH1D *h = hMassOmegaMinus->ProjectionX("hMassOmegaMinus_"+ pTlabels[pTbin + 1]);
    h->GetXaxis()->SetRangeUser(1.6, 1.9);
    h->SetTitle("#Omega^{-} inv. mass for " + pTlabels[pTbin] + " GeV < p_{T} < " + pTlabels[pTbin+1] + " GeV");
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
    h->SetLineWidth(3);

    if(mu + 10*sigma > 1.76) cout << "Warning: Bkg region exceeds limit of 1.76 GeV in pT bin " << pTbin << endl;
    sigOm[pTbin][0] = mu - 3*sigma; sigOm[pTbin][1] = mu + 3*sigma;
    bkgOm[pTbin][0] = mu + 4*sigma; bkgOm[pTbin][1] = mu + 10*sigma;
    // cout << sigOm[pTbin][0] << " " << sigOm[pTbin][1] << " " << bkgOm[pTbin][0] << " " << bkgOm[pTbin][1] << endl;
  }
  outputFile->cd();
}

std::vector<TH1D*> doMixedEvents(TString trainnr, TString filename = "AnalysisResults.root"){
  axranges aMEXiXiOS{{corr::ptTrigg, {1.0, 10.}}, {corr::ptAssoc, {1.0, 10.}}, 
                     {corr::dY, {-1., 1.}},
                     {corr::invMassTrigg, {1.31, 1.33}}, {corr::invMassAssoc, {1.31, 1.33}} // fixme read invmass from fits (dedicated function?)
  };
  // do ME here
  THnSparse *hXiXiOS, *hXiXiSS;
  inputDir->GetObject("MixedEvents/hMEXiXiOS", hXiXiOS);
  inputDir->GetObject("MixedEvents/hMEXiXiSS", hXiXiSS);
  THnSparse *hOmXiOS, *hOmXiSS;
  inputDir->GetObject("MixedEvents/hMEOmXiOS", hOmXiOS);
  inputDir->GetObject("MixedEvents/hMEOmXiSS", hOmXiSS);
  THnSparse *hOmOmOS, *hOmOmSS;
  inputDir->GetObject("MixedEvents/hMEOmOmOS", hOmOmOS);
  inputDir->GetObject("MixedEvents/hMEOmOmSS", hOmOmSS);

  // Okay let's take Xi-Xi OS combination to make ME projection in rapidity and phi (all pT ranges)
  TH1D *hXiXiOSdY = project(hXiXiOS, corr::dY, aMEXiXiOS);
  hXiXiOSdY->SetName("hMEXiOSdY");
  hXiXiOSdY->SetTitle("ME #Xi-#Xi OS (pT integrated)");
  hXiXiOSdY->SetYTitle("pairs OS");
  hXiXiOSdY->SetLineWidth(3);
  hXiXiOSdY->Scale(1. / hXiXiOSdY->GetMaximum());

  TH1D *hXiXiSSdY = project(hXiXiSS, corr::dY, aMEXiXiOS);
  hXiXiSSdY->SetName("hMEXiSSdY");
  hXiXiSSdY->SetTitle("ME #Xi-#Xi SS (pT integrated)");
  hXiXiSSdY->SetYTitle("pairs SS");
  hXiXiSSdY->SetLineWidth(3);
  hXiXiSSdY->Scale(1. / hXiXiSSdY->GetMaximum());

  TH1D *hXiXiOSdPhi = project(hXiXiOS, corr::dPhi, aMEXiXiOS);
  hXiXiOSdPhi->SetName("hMEXiOSdPhi");
  hXiXiOSdPhi->SetTitle("ME #Xi-#Xi OS (pT integrated)");
  hXiXiOSdPhi->SetYTitle("pairs OS");
  hXiXiOSdPhi->SetLineWidth(3);
  hXiXiOSdPhi->Rebin(9);
  hXiXiOSdPhi->Scale(1. / hXiXiOSdPhi->GetMaximum());

  TH1D *hXiXiSSdPhi = project(hXiXiSS, corr::dPhi, aMEXiXiOS);
  hXiXiSSdPhi->SetName("hMEXiSSdPhi");
  hXiXiSSdPhi->SetTitle("ME #Xi-#Xi SS (pT integrated)");
  hXiXiSSdPhi->SetYTitle("pairs SS");
  hXiXiSSdPhi->SetLineWidth(3);
  hXiXiSSdPhi->Rebin(9);
  hXiXiSSdPhi->Scale(1. / hXiXiSSdPhi->GetMaximum());

  axranges aMEOmXiOS{{corr::ptTrigg, {1., 5.0}}, {corr::ptAssoc, {1., 5.0}},
                    //  {corr::dY, {-1., 1.}},
                      {corr::invMassTrigg, {1.66, 1.685}}, {corr::invMassAssoc, {1.31, 1.33}}
  };
  // Okay let's take Om-Xi OS combination to make ME projection in rapidity and phi (all pT ranges)
  TH1D *hOmXiOSdY = project(hOmXiOS, corr::dY, aMEOmXiOS);
  hOmXiOSdY->SetName("hMEOmXiOSdY");
  hOmXiOSdY->SetTitle("ME #Omega-#Xi OS (pT integrated)");
  hOmXiOSdY->SetYTitle("pairs OS");
  hOmXiOSdY->SetLineWidth(3);
  hOmXiOSdY->Scale(1. / hOmXiOSdY->GetMaximum());

  TH1D *hOmXiSSdY = project(hOmXiSS, corr::dY, aMEOmXiOS);
  hOmXiSSdY->SetName("hMEOmXiSSdY");
  hOmXiSSdY->SetTitle("ME #Omega-#Xi SS (pT integrated)");
  hOmXiSSdY->SetYTitle("pairs SS");
  hOmXiSSdY->SetLineWidth(3);
  hOmXiSSdY->Scale(1. / hOmXiSSdY->GetMaximum());

  TH1D *hOmXiOSdPhi = project(hOmXiOS, corr::dPhi, aMEOmXiOS);
  hOmXiOSdPhi->SetName("hMEOmXiOSdPhi");
  hOmXiOSdPhi->SetTitle("ME #Omega-#Xi OS (pT integrated)");
  hOmXiOSdPhi->SetYTitle("pairs OS");
  hOmXiOSdPhi->SetLineWidth(3);
  hOmXiOSdPhi->Rebin(9);
  hOmXiOSdPhi->Scale(1. / hOmXiOSdPhi->GetMaximum());

  TH1D *hOmXiSSdPhi = project(hOmXiSS, corr::dPhi, aMEOmXiOS);
  hOmXiSSdPhi->SetName("hMEOmXiSSdPhi");
  hOmXiSSdPhi->SetTitle("ME #Omega-#Xi SS (pT integrated)");
  hOmXiSSdPhi->SetYTitle("pairs SS");
  hOmXiSSdPhi->SetLineWidth(3);
  hOmXiSSdPhi->Rebin(9);
  hOmXiSSdPhi->Scale(1. / hOmXiSSdPhi->GetMaximum());
  return {hXiXiOSdPhi, hXiXiSSdPhi, hXiXiOSdY, hXiXiSSdY, hOmXiOSdPhi, hOmXiSSdPhi, hOmXiOSdY, hOmXiSSdY};
}

int postprocessingResults(TString trainnr, TString filename = "AnalysisResults.root", bool makePDF = false) {
  TH1::SetDefaultSumw2(); // Make sure we propagate the errors
  gStyle->SetHistLineWidth(3);
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

  // QA plots
  cout << "-- QA PROJECTIONS --" << endl;
  TH1F *hPhi = inputDir->Get<TH1F>("hPhi");
  hPhi->SetDirectory(outputFile);
  hPhi->Draw();

  doQAprojections(inputFile);

  // test ME
  cout << "-- MIXED EVENTS --" << endl;
  std::vector<TH1D*> hME_vector;
  hME_vector = doMixedEvents(trainnr);

  // test inv mass fitting
  cout << "-- INV MASS FITS --" << endl;
  doXiInvMassFits();
  doOmInvMassFits();
  // todo: assert some properties of the inv mass region boundaries to ensure everything went ok.

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
  // // todo remove these couts
  // for (auto i : NTrigXi){
  //   cout << "test Xi " << i << endl;
  // }
  // for (auto i : NTrigOm){
  //   cout << "test Om " << i << endl;
  // }
  
  // Load the THnSparses, 8 in total (4 particle combinations * 2 sign combinations)
  THnSparse *hXiXiOS, *hXiXiSS;
  inputDir->GetObject("hXiXiOS", hXiXiOS);
  inputDir->GetObject("hXiXiSS", hXiXiSS);
  THnSparse *hOmXiOS, *hOmXiSS;
  inputDir->GetObject("hOmXiOS", hOmXiOS);
  inputDir->GetObject("hOmXiSS", hOmXiSS);
  THnSparse *hOmOmOS, *hOmOmSS;
  inputDir->GetObject("hOmOmOS", hOmOmOS);
  inputDir->GetObject("hOmOmSS", hOmOmSS);

  // test projection of Xi-Xi OS onto 2D dphi-dy
  axranges aPtIntegrated{{corr::ptTrigg, {1., 10}}, {corr::ptAssoc, {1., 10}},
                         {corr::dY, {-1., 1.}},
                         {corr::invMassTrigg, {1.31, 1.33}}, {corr::invMassAssoc, {1.31, 1.33}}
  };
  axranges aPtIntMass{{mass::pT, {1., pTmax}}, 
                      // {mass::y, {-0.5, 0.5}},
                      {mass::invMass, {1.31, 1.33}}
  };
  // make single Xi spectra (charge independent)
  // get Nevents:
  double Nevents = inputFile->Get<TH1D>("event-selection-task/hColCounterAcc")->Integral();
  cout << Nevents << endl;
  TH1D* hXiSpectra = project(hEffCorrXiMass, mass::pT, {{mass::pT, {0., pTmax}}, {mass::y, {-0.5, 0.5}}, {mass::invMass, {1.31, 1.33}}});
  hXiSpectra->SetName("hXiSpectra");
  Double_t Xibinning[14] = {0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5};
  TH1D* hSpectraRebin = (TH1D*) hXiSpectra->Rebin(13, "hSpectraRebin", Xibinning);
  hSpectraRebin->Scale(1./(Nevents), "width"); // divide by an extra factor of 2 because of Xi charges ---- 65536 = factor of 2^n used in run 3 AN
  TH2F *hMassXiMinus = inputDir->Get<TH2F>("hMassXiMinus");
  TH1D *hXiMinusMass1D = hMassXiMinus->ProjectionX();
  hXiMinusMass1D->SetName("hXiMinusMass1D");
  TH1D *hRawXiMass = hMassXiMinus->ProjectionY("hRawXiMass");
  TH1D* hRawSpectraRebin = (TH1D*) hRawXiMass->Rebin(13, "hRawSpectraRebin", Xibinning);
  hRawSpectraRebin->Scale(1./Nevents, "width"); // ---- 65536 = factor of 2^n used in run 3 AN

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

  TH1D *hXiMassPtInt = project(hEffCorrXiMass, mass::invMass, aPtIntMass);
  hXiMassPtInt->SetName("hEffXiMass_pTint");
  double nXiTriggers = hXiMassPtInt->Integral();

  // TEMP remove this
  // make quick pT spectrum of eff corrected xi, to make ratio
  TH1D *hXiEffpT = project(hEffCorrXiMass, mass::pT, aPtIntMass);
  hXiEffpT->SetName("TEMPXiEffpT");

  TH2D *h2DtestOS = project2D(hXiXiOS, corr::dPhi, corr::dY, aPtIntegrated);
  h2DtestOS->SetName("h2DtestOS");
  h2DtestOS->Scale(1. / nXiTriggers);
  // h2DtestOS->Rebin2D(1,5);
  TH2D *h2DtestSS = project2D(hXiXiSS, corr::dPhi, corr::dY, aPtIntegrated);
  h2DtestSS->SetName("h2DtestSS");
  h2DtestSS->Scale(1. / nXiTriggers);
  // h2DtestSS->Rebin2D(1,5);
  TH2D *h2Dtest = new TH2D(*h2DtestOS);
  h2Dtest->Add(h2DtestOS, h2DtestSS, 1, -1);
  h2Dtest->SetName("h2Dtest");
  h2Dtest->ProjectionX("hDYtest");
  h2Dtest->ProjectionY("hDPhitest");

  // inv mass of Xi trig & assoc (OS)
  TH1D *hXiInvMassTrig = project(hXiXiOS, corr::invMassTrigg, aPtIntegrated);
  hXiInvMassTrig->SetName("hXiInvMassTrig");
  TH1D *hXiInvMassAssoc = project(hXiXiOS, corr::invMassAssoc, aPtIntegrated);
  hXiInvMassAssoc->SetName("hXiInvMassAssoc");

  // quick pT integrated ME correction applied
  TH1D* hDYOScorrected = project(hXiXiOS, corr::dY, aPtIntegrated);
  hDYOScorrected->SetName("hDYOScorrected");
  hDYOScorrected->Scale(1. / nXiTriggers);
  hDYOScorrected->Divide(hME_vector[2]);
  TH1D* hDYSScorrected = project(hXiXiSS, corr::dY, aPtIntegrated);
  hDYSScorrected->SetName("hDYSScorrected");
  hDYSScorrected->Scale(1. / nXiTriggers);
  hDYSScorrected->Divide(hME_vector[3]);
  TH1D *hDYcorrected = new TH1D(*hDYOScorrected);
  hDYcorrected->Add(hDYOScorrected, hDYSScorrected, 1, -1);
  hDYcorrected->SetName("hDYcorrected");
  
  TH1D* hDPhiOScorrected = project(hXiXiOS, corr::dPhi, aPtIntegrated);
  hDPhiOScorrected->SetName("hDPhiOScorrected");
  hDPhiOScorrected->Scale(1. / nXiTriggers);
  hDPhiOScorrected->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  hDPhiOScorrected->Rebin(9);
  hDPhiOScorrected->Divide(hME_vector[0]);
  hDPhiOScorrected->SetTitle("OS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hDPhiOScorrected->SetStats(false);
  TH1D* hDPhiSScorrected = project(hXiXiSS, corr::dPhi, aPtIntegrated);
  hDPhiSScorrected->SetName("hDPhiSScorrected");
  hDPhiSScorrected->Scale(1. / nXiTriggers);
  hDPhiSScorrected->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  hDPhiSScorrected->Rebin(9);
  hDPhiSScorrected->Divide(hME_vector[1]);
  hDPhiSScorrected->SetTitle("SS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hDPhiSScorrected->SetStats(false);
  TH1D *hDPhicorrected = new TH1D(*hDPhiOScorrected);
  hDPhicorrected->Add(hDPhiOScorrected, hDPhiSScorrected, 1, -1);
  hDPhicorrected->SetName("hDPhicorrected");
  hDPhicorrected->SetTitle("OS - SS #Xi - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 10 GeV/#it{c});#Delta#varphi;1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hDPhicorrected->SetStats(false);

  ///// RUN 2 COMPARISONS

  const double datapointsOS[21] = {0, 0.00772392, 0.00763424, 0.0102429, 0.00970434, 0.0147818, 0.0131757, 0.0114212, 0.00877517, 0.00667961, 0.00695319, 0.00729105, 0.00602286, 0.00449109, 0.00573827, 0.00629415, 0.00505972, 0.00429942, 0.00724581, 0.00570471, 0.00750361};
  const double staterrorsOS[21] = {0, 0.00140528, 0.00105404, 0.00141499, 0.000748423, 0.00128475, 0.00101308, 0.00136926, 0.00085206, 0.000793002, 0.000902522, 0.000955686, 0.000843041, 0.000948397, 0.000851573, 0.000820212, 0.00101491, 0.00104776, 0.000878481, 0.000767664, 0.00100554};
  TH1D* hRun2OS = new TH1D(*hDPhiOScorrected);
  hRun2OS->SetName("hRun2OS");
  hRun2OS->SetTitle("Run 2 OS #Xi - #Xi");
  hRun2OS->SetContent(datapointsOS);
  hRun2OS->SetError(staterrorsOS);
  hRun2OS->SetStats(false);
  hRun2OS->SetLineColor(kRed);
  TCanvas *c = new TCanvas("c");
  c->cd();
  TRatioPlot *hRun2OSRatio = new TRatioPlot(hDPhiOScorrected, hRun2OS, "divsym");
  hRun2OSRatio->SetH1DrawOpt("E");
  hRun2OSRatio->Draw("nogrid");
  hRun2OSRatio->GetUpperRefYaxis()->SetRangeUser(0., 0.02);
  hRun2OSRatio->GetLowerRefYaxis()->SetRangeUser(0., 1.);
  hRun2OSRatio->GetLowYaxis()->SetNdivisions(505);
  TLegend *legendOS = new TLegend(0.55, 0.75, 0.75, 0.85);
  hRun2OSRatio->GetUpperPad()->cd(); // draw legend in upper pad
  legendOS->AddEntry(hDPhiOScorrected, "this analysis");
  legendOS->AddEntry(hRun2OS, "run 2");
  legendOS->Draw();
  c->Write("cOSRatio");
  if(makePDF) c->Print("figures/OSRatio.pdf");
  c->Clear();
  
  const double datapointsSS[21] = {0, 0.00314547, 0.00574206, 0.0049273, 0.00321404, 0.00306671, 0.00141205, 0.00276105, 0.00273416, 0.00406806, 0.00376302, 0.00390338, 0.0048735, 0.00529755, 0.00359733, 0.0080713, 0.00588222, 0.00452432, 0.00337606, 0.00428579, 0.00485876};
  const double staterrorsSS[21] = {0, 0.000688171, 0.00119578, 0.000846638, 0.00051905, 0.000517229, 0.00148369, 0.000717417, 0.000644928, 0.000621341, 0.000561713, 0.000813817, 0.00107704, 0.000920041, 0.00073849, 0.00245958, 0.000746747, 0.000804554, 0.000785611, 0.000791499, 0.000860004};
  TH1D* hRun2SS = new TH1D(*hDPhiSScorrected);
  hRun2SS->SetName("hRun2SS");
  hRun2SS->SetTitle("Run 2 SS #Xi - #Xi");
  hRun2SS->SetContent(datapointsSS);
  hRun2SS->SetError(staterrorsSS);
  hRun2SS->SetStats(false);
  hRun2SS->SetLineColor(kRed);
  TRatioPlot *hRun2SSRatio = new TRatioPlot(hDPhiSScorrected, hRun2SS, "divsym");
  hRun2SSRatio->SetH1DrawOpt("E");
  hRun2SSRatio->Draw("nogrid");
  hRun2SSRatio->GetUpperRefYaxis()->SetRangeUser(0., 0.01);
  hRun2SSRatio->GetLowerRefYaxis()->SetRangeUser(0., 1.);
  hRun2SSRatio->GetLowYaxis()->SetNdivisions(505);
  TLegend *legendSS = new TLegend(0.2, 0.75, 0.4, 0.85);
  hRun2SSRatio->GetUpperPad()->cd(); // draw legend in upper pad
  legendSS->AddEntry(hDPhiSScorrected, "this analysis");
  legendSS->AddEntry(hRun2SS, "run 2");
  legendSS->Draw();
  c->Write("cSSRatio");
  if(makePDF) c->Print("figures/SSRatio.pdf");
  c->Clear();
  
  const double datapoints[21] = {0, 0.00457845, 0.00189219, 0.00531563, 0.00649031, 0.0117151, 0.0117636, 0.00866012, 0.00604101, 0.00261155, 0.00319017, 0.00338767, 0.00114936, -0.000806458, 0.00214094, -0.00177715, -0.000822495, -0.000224905, 0.00386975, 0.00141892, 0.00264486};
  const double staterrors[21] = {0, 0.00156473, 0.00159401, 0.00164894, 0.000910796, 0.00138496, 0.00179657, 0.00154582, 0.00106861, 0.00100743, 0.00106305, 0.00125524, 0.00136775, 0.00132134, 0.00112718, 0.00259274, 0.00126003, 0.00132103, 0.00117852, 0.00110262, 0.00132315};
  TH1D* hRun2 = new TH1D(*hDPhicorrected);
  hRun2->SetName("hRun2");
  hRun2->SetTitle("Run 2 #Xi - #Xi");
  hRun2->SetContent(datapoints);
  hRun2->SetError(staterrors);
  hRun2->SetStats(false);
  hRun2->SetLineColor(kRed);
  TRatioPlot *hRun2Ratio = new TRatioPlot(hDPhicorrected, hRun2, "divsym");
  // TH1D* hRun2Ratio = new TH1D(*hDPhicorrected);
  hRun2Ratio->SetH1DrawOpt("E");
  hRun2Ratio->Draw("nogrid");
  hRun2Ratio->GetUpperRefYaxis()->SetRangeUser(0., 0.015);
  hRun2Ratio->GetLowerRefYaxis()->SetRangeUser(0., 1.);
  hRun2Ratio->GetLowYaxis()->SetNdivisions(505);
  TLegend *legend2 = new TLegend(0.55, 0.75, 0.75, 0.85);
  hRun2Ratio->GetUpperPad()->cd(); // draw legend in upper pad
  legend2->AddEntry(hDPhicorrected, "this analysis");
  legend2->AddEntry(hRun2, "run 2");
  legend2->Draw();
  c->Write("cRun2Ratio");
  if(makePDF) c->Print("figures/Run2Ratio.pdf");
  c->Clear();

  ///// END RUN 2 COMPARISONS

  cout << "Xi-Xi yield " << hDPhicorrected->Integral() << endl;

  // test projection of Xi-Xi OS onto 2D dphi-dy
  axranges aOmPtIntegrated{{corr::ptTrigg, {1., 5.0}}, {corr::ptAssoc, {1., 5.0}},
                        //  {corr::dY, {-1., 1.}},
                         {corr::invMassTrigg, {1.66, 1.685}}, {corr::invMassAssoc, {1.31, 1.33}}
  };
  axranges aOmPtIntMass{{mass::pT, {1., 5.0}}, 
                      // {mass::y, {-0.5, 0.5}},
                      {mass::invMass, {1.66, 1.685}}
  };

  TH1D *hOmMassPtInt = project(hEffCorrOmegaMass, mass::invMass, aOmPtIntMass);
  hOmMassPtInt->SetName("hEffOmMass_pTint");
  double nOmTriggers = hOmMassPtInt->Integral();

  TH2D *h2DOmXiOS = project2D(hOmXiOS, corr::dPhi, corr::dY, aOmPtIntegrated);
  h2DOmXiOS->SetName("h2DOmXiOS");
  h2DOmXiOS->Scale(1. / nOmTriggers);
  h2DOmXiOS->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  // h2DOmXiOS->Rebin2D(1,5);
  TH2D *h2DOmXiSS = project2D(hOmXiSS, corr::dPhi, corr::dY, aOmPtIntegrated);
  h2DOmXiSS->SetName("h2DOmXiSS");
  h2DOmXiSS->Scale(1. / nOmTriggers);
  h2DOmXiSS->Scale(20. / (2 * M_PI)); // scale by bin width for correct dN/d(dphi)
  // h2DOmXiSS->Rebin2D(1,5);
  TH2D *h2DOmXi = new TH2D(*h2DOmXiOS);
  h2DOmXi->Add(h2DOmXiOS, h2DOmXiSS, 1, -1);
  h2DOmXi->SetName("h2DOmXi");
  h2DOmXi->ProjectionX("hDYOmXi");
  h2DOmXi->ProjectionY("hDPhiOmXi");

  // quick pT integrated ME correction applied
  TH1D* hDYOSOmXicorrected = project(hOmXiOS, corr::dY, aOmPtIntegrated);
  hDYOSOmXicorrected->SetName("hDYOSOmXicorrected");
  hDYOSOmXicorrected->Scale(1. / nOmTriggers);
  hDYOSOmXicorrected->Divide(hME_vector[2]);
  TH1D* hDYSSOmXicorrected = project(hOmXiSS, corr::dY, aOmPtIntegrated);
  hDYSSOmXicorrected->SetName("hDYSSOmXicorrected");
  hDYSSOmXicorrected->Scale(1. / nOmTriggers);
  hDYSSOmXicorrected->Divide(hME_vector[3]);
  TH1D *hDYOmXicorrected = new TH1D(*hDYOSOmXicorrected);
  hDYOmXicorrected->Add(hDYOSOmXicorrected, hDYSSOmXicorrected, 1, -1);
  hDYOmXicorrected->SetName("hDYOmXicorrected");
  
  TH1D* hDPhiOSOmXicorrected = project(hOmXiOS, corr::dPhi, aOmPtIntegrated);
  hDPhiOSOmXicorrected->SetName("hDPhiOSOmXicorrected");
  hDPhiOSOmXicorrected->Scale(1. / nOmTriggers);
  hDPhiOSOmXicorrected->Rebin(9);
  hDPhiOSOmXicorrected->Divide(hME_vector[0]);
  TH1D* hDPhiSSOmXicorrected = project(hOmXiSS, corr::dPhi, aOmPtIntegrated);
  hDPhiSSOmXicorrected->SetName("hDPhiSSOmXicorrected");
  hDPhiSSOmXicorrected->Scale(1. / nOmTriggers);
  hDPhiSSOmXicorrected->Rebin(9);
  hDPhiSSOmXicorrected->Divide(hME_vector[1]);
  TH1D *hDPhiOmXicorrected = new TH1D(*hDPhiOSOmXicorrected);
  hDPhiOmXicorrected->Add(hDPhiOSOmXicorrected, hDPhiSSOmXicorrected, 1, -1);
  hDPhiOmXicorrected->SetName("hDPhiOmXicorrected");
  hDPhiOmXicorrected->SetTitle("OS - SS #Omega - #Xi (1 GeV/#it{c} < #it{p}_{T,trigger} < 5 GeV/#it{c})");
  hDPhiOmXicorrected->GetYaxis()->SetTitle("1/#it{N}_{triggers} d#it{N}/d(#Delta#varphi) (rad^{-1})");
  hDPhiOmXicorrected->SetStats(false);
  cout << "Om-Xi yield " << hDPhiOmXicorrected->Integral() << endl;

  // Let's do inv-mass projections
  axranges massXi{};
  axranges massBoth{};
  axranges massOm{};
  TH1D *hInvMassXiXi = project(hXiXiOS, corr::invMassTrigg, massXi);
  hInvMassXiXi->SetName("hInvMassXiXi");
  hInvMassXiXi->SetTitle("inv Mass of Xi (bachelor == pion and != kaon)");
  TH1D *hInvMassXiOm = project(hXiXiOS, corr::invMassTrigg, massOm);
  hInvMassXiOm->SetName("hInvMassXiOm");
  hInvMassXiOm->SetTitle("inv Mass of Xi (bachelor == kaon and != pion)");
  TH1D *hInvMassXiBoth = project(hXiXiOS, corr::invMassTrigg, massBoth);
  hInvMassXiBoth->SetName("hInvMassXiBoth");
  hInvMassXiBoth->SetTitle("inv Mass of Xi (bachelor consistent with both pion, kaon)");

  TH1D *hInvMassOmXi = project(hOmXiOS, corr::invMassTrigg, massXi);
  hInvMassOmXi->SetName("hInvMassOmXi");
  hInvMassOmXi->SetTitle("inv Mass of Omega (bachelor == pion and != kaon)");
  TH1D *hInvMassOmOm = project(hOmXiOS, corr::invMassTrigg, massOm);
  hInvMassOmOm->SetName("hInvMassOmOm");
  hInvMassOmOm->SetTitle("inv Mass of Omega (bachelor == kaon and != pion)");
  TH1D *hInvMassOmBoth = project(hOmXiOS, corr::invMassTrigg, massBoth);
  hInvMassOmBoth->SetName("hInvMassOmBoth");
  hInvMassOmBoth->SetTitle("inv Mass of Omega (bachelor consistent with both pion, kaon)");

  // TODO pT spectra inv mass afhankelijk

  // Xi invMass for different pT:
  axranges massLowpT{{corr::ptTrigg, {pTmin, 4.0}}};//, {selflagTrigg, {0.5, 2.5}}};
  axranges massMedpT{{corr::ptTrigg, {4.0, 8.0}}};//, {selflagTrigg, {0.5, 2.5}}};
  axranges massHigpT{{corr::ptTrigg, {8.0, 15.0}}};//, {selflagTrigg, {0.5, 2.5}}};
  TH1D *hInvMassXiLow = project(hXiXiOS, corr::invMassTrigg, massLowpT);
  hInvMassXiLow->SetName("hInvMassXiLow");
  hInvMassXiLow->SetTitle("inv Mass of Xi trigger (w PID response, 0 < pT < 4)");
  TH1D *hInvMassXiMed = project(hXiXiOS, corr::invMassTrigg, massMedpT);
  hInvMassXiMed->SetName("hInvMassXiMed");
  hInvMassXiMed->SetTitle("inv Mass of Xi trigger (w PID response, 4 < pT < 8)");
  TH1D *hInvMassXiHig = project(hXiXiOS, corr::invMassTrigg, massHigpT);
  hInvMassXiHig->SetName("hInvMassXiHig");
  hInvMassXiHig->SetTitle("inv Mass of Xi trigger (w PID response, 8 < pT < 15)");

  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    // old axranges, save for backup
    // axranges aXiSig{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {1, pTmax}}, 
    //                 {corr::invMassTrigg, {siglowXi[pTbin], sighighXi[pTbin]}}, {corr::invMassAssoc, {siglowXi[pTbin], sighighXi[pTbin]}}
    // };
    // axranges aXiBkg{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {1, pTmax}}, 
    //                 {corr::invMassTrigg, {bkglowXi[pTbin], bkghighXi[pTbin]}}, {corr::invMassAssoc, {bkglowXi[pTbin], bkghighXi[pTbin]}}
    // };
    // axranges aOmSig{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {1, pTmax}}, 
    //                 {corr::invMassTrigg, {siglowOm[pTbin], sighighOm[pTbin]}}, {corr::invMassAssoc, {siglowXi[pTbin], sighighXi[pTbin]}}
    // };
    // axranges aOmBkg{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {1, pTmax}}, 
    //                 {corr::invMassTrigg, {bkglowOm[pTbin], bkghighOm[pTbin]}}, {corr::invMassAssoc, {bkglowXi[pTbin], bkghighXi[pTbin]}}
    // };

    axranges aXiSig{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {1, pTmax}}, 
                    {corr::invMassTrigg, {sigXi[pTbin][0], sigXi[pTbin][1]}}, {corr::invMassAssoc, {sigXi[pTbin][0], sigXi[pTbin][1]}}
    };
    axranges aXiBkg{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {1, pTmax}}, 
                    {corr::invMassTrigg, {bkgXi[pTbin][0], bkgXi[pTbin][1]}}, {corr::invMassAssoc, {bkgXi[pTbin][0], bkgXi[pTbin][1]}}
    };
    axranges aOmSig{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {1, pTmax}}, 
                    {corr::invMassTrigg, {sigOm[pTbin][0], sigOm[pTbin][1]}}, {corr::invMassAssoc, {sigXi[pTbin][0], sigXi[pTbin][1]}}
    };
    axranges aOmBkg{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {1, pTmax}}, 
                    {corr::invMassTrigg, {bkgOm[pTbin][0], bkgOm[pTbin][1]}}, {corr::invMassAssoc, {bkgXi[pTbin][0], bkgXi[pTbin][1]}}
    };

    XiXidir->cd();
    // ad hoc Y projection
    TH1D *hXiYOS = project(hXiXiOS, corr::dY, aXiSig);
    hXiYOS->SetName("hXiYOS"+pTlabels[pTbin]);
    hXiYOS->Scale(1. / NTrigXi[pTbin]); 
    hXiYOS->SetTitle("#Xi-#Xi for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hXiYOS->SetYTitle("pairs OS");
    hXiYOS->SetStats(kFALSE);
    hXiYOS->Rebin(5);
    hXiYOS->GetYaxis()->SetRangeUser(0, 1.1*hXiYOS->GetMaximum());
    hXiYOS->SetLineWidth(3);
    // Xi-Xi sig
    TH1D *hXiOS = project(hXiXiOS, corr::dPhi, aXiSig);
    hXiOS->SetName("hXiOS"+pTlabels[pTbin]);
    hXiOS->Scale(1. / NTrigXi[pTbin]);
    hXiOS->SetTitle("#Xi-#Xi for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hXiOS->SetYTitle("pairs OS");
    hXiOS->SetStats(kFALSE);
    hXiOS->Rebin(5);
    hXiOS->GetYaxis()->SetRangeUser(0, 1.1*hXiOS->GetMaximum());
    hXiOS->SetLineWidth(3);
    TH1D *hXiSS = project(hXiXiSS, corr::dPhi, aXiSig);
    hXiSS->SetName("hXiSS"+pTlabels[pTbin]);
    hXiSS->Scale(1. / NTrigXi[pTbin]);
    hXiSS->SetTitle("#Xi-#Xi for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hXiSS->SetYTitle("pairs SS");
    hXiSS->SetStats(kFALSE);
    hXiSS->Rebin(5);
    hXiSS->GetYaxis()->SetRangeUser(0, 1.1*hXiSS->GetMaximum());
    hXiSS->SetLineWidth(3);
    TH1D *hXi = new TH1D(*hXiOS);
    hXi->Add(hXiOS, hXiSS, 1, -1);
    // if (pTbin>0) h112->Add(hXi);
    hXi->SetName("hXi"+pTlabels[pTbin]);
    hXi->SetTitle("#Xi-#Xi for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hXi->SetYTitle("pairs (OS - SS)");
    hXi->SetStats(kFALSE);
    hXi->SetLineWidth(3);
    // Xi-Xi bkg
    // TH1D *hXiOSbkg = project(hXiXiOS, corr::dPhi, aXiBkg);
    // hXiOSbkg->SetName("hXiOSbkg"+pTlabels[pTbin]);
    // TH1D *hXiSSbkg = project(hXiXiSS, corr::dPhi, aXiBkg);
    // hXiSSbkg->SetName("hXiSSbkg"+pTlabels[pTbin]);
    // TH1D *hXibkg = new TH1D(*hXiOSbkg);
    // hXibkg->Add(hXiOSbkg, hXiSSbkg, 1, -1);
    // hXibkg->SetName("hXibkg"+pTlabels[pTbin]);

    // Om-Xi sig
    OmXidir->cd();
    TH1D *hOmOS = project(hOmXiOS, corr::dPhi, aOmSig);
    hOmOS->SetName("hOmOS"+pTlabels[pTbin]);
    hOmOS->Scale(1./NTrigOm[pTbin]);
    hOmOS->SetTitle("#Omega-#Xi for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hOmOS->SetYTitle("pairs OS");
    hOmOS->SetStats(kFALSE);
    hOmOS->Rebin(5);
    hOmOS->GetYaxis()->SetRangeUser(0, 1.1*hOmOS->GetMaximum());
    hOmOS->SetLineWidth(3);
    TH1D *hOmSS = project(hOmXiSS, corr::dPhi, aOmSig);
    hOmSS->SetName("hOmSS"+pTlabels[pTbin]);
    hOmSS->Scale(1./NTrigOm[pTbin]);
    hOmSS->SetTitle("#Omega-#Xi for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hOmSS->SetYTitle("pairs SS");
    hOmSS->SetStats(kFALSE);
    hOmSS->Rebin(5);
    hOmSS->GetYaxis()->SetRangeUser(0, 1.1*hOmSS->GetMaximum());
    hOmSS->SetLineWidth(3);
    TH1D *hOm = new TH1D(*hOmOS);
    hOm->Add(hOmOS, hOmSS, 1, -1);
    hOm->SetName("hOm"+pTlabels[pTbin]);
    hOm->SetTitle("#Omega-#Xi for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hOm->SetYTitle("pairs (OS - SS)");
    hOm->SetStats(kFALSE);
    hOm->SetLineWidth(3);
    // Om-Xi bkg
    // TH1D *hOmOSbkg = project(hOmXiOS, corr::dPhi, aOmBkg);
    // hOmOSbkg->SetName("hOmOSbkg"+pTlabels[pTbin]);
    // TH1D *hOmSSbkg = project(hOmXiSS, corr::dPhi, aOmBkg);
    // hOmSSbkg->SetName("hOmSSbkg"+pTlabels[pTbin]);
    // TH1D *hOmbkg = new TH1D(*hOmOSbkg);
    // hOmbkg->Add(hOmOSbkg, hOmSSbkg, 1, -1);
    // hOmbkg->SetName("hOmbkg"+pTlabels[pTbin]);

    // Let's do Om Om
    axranges aOmOmSig{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}}, 
                    {corr::invMassTrigg, {sigOm[pTbin][0], sigOm[pTbin][1]}}, {corr::invMassAssoc, {sigOm[pTbin][0], sigOm[pTbin][1]}}
    };
    axranges aOmOmBkg{{corr::ptTrigg, {pTbins[pTbin], pTbins[pTbin+1]}}, {corr::ptAssoc, {pTmin, pTmax}}, 
                    {corr::invMassTrigg, {bkgOm[pTbin][0], bkgOm[pTbin][1]}}, {corr::invMassAssoc, {bkgOm[pTbin][0], bkgOm[pTbin][1]}}
    };
    // Om-Om sig
    OmOmdir->cd();
    TH1D *hOmOmOSsig = project(hOmOmOS, corr::dPhi, aOmOmSig);
    hOmOmOSsig->SetName("hOmOmOSsig"+pTlabels[pTbin]);
    hOmOmOSsig->Scale(1. / NTrigOm[pTbin]);
    hOmOmOSsig->SetTitle("#Omega-#Omega for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hOmOmOSsig->SetYTitle("pairs OS");
    hOmOmOSsig->SetStats(kFALSE);
    hOmOmOSsig->Rebin(5);
    hOmOmOSsig->GetYaxis()->SetRangeUser(0, 1.1*hOmOmOSsig->GetMaximum());
    hOmOmOSsig->SetLineWidth(3);
    TH1D *hOmOmSSsig = project(hOmOmSS, corr::dPhi, aOmOmSig);
    hOmOmSSsig->SetName("hOmOmSSsig"+pTlabels[pTbin]);
    hOmOmSSsig->Scale(1. / NTrigOm[pTbin]);
    hOmOmSSsig->SetTitle("#Omega-#Omega for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hOmOmSSsig->SetYTitle("pairs SS");
    hOmOmSSsig->SetStats(kFALSE);
    hOmOmSSsig->Rebin(5);
    hOmOmSSsig->GetYaxis()->SetRangeUser(0, 1.1*hOmOmSSsig->GetMaximum());
    hOmOmSSsig->SetLineWidth(3);
    TH1D *hOmOmsig = new TH1D(*hOmOmOSsig);
    hOmOmsig->Add(hOmOmOSsig, hOmOmSSsig, 1, -1);
    hOmOmsig->SetName("hOmOmsig"+pTlabels[pTbin]);
    hOmOmsig->SetTitle("#Omega-#Omega for " + pTlabels[pTbin] + " < p_{T,Trigger} < " + pTlabels[pTbin + 1]);
    hOmOmsig->SetYTitle("pairs (OS - SS)");
    hOmOmsig->SetStats(kFALSE);
    // hOmOmsig->Rebin(5);
    hOmOmsig->SetLineWidth(3);
    // Om-Om bkg
    TH1D *hOmOmOSbkg = project(hOmOmOS, corr::dPhi, aOmOmBkg);
    hOmOmOSbkg->SetName("hOmOmOSbkg"+pTlabels[pTbin]);
    TH1D *hOmOmSSbkg = project(hOmOmSS, corr::dPhi, aOmOmBkg);
    hOmOmSSbkg->SetName("hOmOmSSbkg"+pTlabels[pTbin]);
    TH1D *hOmOmbkg = new TH1D(*hOmOmOSbkg);
    hOmOmbkg->Add(hOmOmOSbkg, hOmOmSSbkg, 1, -1);
    hOmOmbkg->SetName("hOmOmbkg"+pTlabels[pTbin]);
    // sig-sig
      // do Xi, Om, both (OS, SS, subtracted)
      // then do logic with Xi + both, Om + both
    // repeat exercise for bkg-bkg
    // remember to save histos in relevant dirs
  }
 
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  return 0;
}