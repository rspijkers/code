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
const std::vector<double> pTbins = {0.6, 1.0, 2.0, 3.0, 3.5, 5.0, 10.0};
const std::vector<TString> pTlabels = {"0.6", "1.0", "2.0", "3.0", "5.0", "10.0"};
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
    ax->SetRangeUser(bounds[0], bounds[1]);
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
    ax->SetRangeUser(bounds[0], bounds[1]);
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

// define function that projects all the QA histograms?
void doQAprojections(TFile* infile) { // todo fix infile
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

void doXiInvMassFits(){
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
  axranges aMEXiXiOS{{corr::ptTrigg, {1.0, pTmax}}, {corr::ptAssoc, {1.0, pTmax}}, 
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
  hXiXiOSdPhi->Rebin(5);
  hXiXiOSdPhi->Scale(1. / hXiXiOSdPhi->GetMaximum());

  TH1D *hXiXiSSdPhi = project(hXiXiSS, corr::dPhi, aMEXiXiOS);
  hXiXiSSdPhi->SetName("hMEXiSSdPhi");
  hXiXiSSdPhi->SetTitle("ME #Xi-#Xi SS (pT integrated)");
  hXiXiSSdPhi->SetYTitle("pairs SS");
  hXiXiSSdPhi->SetLineWidth(3);
  hXiXiSSdPhi->Rebin(5);
  hXiXiSSdPhi->Scale(1. / hXiXiSSdPhi->GetMaximum());

  return {hXiXiOSdPhi, hXiXiSSdPhi, hXiXiOSdY, hXiXiSSdY};
}

int postprocessingResults(TString trainnr, TString filename = "AnalysisResults.root") {
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
  TH1F *hPhi = inputDir->Get<TH1F>("hPhi");
  hPhi->SetDirectory(outputFile);
  hPhi->Draw();

  doQAprojections(inputFile);
  // test ME
  std::vector<TH1D*> hME_vector;
  hME_vector = doMixedEvents(trainnr);

  // test inv mass fitting
  doXiInvMassFits();
  doOmInvMassFits();
  // todo: assert some properties of the inv mass region boundaries to ensure everything went ok.

  // now that we have the inv mass regions, let's project and integrate the eff corrected inv mass plots to determine the number of triggers
  TDirectory* XiEffdir = outputFile->mkdir("XiEffdir");
  TDirectory* OmEffdir = outputFile->mkdir("OmEffdir");
  THnSparse *hEffCorrXiMass, *hEffCorrOmegaMass;
  inputDir->GetObject("hMassXiEffCorrected", hEffCorrXiMass);
  inputDir->GetObject("hMassOmegaEffCorrected", hEffCorrOmegaMass);
  // put the projections in a vector so we can access them later
  // array of Ntriggers for Xi, Omega
  double NTrigXi[maxPtBins - 1], NTrigOm[maxPtBins - 1];
  std::vector<TH1D*> v_effCorrXiMass, v_effCorrOmMass; // obsolete?
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
    // v_effCorrXiMass.push_back(hXi);
    TH1D *hOm = project(hEffCorrOmegaMass, mass::invMass, aOm);
    hOm->SetName("hEffOmMass_" + pTlabels[pTbin + 1]);
    hOm->SetDirectory(OmEffdir);
    double omInt = hOm->Integral();
    if(std::isinf(omInt)) {
      omInt = 1.; // set to 1 in case of infinity, warn the user
      cout << "WARNING: Ntrig Om in pT bin " << pTlabels[pTbin] << " - " << pTlabels[pTbin + 1] << " is inf - Ntrig set to 1." << endl;
    }
    NTrigOm[pTbin] = omInt;
    // v_effCorrOmMass.push_back(*hOm);
  }
  // // todo remove these couts
  // for (auto i : NTrigXi){
  //   cout << "test Xi " << i << endl;
  // }
  // for (auto i : NTrigOm){
  //   cout << "test Om " << i << endl;
  // }
  
  /*
  // TODO MAKE DEDICATED FUNCTION FOR INV MASS FITTING, CALL 4 TIMES (xi+/-, omega+/-)
  // Xi mass
  std::vector<double> siglowXi, sighighXi;
  std::vector<double> bkglowXi, bkghighXi;
  std::vector<double> NTrigXiSig, NTrigXiBkg;
  TF1 *f1 = new TF1("f1", "pol2(0) + gaus(3) + gaus(6)", 1.29, 1.42); 
  f1->SetParameters(0, 0, 10, 0, 1.321, 0.005, 0, 1.321, 0.001); 
  // f1->SetParLimits(0, -1000, 0);
  f1->SetParLimits(4, 1.31, 1.33);
  f1->SetParLimits(7, 1.31, 1.33);
  f1->SetParLimits(5, 0, 0.01);
  f1->SetParLimits(8, 0, 0.01);

  TF1 *fBKGXi = new TF1("fBKGXi", pol2bkgXi, 1.29, 1.42, 3);
  
  // do inv mass fit in pT bins
  TH2F *hMassXiMinus = inputDir->Get<TH2F>("hMassXiMinus");
  hMassXiMinus->SetDirectory(outputFile);
  cout << "Start Xi mass fitting..." << endl;
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hMassXiMinus->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin+1]);
    TH1D *h = hMassXiMinus->ProjectionX("hMassXiMinus_"+ pTlabels[pTbin + 1]);
    h->GetXaxis()->SetRangeUser(1.2, 1.5);
    h->SetTitle("#Xi^{-} inv. mass for " + pTlabels[pTbin] + " GeV < p_{T} < " + pTlabels[pTbin+1] + " GeV");
    // Double_t nentries = h->Integral();
    // h->Scale(1. / nentries, "width");
    f1->SetParameter(3, .4* h->GetMaximum());
    f1->SetParameter(6, .4* h->GetMaximum());
    h->Fit(fBKGXi, "SLBQRO");
    f1->SetParameter(0, fBKGXi->GetParameter(0));
    f1->SetParameter(1, fBKGXi->GetParameter(1));
    f1->SetParameter(2, fBKGXi->GetParameter(2));
    TFitResultPtr r1 = h->Fit("f1", "SLBQR", "", 1.29, 1.42);
    if (r1->Chi2() > 1) cout << "Warning: Fit Chi2(" << r1->Chi2() << ") > 1 in pT bin " << pTbin << endl;
    // calculate #triggers in region:
    double mu = (f1->GetParameter(4) + f1->GetParameter(7)) / 2.;
    double sigma = (f1->GetParameter(5) + f1->GetParameter(8)) / 2.;
    NTrigXiSig.push_back(h->Integral(h->FindBin(mu-3*sigma), h->FindBin(mu+3*sigma)));
    NTrigXiBkg.push_back(h->Integral(h->FindBin(mu+4*sigma), h->FindBin(mu+10*sigma)));
    // cout << v_effCorrXiMass[pTbin]->Integral(v_effCorrXiMass[pTbin]->FindBin(mu-3*sigma), v_effCorrXiMass[pTbin]->FindBin(mu+3*sigma)) << endl;
    // NTrigXiSig.push_back(v_effCorrXiMass[pTbin]->Integral(h->FindBin(mu-3*sigma), h->FindBin(mu+3*sigma)));
    // NTrigXiBkg.push_back(v_effCorrXiMass[pTbin]->Integral(h->FindBin(mu+4*sigma), h->FindBin(mu+10*sigma)));
    h->GetXaxis()->SetRangeUser(1.28, 1.38);
    h->SetStats(kFALSE);
    // h->SetLineWidth(3);
    // save upper/lower limits of sig/bkg regions
    siglowXi.push_back(mu - 3*sigma);
    sighighXi.push_back(mu + 3*sigma);
    bkglowXi.push_back(mu + 4*sigma);
    bkghighXi.push_back(mu + 10*sigma);
  }

  // AD HOC XI PLUS FIXME FIXME FIXME
  TH2F *hMassXiPlus = inputDir->Get<TH2F>("hMassXiPlus");
  hMassXiPlus->SetDirectory(outputFile);
  cout << "Start Xi mass fitting..." << endl;
  for (int pTbin = 0; pTbin < maxPtBins - 1; pTbin++){
    hMassXiPlus->GetYaxis()->SetRangeUser(pTbins[pTbin], pTbins[pTbin+1]);
    TH1D *h = hMassXiPlus->ProjectionX("hMassXiPlus_"+ pTlabels[pTbin + 1]);
    h->GetXaxis()->SetRangeUser(1.2, 1.5);
    h->SetTitle("#Xi^{-} inv. mass for " + pTlabels[pTbin] + " GeV < p_{T} < " + pTlabels[pTbin+1] + " GeV");
    // Double_t nentries = h->Integral();
    // h->Scale(1. / nentries, "width");
    f1->SetParameter(3, .4* h->GetMaximum());
    f1->SetParameter(6, .4* h->GetMaximum());
    h->Fit(fBKGXi, "SLBQRO");
    f1->SetParameter(0, fBKGXi->GetParameter(0));
    f1->SetParameter(1, fBKGXi->GetParameter(1));
    f1->SetParameter(2, fBKGXi->GetParameter(2));
    TFitResultPtr r1 = h->Fit("f1", "SLBQR", "", 1.29, 1.42);
    if (r1->Chi2() > 1) cout << "Warning: Fit Chi2(" << r1->Chi2() << ") > 1 in pT bin " << pTbin << endl;
    // calculate #triggers in region:
    double mu = (f1->GetParameter(4) + f1->GetParameter(7)) / 2.;
    double sigma = (f1->GetParameter(5) + f1->GetParameter(8)) / 2.;
    // NTrigXiSig.push_back(h->Integral(h->FindBin(mu-3*sigma), h->FindBin(mu+3*sigma)));
    // NTrigXiBkg.push_back(h->Integral(h->FindBin(mu+4*sigma), h->FindBin(mu+10*sigma)));
    // cout << v_effCorrXiMass[pTbin]->Integral(v_effCorrXiMass[pTbin]->FindBin(mu-3*sigma), v_effCorrXiMass[pTbin]->FindBin(mu+3*sigma)) << endl;
    // NTrigXiSig.push_back(v_effCorrXiMass[pTbin]->Integral(h->FindBin(mu-3*sigma), h->FindBin(mu+3*sigma)));
    // NTrigXiBkg.push_back(v_effCorrXiMass[pTbin]->Integral(h->FindBin(mu+4*sigma), h->FindBin(mu+10*sigma)));
    h->GetXaxis()->SetRangeUser(1.28, 1.38);
    h->SetStats(kFALSE);
    // h->SetLineWidth(3);
    // save upper/lower limits of sig/bkg regions
    // siglowXi.push_back(mu - 3*sigma);
    // sighighXi.push_back(mu + 3*sigma);
    // bkglowXi.push_back(mu + 4*sigma);
    // bkghighXi.push_back(mu + 10*sigma);
  }

  // Omega mass
  std::vector<double> siglowOm, sighighOm;
  std::vector<double> bkglowOm, bkghighOm;
  std::vector<double> NTrigOmSig, NTrigOmBkg;
  TF1 *f2 = new TF1("f2", "pol2(0) + gaus(3) + gaus(6)", 1.64, 1.74); 
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
    // Double_t nentries = h->Integral();
    // h->Scale(1. / nentries, "width");
    f2->SetParameter(3, 0.3* h->GetMaximum());
    f2->SetParameter(6, 0.4* h->GetMaximum());
    h->Fit(fBKG, "SLBQRO");
    f2->SetParameter(0, fBKG->GetParameter(0));
    f2->SetParameter(1, fBKG->GetParameter(1));
    f2->SetParameter(2, fBKG->GetParameter(2));

    TFitResultPtr r2 = h->Fit("f2", "SLQBR", "", 1.64, 1.74);
    if (r2->Chi2() > 1) cout << "Warning: Fit Chi2(" << r2->Chi2() << ") > 1 in pT bin " << pTbin << endl;

    // calculate #triggers in region:
    double mu = (f2->GetParameter(4) + f2->GetParameter(7)) / 2.;
    double sigma = (f2->GetParameter(5) + f2->GetParameter(8)) / 2.;
    NTrigOmSig.push_back(h->Integral(h->FindBin(mu-3*sigma), h->FindBin(mu+3*sigma)));
    NTrigOmBkg.push_back(h->Integral(h->FindBin(mu+4*sigma), h->FindBin(mu+10*sigma)));

    h->GetXaxis()->SetRangeUser(1.63, 1.72);
    h->SetStats(kFALSE);
    h->SetLineWidth(3);

    // save upper/lower limits of sig/bkg regions
    siglowOm.push_back(mu - 3*sigma);
    sighighOm.push_back(mu + 3*sigma);
    bkglowOm.push_back(mu + 4*sigma);
    bkghighOm.push_back(mu + 10*sigma);
    if(mu + 10*sigma > 1.76) cout << "Warning: Bkg region exceeds limit of 1.76 GeV in pT bin " << pTbin << endl;
  }
  for(int i = 0; i < maxPtBins - 1; i++){
    cout << NTrigXiSig[i] << " <- Xi, Om -> " << NTrigOmSig[i] << endl;
  }
  */

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
  axranges aPtIntegrated{{corr::ptTrigg, {1., pTmax}}, {corr::ptAssoc, {1., pTmax}},
                 {corr::invMassTrigg, {1.31, 1.33}}, {corr::invMassAssoc, {1.31, 1.33}}
  };
  TH2D *h2DtestOS = project2D(hXiXiOS, corr::dPhi, corr::dY, aPtIntegrated);
  h2DtestOS->SetName("h2DtestOS");
  TH2D *h2DtestSS = project2D(hXiXiSS, corr::dPhi, corr::dY, aPtIntegrated);
  h2DtestSS->SetName("h2DtestSS");
  TH2D *h2Dtest = new TH2D(*h2DtestOS);
  h2Dtest->Add(h2DtestOS, h2DtestSS, 1, -1);
  h2Dtest->SetName("h2Dtest");
  h2Dtest->ProjectionX("hDYtest");
  h2Dtest->ProjectionY("hDPhitest");

  // quick pT integrated ME correction applied
  TH1D* hDYOScorrected = project(hXiXiOS, corr::dY, aPtIntegrated);
  hDYOScorrected->SetName("hDYOScorrected");
  hDYOScorrected->Divide(hME_vector[2]);
  TH1D* hDYSScorrected = project(hXiXiSS, corr::dY, aPtIntegrated);
  hDYSScorrected->SetName("hDYSScorrected");
  hDYSScorrected->Divide(hME_vector[3]);
  TH1D *hDYcorrected = new TH1D(*hDYOScorrected);
  hDYcorrected->Add(hDYOScorrected, hDYSScorrected, 1, -1);
  hDYcorrected->SetName("hDYcorrected");
  
  TH1D* hDPhiOScorrected = project(hXiXiOS, corr::dPhi, aPtIntegrated);
  hDPhiOScorrected->SetName("hDPhiOScorrected");
  hDPhiOScorrected->Divide(hME_vector[0]);
  TH1D* hDPhiSScorrected = project(hXiXiSS, corr::dPhi, aPtIntegrated);
  hDPhiSScorrected->SetName("hDPhiSScorrected");
  hDPhiSScorrected->Divide(hME_vector[1]);
  TH1D *hDPhicorrected = new TH1D(*hDPhiOScorrected);
  hDPhicorrected->Add(hDPhiOScorrected, hDPhiSScorrected, 1, -1);
  hDPhicorrected->SetName("hDPhicorrected");

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

  TH1D *hTest = project(hXiXiOS, corr::dPhi, massXi);
  TH1D *h112 = new TH1D(*hTest);
  h112->Reset();
  h112->Rebin(5);
  h112->SetName("h112");
  
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
    if (pTbin>0) h112->Add(hXi);
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
  /* 
  TODO: what is this again and what do we do with it?
  double totaltriggers = 0;
  for(int i = 1; i < maxPtBins -1; i++) totaltriggers += NTrigXiSig[i];
  cout << totaltriggers << endl;
  h112->Scale(1. / totaltriggers);
  h112->SetTitle("#Xi-#Xi correlations for 1 < p_{T,Trig} < 15 & 0.15 < p_{T,Assoc} < p_{T,Trig}");
  h112->SetYTitle("OS - SS");
  h112->SetStats(kFALSE);
  h112->Rebin(5);
  h112->GetYaxis()->SetRangeUser(0, 0.0002);
  h112->SetLineWidth(3);
  */
 
  outputFile->cd();
  outputFile->Write();
  outputFile->Close();
  return 0;
}