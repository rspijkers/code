// This script makes the relevant plots/projections from a given AnalysisResults.root.
// Run me like `root 'MCEfficiency.cxx("trainnr", "outfilename")' -q`, where the second 
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
#include "TStyle.h"
#include "TCanvas.h"

using std::cout; using std::endl;

// temp cut variations
// default cascradius 1.01, casccpa 0.9947, dcacascdau 0.25
// id24075: def 
// id24076: cascradius 0.5
// id24077: casccpa 0.99
// id24078: dcacascdau 0.8
// id24079: all

TH2D* RebinX2D(TH2D* h, int nbins, double* bins, TString name = "hrebin") {

  TAxis *xaxis = h->GetXaxis();
  TAxis *yaxis = h->GetYaxis();
  TH2D *hnew = new TH2D(name, name, nbins-1, bins, h->GetNbinsY(), yaxis->GetXmin(), yaxis->GetXmax());

  for (int j=1; j<=yaxis->GetNbins(); j++) {
    for (int i=1; i<=xaxis->GetNbins(); i++) {
      hnew->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h->GetBinContent(i,j));
    }
  }
  // this does not account for proper error propagation, but it should just be sqrt(N) so set the errors manually:
  for (int xbin = 0; xbin <= hnew->GetNbinsX(); xbin++){
    for (int ybin = 0; ybin <= hnew->GetNbinsY(); ybin++){
      int bin = hnew->GetBin(xbin, ybin);
      hnew->SetBinError(bin, std::sqrt(hnew->GetBinContent(bin)));
    }
  }
  return hnew;
}

int MCEfficiency(TString trainnr, bool makePDF = false, TString option = "correlations", TString id = "", TString filename = "AnalysisResults.root") {
  TH1::SetDefaultSumw2(); // Make sure we propagate the errors
  gStyle->SetHistLineWidth(3);
  gROOT->ForceStyle();

  TFile *infile;
  if(trainnr == "test"){ // if not a trainnr but just a test, look for AnalysisResults.root in current dir
    infile = new TFile(filename, "READ");
  } else {
    infile = new TFile("results/" + trainnr + "/" + filename, "READ");
  }
  TDirectory *rec, *gen;
  infile->GetObject("cascade-analysis-m-c" + id, rec);
  infile->GetObject("cascade-generated" + id, gen);
  // ROOT is being a bitch about TBrowser in batch mode, so make an outfile:
  TFile *outfile = new TFile("plots/MC" + trainnr + id + ".root", "RECREATE");

  // ad hoc pT bins:
  // const int maxPtBins = 4;
  // double pTbins[maxPtBins] = {pTmin, 4., 8., pTmax}; // edges of the different (trigger) pT ranges
  // const TString pTlabels[maxPtBins] = {"Low", "Med", "High", "ptmax"};
  const int maxPtBins = 23;
  double pTbins[maxPtBins] = {0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0}; // edges of the different (trigger) pT ranges
  // const TString pTlabels[maxPtBins] = {"pTmin", "1", "2", "3", "4", "5", "6", "8", "pTmax"};

  // require reconstructed counterpart for each generated collision in case of efficiency for correlations,
  // otherwise simply do rec/gen for spectra and calculate event factor afterwards.
  TString gendir;
  if (option == "correlations"){
    gendir = "genwithrec";
  } else if (option == "spectra"){
    gendir = "gen";
  } else {
    cout << "Invalid option! expects either 'correlations' or 'spectra', not '" << option << "'" << endl;
    return 1;
  }

  //// 2D Xi Minus
  // get the relevant histo's from file
  TH2D *hXiMin2DGen = (TH2D*) infile->Get<TH2F>("cascade-selector/"+gendir+"/hXiMinus");
  TH2D *hXiMin2DRec = (TH2D*) infile->Get<TH2F>("cascade-selector/truerec/hXiMinus");

  // gen
  TH2D* h2DXiMinGen = RebinX2D(hXiMin2DGen, maxPtBins, pTbins, "h2DXiMinGen");
  h2DXiMinGen = (TH2D*) h2DXiMinGen->RebinY(5);
  // projections for 1D eff sanity checks
  TH1D* hPtXiMinGen = (TH1D*) h2DXiMinGen->ProjectionX("hPtXiMinGen");
  TH1D* hYXiMinGen = (TH1D*) h2DXiMinGen->ProjectionY("hYXiMinGen");

  // rec
  TH2D* h2DXiMinRec = RebinX2D(hXiMin2DRec, maxPtBins, pTbins, "h2DXiMinRec");
  h2DXiMinRec = (TH2D*) h2DXiMinRec->RebinY(5);
  // projections for 1D eff sanity checks
  TH1D* hPtXiMinRec = (TH1D*) h2DXiMinRec->ProjectionX("hPtXiMinRec");
  TH1D* hYXiMinRec = (TH1D*) h2DXiMinRec->ProjectionY("hYXiMinRec");

  // eff
  TH2D *h2DXiMinEff = new TH2D(*h2DXiMinGen);
  h2DXiMinEff->Divide(h2DXiMinRec, h2DXiMinGen);
  h2DXiMinEff->SetName("hXiMinEff");
  h2DXiMinEff->SetTitle("2D efficiency of #Xi^{-}");

  // 1D eff sanity checks
  TH1D *hPtXiMinEff = new TH1D(*hPtXiMinGen);
  hPtXiMinEff->Divide(hPtXiMinRec, hPtXiMinGen);
  hPtXiMinEff->SetName("hPtXiMinEff");
  hPtXiMinEff->SetTitle("Efficiency of #Xi^{-} vs. p_T");
  TH1D *hYXiMinEff = new TH1D(*hYXiMinGen);
  hYXiMinEff->Divide(hYXiMinRec, hYXiMinGen);
  hYXiMinEff->SetName("hYXiMinEff");
  hYXiMinEff->SetTitle("Efficiency of #Xi^{-} vs. y");

  //// 2D Xi Plus
  // get the relevant histo's from file
  TH2D *hXiPlus2DGen = (TH2D*)infile->Get<TH2F>("cascade-selector/"+gendir+"/hXiPlus");
  TH2D *hXiPlus2DRec = (TH2D*)infile->Get<TH2F>("cascade-selector/truerec/hXiPlus");

  // gen
  TH2D* h2DXiPlusGen = RebinX2D(hXiPlus2DGen, maxPtBins, pTbins, "h2DXiPlusGen");
  h2DXiPlusGen = (TH2D*) h2DXiPlusGen->RebinY(5);
  // rec
  TH2D* h2DXiPlusRec = RebinX2D(hXiPlus2DRec, maxPtBins, pTbins, "h2DXiPlusRec");
  h2DXiPlusRec = (TH2D*) h2DXiPlusRec->RebinY(5);
  // eff
  TH2D *h2DXiPlusEff = new TH2D(*h2DXiPlusGen);
  h2DXiPlusEff->Divide(h2DXiPlusRec, h2DXiPlusGen);
  h2DXiPlusEff->SetName("hXiPlusEff");
  h2DXiPlusEff->SetTitle("2D efficiency of #Xi^{+}");

  //// 2D Omega Minus
  // get the relevant histo's from file
  TH2D *hOmegaMin2DGen = (TH2D*)infile->Get<TH2F>("cascade-selector/"+gendir+"/hOmegaMinus");
  TH2D *hOmegaMin2DRec = (TH2D*)infile->Get<TH2F>("cascade-selector/truerec/hOmegaMinus");

  // gen
  TH2D* h2DOmegaMinGen = RebinX2D(hOmegaMin2DGen, maxPtBins, pTbins, "h2DOmegaMinGen");
  h2DOmegaMinGen = (TH2D*) h2DOmegaMinGen->RebinY(5);
  // rec
  TH2D* h2DOmegaMinRec = RebinX2D(hOmegaMin2DRec, maxPtBins, pTbins, "h2DOmegaMinRec");
  h2DOmegaMinRec = (TH2D*) h2DOmegaMinRec->RebinY(5);
  // eff
  TH2D *h2DOmegaMinEff = new TH2D(*h2DOmegaMinGen);
  h2DOmegaMinEff->Divide(h2DOmegaMinRec, h2DOmegaMinGen);
  h2DOmegaMinEff->SetName("hOmegaMinEff");
  h2DOmegaMinEff->SetTitle("2D efficiency of #Omega^{-}");

  //// 2D Omega Plus
  // get the relevant histo's from file
  TH2D *hOmegaPlus2DGen = (TH2D*)infile->Get<TH2F>("cascade-selector/"+gendir+"/hOmegaPlus");
  TH2D *hOmegaPlus2DRec = (TH2D*)infile->Get<TH2F>("cascade-selector/truerec/hOmegaPlus");

  // gen
  TH2D* h2DOmegaPlusGen = RebinX2D(hOmegaPlus2DGen, maxPtBins, pTbins, "h2DOmegaPlusGen");
  h2DOmegaPlusGen = (TH2D*) h2DOmegaPlusGen->RebinY(5);
  // rec
  TH2D* h2DOmegaPlusRec = RebinX2D(hOmegaPlus2DRec, maxPtBins, pTbins, "h2DOmegaPlusRec");
  h2DOmegaPlusRec = (TH2D*) h2DOmegaPlusRec->RebinY(5);
  // eff
  TH2D *h2DOmegaPlusEff = new TH2D(*h2DOmegaPlusGen);
  h2DOmegaPlusEff->Divide(h2DOmegaPlusRec, h2DOmegaPlusGen);
  h2DOmegaPlusEff->SetName("hOmegaPlusEff");
  h2DOmegaPlusEff->SetTitle("2D efficiency of #Omega^{+}");

  // write canvasses to outfile and make pdf's
  TDirectory* figures = outfile->mkdir("figures");
  figures->cd();
  TCanvas *c = new TCanvas("c");
  c->cd();
  for (auto h : {h2DXiMinEff, h2DXiPlusEff, h2DOmegaMinEff, h2DOmegaPlusEff}){
    h->SetStats(kFALSE);
    h->SetXTitle("#it{p}_{T}");
    h->SetYTitle("y");
    h->Draw("colz");
    TString name = h->GetName();
    c->Write(name);
    if(makePDF) c->Print("figures/eff/2D/" + name + ".pdf");
    c->Clear();
  }

  // write a list with the efficiency histo's to a file for ccdb
  TFile* ccdbFile = new TFile("CascadeEfficiencies.root", "RECREATE");
  ccdbFile->cd();
  TList* effList = new TList();
  effList->SetName("ccdb_object");
  effList->Add(h2DXiMinEff);
  effList->Add(h2DXiPlusEff);
  effList->Add(h2DOmegaMinEff);
  effList->Add(h2DOmegaPlusEff);
  effList->Write("ccdb_object", 1);
  ccdbFile->Write();
  ccdbFile->Close();

  outfile->Write();
  outfile->Close();

  return 0;
}