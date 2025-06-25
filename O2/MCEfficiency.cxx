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

int MCEfficiency(TString trainnr, TString id = "", TString filename = "AnalysisResults.root") {
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
  
  // //// 1D Xi Minus
  // TH1F *hXiGenF = gen->Get<TH1F>("hPtXiMinus");
  // TH1D *hXiGen = new TH1D();
  // hXiGenF->Copy(*hXiGen);
  // hXiGen->SetDirectory(outfile);
  // cout << hXiGen->GetXaxis()->GetNbins() << endl;
  // TH1D *hXiGenRebin = dynamic_cast<TH1D*>(hXiGen->Rebin(maxPtBins-1, "hXiGenRebin", pTbins));

  // TH2F *hXiRec2D = rec->Get<TH2F>("h2dMassXiMinus");
  // hXiRec2D->SetDirectory(outfile);
  // TH1D *hXiRec = hXiRec2D->ProjectionX();
  // cout << hXiRec->GetXaxis()->GetNbins() << endl;
  // hXiRec->SetName("hXiRec");
  // TH1D *hXiRecRebin = dynamic_cast<TH1D*>(hXiRec->Rebin(maxPtBins-1, "hXiRecRebin", pTbins));

  // // quick check to count generated Ximin:
  // cout << "# Xi gen. = " << hXiGen->Integral(hXiGen->FindBin(1.0), hXiGen->FindBin(10.)) << endl;

  // // TCanvas *c = new TCanvas();
  // TH1D *hXiMinEff = new TH1D(*hXiGenRebin);
  // hXiMinEff->Divide(hXiRecRebin, hXiGenRebin);
  // hXiMinEff->SetName("hXiMinEff1D");
  // hXiMinEff->SetStats(kFALSE);
  // hXiMinEff->SetTitle("#Xi^{-} efficiency");
  // hXiMinEff->GetYaxis()->SetTitle("#varepsilon");
  // // hXiMinEff->Draw();

  // //// 2D Xi Minus
  // // get the relevant histo's from file
  // TH2D *hXiMin2DGen = (TH2D*) gen->Get<TH2F>("h2DXiMinus");
  // TH3D *hXiMin3DRec = (TH3D*) rec->Get<TH3F>("hPtYMassXiMinus");
  // TH2D *hXiMin2DRec = (TH2D*) hXiMin3DRec->Project3D("yx");

  // // gen
  // TH2D* h2DXiMinGen = RebinX2D(hXiMin2DGen, maxPtBins, pTbins, "h2DXiMinGen");
  // h2DXiMinGen = (TH2D*) h2DXiMinGen->RebinY(10);
  // // projections for 1D eff sanity checks
  // TH1D* hPtXiMinGen = (TH1D*) h2DXiMinGen->ProjectionX("hPtXiMinGen");
  // TH1D* hYXiMinGen = (TH1D*) h2DXiMinGen->ProjectionY("hYXiMinGen");

  // // rec
  // TH2D* h2DXiMinRec = RebinX2D(hXiMin2DRec, maxPtBins, pTbins, "h2DXiMinRec");
  // h2DXiMinRec = (TH2D*) h2DXiMinRec->RebinY(10);
  // // projections for 1D eff sanity checks
  // TH1D* hPtXiMinRec = (TH1D*) h2DXiMinRec->ProjectionX("hPtXiMinRec");
  // TH1D* hYXiMinRec = (TH1D*) h2DXiMinRec->ProjectionY("hYXiMinRec");

  // // eff
  // TH2D *h2DXiMinEff = new TH2D(*h2DXiMinGen);
  // h2DXiMinEff->Divide(h2DXiMinRec, h2DXiMinGen);
  // h2DXiMinEff->SetName("hXiMinEff");
  // h2DXiMinEff->SetTitle("2D efficiency of #Xi^{-}");
  // // cout << "test: 5.830e-6 = " << h2DXiMinEff->GetBinContent(h2DXiMinEff->FindFixBin(0.3, 0.05)) << endl; 

  // // 1D eff sanity checks
  // TH1D *hPtXiMinEff = new TH1D(*hPtXiMinGen);
  // hPtXiMinEff->Divide(hPtXiMinRec, hPtXiMinGen);
  // hPtXiMinEff->SetName("hPtXiMinEff");
  // hPtXiMinEff->SetTitle("Efficiency of #Xi^{-} vs. p_T");
  // TH1D *hYXiMinEff = new TH1D(*hYXiMinGen);
  // hYXiMinEff->Divide(hYXiMinRec, hYXiMinGen);
  // hYXiMinEff->SetName("hYXiMinEff");
  // hYXiMinEff->SetTitle("Efficiency of #Xi^{-} vs. y");

  //// 2D Xi Minus
  // get the relevant histo's from file
  TH2D *hXiMin2DGen = (TH2D*) infile->Get<TH2F>("cascade-selector/gen/hXiMinus");
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
  // cout << "test: 5.830e-6 = " << h2DXiMinEff->GetBinContent(h2DXiMinEff->FindFixBin(0.3, 0.05)) << endl; 

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
  TH2D *hXiPlus2DGen = (TH2D*)infile->Get<TH2F>("cascade-selector/gen/hXiPlus");
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
  TH2D *hOmegaMin2DGen = (TH2D*)infile->Get<TH2F>("cascade-selector/gen/hOmegaMinus");
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
  TH2D *hOmegaPlus2DGen = (TH2D*)infile->Get<TH2F>("cascade-selector/gen/hOmegaPlus");
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

  // //// Omega Minus
  // TH1F *hOmegaGenF = gen->Get<TH1F>("hPtOmegaMinus");
  // TH1D *hOmegaGen = new TH1D();
  // hOmegaGenF->Copy(*hOmegaGen);
  // hOmegaGen->SetDirectory(outfile);
  // cout << hOmegaGen->GetXaxis()->GetNbins() << endl;
  // TH1D *hOmegaGenRebin = dynamic_cast<TH1D*>(hOmegaGen->Rebin(maxPtBins-1, "hOmegaGenRebin", pTbins));

  // TH2F *hOmegaRec2D = rec->Get<TH2F>("h2dMassOmegaMinus");
  // hOmegaRec2D->SetDirectory(outfile);
  // TH1D *hOmegaRec = hOmegaRec2D->ProjectionX();
  // cout << hOmegaRec->GetXaxis()->GetNbins() << endl;
  // hOmegaRec->SetName("hOmegaRec");
  // TH1D *hOmegaRecRebin = dynamic_cast<TH1D*>(hOmegaRec->Rebin(maxPtBins-1, "hOmegaRecRebin", pTbins));

  // TH1D *hOmegaMinEff = new TH1D(*hOmegaGenRebin);
  // hOmegaMinEff->Divide(hOmegaRecRebin, hOmegaGenRebin);
  // hOmegaMinEff->SetName("hOmegaMinEff");
  // hOmegaMinEff->SetStats(kFALSE);
  // hOmegaMinEff->SetTitle("#Omega^{-} efficiency");
  // hOmegaMinEff->GetYaxis()->SetTitle("#varepsilon");
  // hOmegaMinEff->Draw();
  // c->Print("figures/eff/EffOmegaMin.pdf");
  // c->Clear();

  // //// Xi Plus
  // TH1F *hXiPlusGenF = gen->Get<TH1F>("hPtXiPlus");
  // TH1D *hXiPlusGen = new TH1D();
  // hXiPlusGenF->Copy(*hXiPlusGen);
  // hXiPlusGen->SetDirectory(outfile);
  // cout << hXiPlusGen->GetXaxis()->GetNbins() << endl;
  // TH1D *hXiPlusGenRebin = dynamic_cast<TH1D*>(hXiPlusGen->Rebin(maxPtBins-1, "hXiPlusGenRebin", pTbins));

  // TH2F *hXiPlusRec2D = rec->Get<TH2F>("h2dMassXiPlus");
  // hXiPlusRec2D->SetDirectory(outfile);
  // TH1D *hXiPlusRec = hXiPlusRec2D->ProjectionX();
  // cout << hXiPlusRec->GetXaxis()->GetNbins() << endl;
  // hXiPlusRec->SetName("hXiPlusRec");
  // TH1D *hXiPlusRecRebin = dynamic_cast<TH1D*>(hXiPlusRec->Rebin(maxPtBins-1, "hXiPlusRecRebin", pTbins));

  // TH1D *hXiPlusEff = new TH1D(*hXiPlusGenRebin);
  // hXiPlusEff->Divide(hXiPlusRecRebin, hXiPlusGenRebin);
  // hXiPlusEff->SetName("hXiPlusEff");
  // hXiPlusEff->SetStats(kFALSE);
  // hXiPlusEff->SetTitle("#Xi^{+} efficiency");
  // hXiPlusEff->GetYaxis()->SetTitle("#varepsilon");
  // hXiPlusEff->Draw();
  // c->Print("figures/eff/EffXiPlus.pdf");
  // c->Clear();

  // //// Omega Plus
  // TH1F *hOmegaPlusGenF = gen->Get<TH1F>("hPtOmegaPlus");
  // TH1D *hOmegaPlusGen = new TH1D();
  // hOmegaPlusGenF->Copy(*hOmegaPlusGen);
  // hOmegaPlusGen->SetDirectory(outfile);
  // cout << hOmegaPlusGen->GetXaxis()->GetNbins() << endl;
  // TH1D *hOmegaPlusGenRebin = dynamic_cast<TH1D*>(hOmegaPlusGen->Rebin(maxPtBins-1, "hOmegaPlusGenRebin", pTbins));

  // TH2F *hOmegaPlusRec2D = rec->Get<TH2F>("h2dMassOmegaPlus");
  // hOmegaPlusRec2D->SetDirectory(outfile);
  // TH1D *hOmegaPlusRec = hOmegaPlusRec2D->ProjectionX();
  // cout << hOmegaPlusRec->GetXaxis()->GetNbins() << endl;
  // hOmegaPlusRec->SetName("hOmegaPlusRec");
  // TH1D *hOmegaPlusRecRebin = dynamic_cast<TH1D*>(hOmegaPlusRec->Rebin(maxPtBins-1, "hOmegaPlusRecRebin", pTbins));

  // TH1D *hOmegaPlusEff = new TH1D(*hOmegaPlusGenRebin);
  // hOmegaPlusEff->Divide(hOmegaPlusRecRebin, hOmegaPlusGenRebin);
  // hOmegaPlusEff->SetName("hOmegaPlusEff");
  // hOmegaPlusEff->SetStats(kFALSE);
  // hOmegaPlusEff->SetTitle("#Omega^{+} efficiency");
  // hOmegaPlusEff->GetYaxis()->SetTitle("#varepsilon");
  // hOmegaPlusEff->Draw();
  // c->Print("figures/eff/EffOmegaPlus.pdf");
  // c->Clear();

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