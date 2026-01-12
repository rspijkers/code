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

// using json = nlohmann::json;
using std::cout; using std::endl;
using axranges = std::map<int, std::vector<double>>;


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
    if(bounds[1] == ax->GetBinLowEdge(ax->FindFixBin(bounds[1]))){
      // if the upper boundary is at the lower edge of a bin, we actually want to integrate up to the bin before
      ax->SetRange(ax->FindFixBin(bounds[0]), ax->FindFixBin(bounds[1]) - 1);
    } else {
      ax->SetRangeUser(bounds[0], bounds[1]);
    }
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
    if(bounds[1] == ax->GetBinLowEdge(ax->FindFixBin(bounds[1]))){
      // if the upper boundary is at the lower edge of a bin, we actually want to integrate up to the bin before
      ax->SetRange(ax->FindFixBin(bounds[0]), ax->FindFixBin(bounds[1]) - 1);
    } else {
      ax->SetRangeUser(bounds[0], bounds[1]);
    }
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

void DivideAndComputeRogerBarlow(TH1D* h1, TH1D *h2){ 
  // Use Roger Barlow "sigma_{delta}" as errors for ratios
  Double_t lh1NBins = h1->GetNbinsX(); 
  Double_t lh2NBins = h2->GetNbinsX(); 

  if( lh1NBins != lh2NBins ){ 
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return;
  }

  Double_t lSigmaDelta[100]; 
  for(Int_t i = 1; i < h1->GetNbinsX() + 1; i++){ 
    // Computation of roger barlow sigma_{delta} 
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
    // Computation of relationship to h2 for plotting in ratio plot 
    if (h2->GetBinContent(i) > 1e-12) { 
      lSigmaDelta[i] /= h2->GetBinContent(i); 
    } else { 
      lSigmaDelta[i] = 0; 
    }
  }
  // Regular Division 
  h1->Divide(h2);
  // Replace Errors
  for(Int_t i = 1; i < h1->GetNbinsX() + 1; i++){
    h1->SetBinError(i, lSigmaDelta[i]);
  }
}

void SubtractAndComputeRogerBarlow(TH1D* h1, TH1D *h2){ 
  // Use Roger Barlow "sigma_{delta}" as errors for differences
  Double_t lh1NBins = h1->GetNbinsX(); 
  Double_t lh2NBins = h2->GetNbinsX(); 

  if( lh1NBins != lh2NBins ){ 
    cout<<"Problem! Number of bins doesn't match! "<<endl;
    return;
  }

  Double_t lSigmaDelta[100]; 
  for(Int_t i = 1; i < h1->GetNbinsX() + 1; i++){ 
    // Computation of roger barlow sigma_{delta} 
    lSigmaDelta[i] = TMath::Sqrt( TMath::Abs( TMath::Power(h1->GetBinError(i),2) - TMath::Power(h2->GetBinError(i),2) ) );
  }
  // Regular Subtraction 
  h1->Add(h2, -1);
  // Replace Errors
  for(Int_t i = 1; i < h1->GetNbinsX() + 1; i++){
    h1->SetBinError(i, lSigmaDelta[i]);
  }
}

void checkSignificance(TH1D* h, TString method = "ratio"){
  // this function takes a 1D histogram that consists of a ratio (or difference) with Roger Barlow errors
  // it integrates the ratio/difference, taking errors into account
  assert(method == "ratio" || method == "difference" && "method must be 'ratio' or 'difference'");

  double error;
  double avg = h->IntegralAndError(1, h->GetNbinsX(), error);
  avg /= h->GetNbinsX(); error /= h->GetNbinsX();
  cout << avg << " +/- " << error << endl;
  double sig;
  if(method == "difference")
    sig = std::abs(avg)/error;
  else if(method == "ratio")
    sig = std::abs(1.-avg)/error;
  cout << "Average " << method << " is " << sig << " sigma away from 1." << endl;
  if (sig < 1)
    cout << "=> Variation is not significant." << endl;
  else 
    cout << "=> Variation is significant!" << endl;
  
  // maybe do some systematics here based on significance?
}