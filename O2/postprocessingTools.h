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
