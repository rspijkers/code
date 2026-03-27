// Analysis macro for Xi-Xi trigger-normalised delta-phi correlations.
// Reads the tree produced by PythiaEventGen/cascade_correlations.cpp.
//
// Usage: ./xi_xi_correlations inputfile(s).root [outputfile.root]
//   inputfile(s) can be a glob, e.g. "data/*.root"
//   outputfile defaults to xi_xi_correlations.root
//
// Kinematic selection: 1 < pT < 8 GeV/c, |eta| < 0.8
// Produces trigger-normalised delta-phi histograms for same-sign and
// opposite-sign Xi-Xi pairs. To avoid charge bias (unequal Xi- vs Xi+
// populations), each trigger charge is normalised by its own trigger count
// and the two are averaged before writing.

// std
#include <iostream>
#include <vector>
#include <cmath>
// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TChain.h"
#include "TSystem.h"
// custom
#include "include/helperfunctions.h"
#include "include/SmallEvent.h"

using std::cout; using std::cerr; using std::endl;

int main(int argc, char** argv)
{
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " inputfile(s).root [outputfile.root]" << endl;
    return 1;
  }

  const char* inputFiles = argv[1];
  const char* outputPath = (argc >= 3) ? argv[2] : "xi_xi_correlations.root";

  gSystem->Load("lib/libEvent.so");
  TH1::SetDefaultSumw2();

  TChain* chain = new TChain("tree");
  chain->Add(inputFiles);

  SmallEvent* event = new SmallEvent();
  chain->SetBranchAddress("event", &event);

  // --- Kinematic cuts ---
  const Double_t pTmin  = 1.0;
  const Double_t pTmax  = 8.0;
  const Double_t etaMax = 0.8;
  const Int_t    xiNeg  =  3312; // Xi-
  const Int_t    xiPos  = -3312; // Xi+ (antiXi)

  // Four histograms split by trigger charge, so we can normalise each
  // by its own trigger count before combining (avoids charge-population bias).
  const Int_t  nBins  = 64;
  const Double_t lo   = -0.5*PI;
  const Double_t hi   =  1.5*PI;
  const char*  xTitle = "#Delta#varphi (rad)";
  const char*  yTitle = "(1/N_{trig}) dN/d#Delta#varphi";

  TH1D* hDPhiXiMinXiMin = new TH1D("hDPhiXiMinXiMin", TString::Format("#Xi^{-} trig, #Xi^{-} assoc;%s;%s", xTitle, yTitle), nBins, lo, hi);
  TH1D* hDPhiXiMinXiPlus = new TH1D("hDPhiXiMinXiPlus", TString::Format("#Xi^{-} trig, #Xi^{+} assoc;%s;%s", xTitle, yTitle), nBins, lo, hi);
  TH1D* hDPhiXiPlusXiPlus = new TH1D("hDPhiXiPlusXiPlus", TString::Format("#Xi^{+} trig, #Xi^{+} assoc;%s;%s", xTitle, yTitle), nBins, lo, hi);
  TH1D* hDPhiXiPlusXiMin = new TH1D("hDPhiXiPlusXiMin", TString::Format("#Xi^{+} trig, #Xi^{-} assoc;%s;%s", xTitle, yTitle), nBins, lo, hi);

  Double_t nTrig_neg = 0;
  Double_t nTrig_pos = 0;

  const Long64_t nEvents = chain->GetEntries();
  cout << "Processing " << nEvents << " events..." << endl;

  for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
    chain->GetEntry(iEvent);
    const std::vector<SmallTrack> cands = event->getCandidates();

    // collect Xi candidates passing kinematic cuts
    std::vector<SmallTrack> xi;
    for (const SmallTrack& t : cands) {
      if (std::abs(t.getPDG()) != std::abs(xiNeg)) continue;
      if (t.getpT() < pTmin || t.getpT() > pTmax)  continue;
      if (std::abs(t.getEta()) > etaMax)             continue;
      xi.push_back(t);
    }

    const Int_t n = xi.size();

    // count triggers per charge
    for (Int_t i = 0; i < n; i++) {
      if (xi[i].getPDG() == xiNeg) nTrig_neg++;
      else                         nTrig_pos++;
    }

    // trigger-associate double loop; require pT_trig > pT_assoc
    for (Int_t i = 0; i < n; i++) {
      const Int_t pdgTrig = xi[i].getPDG();
      for (Int_t j = 0; j < n; j++) {
        if (i == j) continue;
        if (xi[i].getpT() <= xi[j].getpT()) continue; // enforce pT ordering
        const Int_t pdgAssoc = xi[j].getPDG();
        const Double_t dphi  = DeltaPhi(xi[i].getPhi(), xi[j].getPhi());
        const Bool_t sameSign = (pdgTrig == pdgAssoc);

        if (pdgTrig == xiNeg) {
          if (sameSign) hDPhiXiMinXiMin->Fill(dphi);
          else          hDPhiXiMinXiPlus->Fill(dphi);
        } else {
          if (sameSign) hDPhiXiPlusXiPlus->Fill(dphi);
          else          hDPhiXiPlusXiMin->Fill(dphi);
        }
      }
    }
  } // event loop

  cout << "Xi- triggers: " << nTrig_neg << "  Xi+ triggers: " << nTrig_pos << endl;

  if (nTrig_neg == 0 || nTrig_pos == 0) {
    cerr << "Warning: one trigger charge has zero counts. Check kinematic cuts." << endl;
    return 1;
  }

  // normalise each histogram by its own trigger count
  hDPhiXiMinXiMin->Scale(1. / nTrig_neg);
  hDPhiXiMinXiPlus->Scale(1. / nTrig_neg);
  hDPhiXiPlusXiPlus->Scale(1. / nTrig_pos);
  hDPhiXiPlusXiMin->Scale(1. / nTrig_pos);

  // combine trigger charges: add the two charge-unbiased distributions
  TH1D* hDPhi_SS = (TH1D*) hDPhiXiMinXiMin->Clone("hDPhi_SS");
  hDPhi_SS->SetTitle(TString::Format("#Xi#Xi same-sign #Delta#varphi;%s;%s", xTitle, yTitle));
  hDPhi_SS->Add(hDPhiXiPlusXiPlus);
  hDPhi_SS->Scale(0.5);

  TH1D* hDPhi_OS = (TH1D*) hDPhiXiMinXiPlus->Clone("hDPhi_OS");
  hDPhi_OS->SetTitle(TString::Format("#Xi#bar{#Xi} opposite-sign #Delta#varphi;%s;%s", xTitle, yTitle));
  hDPhi_OS->Add(hDPhiXiPlusXiMin);
  hDPhi_OS->Scale(0.5);

  TFile* outFile = new TFile(outputPath, "RECREATE");
  hDPhi_SS->Write();
  hDPhi_OS->Write();
  // also save the per-charge histograms for cross-checks
  hDPhiXiMinXiMin->Write();
  hDPhiXiMinXiPlus->Write();
  hDPhiXiPlusXiPlus->Write();
  hDPhiXiPlusXiMin->Write();
  outFile->Close();

  cout << "Output written to " << outputPath << endl;
  return 0;
}
