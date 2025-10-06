#include "postprocessingTools.h"
#include "TRatioPlot.h"

void plotMCClosure(TString runnr) {
  TFile *fData = new TFile("plots/"+runnr+".root", "READ");
  TFile *fMC = new TFile("plots/Closure"+runnr+".root", "READ");
  
  if(!fData){
    cout << "data file is null! aborting..." << endl;
    return;
  }
  if(!fMC){
    cout << "MC file is null! aborting..." << endl;
    return;
  }

  int nplots = 3;
  TString dataHistNames[3] = {"hXiXiOS", "hXiXiSS", "hXiXisubtracted"};
  TString mcHistNames[3] = {"hOS", "hSS", "hSubtracted"};

  TCanvas *c = new TCanvas("c");
  c->cd();

  for (int i = 0; i < nplots; i++){
    TH1D *hData, *hMC;
    fData->GetObject(dataHistNames[i], hData);
    fMC->GetObject(mcHistNames[i], hMC);

    hData->SetLineColor(kRed);
    hData->SetStats(false);

    TRatioPlot *hRatio = new TRatioPlot(hData, hMC, "divsym");
    hRatio->SetH1DrawOpt("E");
    hRatio->Draw("nogrid");
    hRatio->GetUpperRefYaxis()->SetRangeUser(0., 0.035);
    hRatio->GetLowerRefYaxis()->SetRangeUser(0., 2.);
    TLegend *l = new TLegend(0.55, 0.75, 0.75, 0.85);
    l->AddEntry(hData, "data");
    l->AddEntry(hMC, "MC truth");
    l->Draw();
    c->Print("figures/closure/"+mcHistNames[i]+".pdf");
    c->Clear();
  }

}

int plottingMacro(TString runnr){
  plotMCClosure(runnr);

	return 0;
}

// hRun2OS->SetName("hRun2OS");
// hRun2OS->SetTitle("Run 2 OS #Xi - #Xi");
// hRun2OS->SetContent(datapointsOS);
// hRun2OS->SetError(staterrorsOS);
// hRun2OS->SetStats(false);
// hRun2OS->SetLineColor(kRed);
// TRatioPlot *hRun2OSRatio = new TRatioPlot(hXiXiOS, hRun2OS, "divsym");
// hRun2OSRatio->SetH1DrawOpt("E");
// hRun2OSRatio->Draw("nogrid");
// hRun2OSRatio->GetUpperRefYaxis()->SetRangeUser(0., 0.02);
// hRun2OSRatio->GetLowerRefYaxis()->SetRangeUser(0., 1.);
// hRun2OSRatio->GetLowYaxis()->SetNdivisions(505);
// TLegend *legendOS = new TLegend(0.55, 0.75, 0.75, 0.85);
// hRun2OSRatio->GetUpperPad()->cd(); // draw legend in upper pad
// legendOS->AddEntry(hXiXiOS, "this analysis");
// legendOS->AddEntry(hRun2OS, "run 2");
// legendOS->Draw();
// c->Write("cOSRatio");
// if(makePDF) c->Print("figures/OSRatio.pdf");
// c->Clear();