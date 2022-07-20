#include "helperfunctions.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include "TRandom.h"

int helperfunctionstester(){
    TH1D* hGaussian1 = new TH1D("hGaussian1", "A Gaussian distribution", 100, -10, 10);
    TH1D* hGaussian2 = new TH1D("hGaussian2", "A Gaussian distribution", 50, 0, 10);

    gRandom->SetSeed(0);
    TF1* fGaussian1 = new TF1("fGaussian1", "gaus");
    TF1* fGaussian2 = new TF1("fGaussian2", "gaus");
    fGaussian1->SetParameters(1,0,2);
    fGaussian2->SetParameters(1,5,1);

    hGaussian1->FillRandom("fGaussian1", 5000);
    hGaussian2->FillRandom("fGaussian2", 5000);
    std::vector<TH1D*> input = {hGaussian1, hGaussian2};

    TCanvas* test = new TCanvas("name", "title");
    TH1D* h = new TH1D("hname", "htitle", 10, 0 ,10);
    h->SetTitle("#Lambda^{-} title");
    h->Fill(1);
    h->Draw();

    return 0;
}