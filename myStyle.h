#include "TStyle.h"
/* 
This header file creates a custom style named "myStyle". This style is not by default selected.
Use `gROOT->SetStyle("myStyle");` or `myStyle->cd();` to set the current style.
After this, `gStyle` will point to `myStyle`, allowing the user to further customize the style
by `gStyle->Set...` without having to edit this header
*/

TStyle* createStyle(Bool_t setstyle = false){
    TStyle* myStyle = new TStyle("myStyle", "A custom style for ROOT figures");

    myStyle->Reset("Plain");
    myStyle->SetOptTitle(0);
    myStyle->SetOptStat(0);
    myStyle->SetPalette(1); // TODO: better palette's are available, such as kBird?
    myStyle->SetCanvasColor(10);
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetFrameLineWidth(1);
    myStyle->SetFrameFillColor(kWhite);
    myStyle->SetPadColor(10);
    myStyle->SetPadTickX(1);
    myStyle->SetPadTickY(1);
    myStyle->SetPadBottomMargin(0.15);
    myStyle->SetPadLeftMargin(0.15);
    myStyle->SetHistLineWidth(1);
    myStyle->SetHistLineColor(kRed);
    myStyle->SetFuncWidth(2);
    myStyle->SetFuncColor(kGreen);
    myStyle->SetLineWidth(1);
    myStyle->SetLabelSize(0.045,"xyz");
    myStyle->SetLabelOffset(0.01,"y");
    myStyle->SetLabelOffset(0.01,"x");
    myStyle->SetLabelColor(kBlack,"xyz");
    myStyle->SetTitleSize(0.05,"xyz");
    myStyle->SetTitleOffset(1.25,"y");
    myStyle->SetTitleOffset(1.2,"x");
    myStyle->SetTitleFillColor(kWhite);
    myStyle->SetTextSizePixels(26);
    myStyle->SetTextFont(42);
    //  myStyle->SetTickLength(0.04,"X");  myStyle->SetTickLength(0.04,"Y"); 
    myStyle->SetLegendBorderSize(1);
    myStyle->SetLegendFillColor(kWhite);
    //  myStyle->SetFillColor(kWhite);
    myStyle->SetLegendFont(42);

    if(setstyle) myStyle->cd();

    return myStyle;
}

// put custom styles in a namespace
namespace customStyle{
    TStyle* myStyle = createStyle();
    const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7}; // for syst bands
    const Int_t colors[]     = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2};
    const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};
}