#ifndef HELPERFUNCTIONS_
#define HELPERFUNCTIONS_

#include <iostream> // cout and stuff
#include <cmath> // sinh, cosh, tanh, abs
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TCanvas.h"
#include "myStyle.h" // check to see if it can find the file first? otherwise the entire header is useless

double rapidityFromEta(double eta, double pt, double m){
    /*Converts pseudorapidity to rapidity.*/
    double cosheta = cosh(eta);
    double m2 = m*m;
    double pt2 = pt*pt;

    return log((sqrt(m2 + pt2*cosheta*cosheta) + pt*sinh(eta)) / sqrt(m2 + pt2));
}

// function that determines the strangeness of a hadron given its pdg code. 
// returns the absolute value of the strangeness (i.e. 2 for ssbar mesons, instead of net zero)
int strangenessFromPDG(const int _pdg){
    int strangeness = 0;
    int pdg = std::abs(_pdg);
    pdg /= 10;
    if (pdg % 10 == 3) strangeness++; // 3rd quark
    pdg /= 10;
	if (pdg % 10 == 3) strangeness++; // 2nd quark
	pdg /= 10;
	if (pdg % 10 == 3) strangeness++; // 1st quark
    return strangeness;
}

// returns the value of the nth digit, counted from right to left starting at 0
// example: getDigitN(123, 0) will return 3
int getDigitN(const int target, int n){
    // create 10^n so we can divide the target by it. 
    int En = 1; 
    while (n--) En *= 10;
    return (target/En) % 10; // divide by 10^n, and return mod % 10 yielding the content of the digit at the nth place
}

// A class that adds a name to a TFile other than the filepath.
// Note that TFile has a title option built in, but this does not work for some reason
class NamedFile : public TFile {
    private:
        TString customName;
    
    public:
        // This constructor passes filename and option to the TFile constructor, while explicitly setting customName
        NamedFile(const char* filename, const char* _customName, Option_t* option = "") : TFile(filename, option) {
            customName = _customName;
        }

        template <class T>
        void SetCustomName(T _name) {customName = _name;}
        TString GetCustomName() {return customName;}
};

// Make a function that creates a default canvas complete with axes, labels, everything
// This canvas can then be used in plots with multiple histo's, simply by doing h->Draw("same")
template <class T> // can be TH1D, TH1F, TH1I, etc.
TCanvas* createMultiTH1(std::vector<T*> TH1vector, Bool_t useMyStyle=true){
    if(useMyStyle){
        customStyle::myStyle->cd(); // set custom style
        gROOT->ForceStyle(); // force custom style on objects created with a different style
    }
    Double_t xmin, xmax, ymin, ymax;
    // Does this also take error bars into account?
    int i = 0;
    for(T* h : TH1vector){ // iteration 0 over all histo's, to determine x and y min and max
        Double_t _xmin, _xmax, _ymin, _ymax;
        h->GetMinimumAndMaximum(_ymin, _ymax);
        _xmin = h->GetXaxis()->GetXmin();
        _xmax = h->GetXaxis()->GetXmax();
        if(i == 0){ // first iteration always set the values
            xmin = _xmin;
            xmax = _xmax;
            ymin = _ymin;
            ymax = _ymax;
        } else { // not the first iteration, check per value
            if(_xmin < xmin) xmin = _xmin;
            if(_xmax > xmax) xmax = _xmax;
            if(_ymin < ymin) ymin = _ymin;
            if(_ymax > ymax) ymax = _ymax;
        }
        i++;
    }
    ymax = 1.1*ymax + std::sqrt(ymax); // y margin statistical error
    if(ymin != 0) ymin -= std::sqrt(std::abs(ymin)); // statistical error
    TCanvas* canvas = new TCanvas("multiplotname", "multiplottitle", 800, 600);
    canvas->cd();
    TH1* frame = canvas->DrawFrame(xmin, ymin, xmax, ymax);
    
    for(T* h : TH1vector){
        h->Draw("E same");
    }
    return canvas; //???
}

// // Try to construct a method that lists all unique values of a branch in a tree. 
// void ListUniqueValues(TBranch* branch){
//     branch->G
// }
// // maybe make a copy where the branch is specified by name? (i.e. a string containing the name)

#endif