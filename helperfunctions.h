#ifndef HELPERFUNCTIONS_
#define HELPERFUNCTIONS_

// std
#include <iostream> // cout and stuff
#include <cmath> // sinh, cosh, tanh, abs
// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TCanvas.h"
// custom
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

// updates the ranges of h0 to include h + margins
template <class T>
void updateRanges(T* h0, T* h) {
    Double_t xmin, xmax, ymin, ymax, x0min, x0max, y0min, y0max;

    // Xrange
    TAxis* h0X = h0->GetXaxis();
    TAxis* hX = h->GetXaxis();
    xmin = hX->GetXmin();
    x0min = h0X->GetXmin();
    xmax = hX->GetXmax();
    x0max = h0X->GetXmax();
    if(xmin < h0X->GetXmin()) x0min = xmin;
    if(xmax > h0X->GetXmax()) x0max = xmax;
    h0X->SetRangeUser(x0min, x0max);

    // Yrange
    h0->GetMinimumAndMaximum(y0min, y0max);
    h->GetMinimumAndMaximum(ymin, ymax);
    if(ymin < y0min) y0min = ymin;
    if(ymax > y0max) y0max = ymax;
    // margins for the y-range
    y0max = 1.1*y0max + std::sqrt(y0max); // upper margin = 10% + statistical error
    if(y0min != 0) y0min -= std::sqrt(std::abs(y0min)); // lower margin = statistical error
    h0->GetYaxis()->SetRangeUser(y0min, y0max);
}

// // Try to construct a method that lists all unique values of a branch in a tree. 
// void ListUniqueValues(TBranch* branch){
//     branch->G
// }
// // maybe make a copy where the branch is specified by name? (i.e. a string containing the name)

#endif