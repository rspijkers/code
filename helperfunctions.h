#ifndef HELPERFUNCTIONS_
#define HELPERFUNCTIONS_

#include <cmath> // sinh, cosh, tanh
#include "TFile.h"

double rapidityFromEta(double eta, double pt, double m){
    /*Converts pseudorapidity to rapidity.*/
    double cosheta = cosh(eta);
    double m2 = m*m;
    double pt2 = pt*pt;

    return log((sqrt(m2 + pt2*cosheta*cosheta) + pt*sinh(eta)) / sqrt(m2 + pt2));
}

class NamedFile : public TFile {
    // A class that adds a name to a TFile other than the filepath.
    private:
        TString customName;
    
    public:
        template <class T>
        void SetCustomName(T _name) {customName = _name;}
        TString GetCustomName() {return customName;}
};

#endif