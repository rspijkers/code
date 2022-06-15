#ifndef HELPERFUNCTIONS_
#define HELPERFUNCTIONS_

#include <iostream> // cout and stuff
#include <cmath> // sinh, cosh, tanh
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

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

// // Try to construct a method that lists all unique values of a branch in a tree. 
// void ListUniqueValues(TBranch* branch){
//     branch->G
// }
// // maybe make a copy where the branch is specified by name? (i.e. a string containing the name)

#endif