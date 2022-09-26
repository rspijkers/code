#ifndef SMALLTRACK_H
#define SMALLTRACK_H
#include "TClass.h"
#include <vector>

class SmallTrack {
    private:
        Int_t pdg; // pdg code
        Double_t pT;
        Double_t eta;
        Double_t phi;

    public:
        SmallTrack();
        ~SmallTrack();
        void Clear();

        // setters
        void setPDG(Int_t _pdg) {pdg = _pdg;}
        void setpT(Double_t _pT) {pT = _pT;}
        void setEta(Double_t _eta) {eta = _eta;}
        void setPhi(Double_t _phi) {phi = _phi;}

        // getters
        Int_t getPDG() const {return pdg;}
        Double_t getpT() const {return pT;}
        Double_t geteta() const {return eta;}
        Double_t getphi() const {return phi;}

        ClassDef(SmallTrack,1); 
};
#endif