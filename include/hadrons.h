#ifndef HADRONS_
#define HADRONS_

#include <array> // stl array allows for iteration over its elements
#include <unordered_map>
#include "TDataType.h"

// container for all info regarding a specific hadron, including it the info on its antiparticle
class Hadron {
    private:
        //// PROPERTIES ////
        TString name;
        Int_t pdg;
        TString latex;
        Double_t mass; // in MeV
        Int_t strangeness;
        Hadron* antiParticle;

    public:
        //// CONSTRUCTORS ////
        // default constructor, doesn't do anything but allows for declaring instances to be "filled in" later
        Hadron() = default;

        template <class T1, class T2> // use template so that we can use any combination of 'const char *' and 'TString'
        Hadron(T1 _name, const Int_t _pdg, T2 _latex, const Double_t _mass, const Int_t _strangeness)
        {
            name = (TString) _name;
            pdg = _pdg;
            latex = (TString) _latex;
            mass = _mass;
            strangeness = _strangeness;
        }

        void setAntiParticle(Hadron* h) {antiParticle = h;}

        //// GETTERS ////
        TString getName() {return name;}
        Int_t getPDG() {return pdg;}
        TString getLatex() {return latex;}
        Double_t getMass() {return mass;}
        Int_t getStrangeness() {return strangeness;}
        Hadron* getAntiParticle() {return antiParticle;}
};

// declare the hadrons needed for analysis
Hadron Kzero = Hadron("K0", 311, "K^{0}", 497.611, -1);
Hadron Kzerobar = Hadron("K0bar", -311, "#bar{K^{0}}", 497.611, 1);

Hadron Kplus = Hadron("K+", 321, "K^{+}", 493.677, -1);
Hadron Kminus = Hadron("K-", -321, "K^{-}", 493.677, 1);

Hadron Lambda = Hadron("Lambda", 3122, "#Lambda^{0}", 1115.683, 1);
Hadron Lambdabar = Hadron("Lambdabar", -3122, "#bar{#Lambda^{0}}", 1115.683, -1);

Hadron Sigmazero = Hadron("Sigma0", 3212, "#Sigma^{0}", 1192.642, 1);
Hadron Sigmazerobar = Hadron("Sigma0bar", -3212, "#bar{#Sigma^{0}}", 1192.642, -1);

Hadron Sigmaminus = Hadron("Sigma-", 3112, "#Sigma^{-}", 1197.449, 1);
Hadron Sigmaminusbar = Hadron("Sigma-bar", -3112, "#bar{#Sigma^{-}}", 1197.449, -1);

Hadron Sigmaplus = Hadron("Sigma+", 3222, "#Sigma^{+}", 1189.37, 1);
Hadron Sigmaplusbar = Hadron("Sigma+bar", -3222, "#bar{#Sigma^{+}}", 1189.37, -1);

Hadron Xizero = Hadron("Xi0", 3322, "#Xi^{0}", 1314.86, 2);
Hadron Xizerobar = Hadron("Xi0bar", -3322, "#bar{#Xi^{0}}", 1314.86, -2);

Hadron Ximinus = Hadron("Xi-", 3312, "#Xi^{-}", 1321.71, 2);
Hadron Xiplus = Hadron("Xi+", -3312, "#Xi^{+}", 1321.71, -2);

Hadron Omegaminus = Hadron("Omega-", 3334, "#Omega^{-}", 1672.45, 3);
Hadron Omegaplus = Hadron("Omega+", -3334, "#Omega^{+}", 1672.45, -3);

Hadron Kzeroshort = Hadron("K0_S", 310, "K^{0}_{S}", 497.611, 0);
Hadron Kzerolong = Hadron("K0_L", 130, "K^{0}_{L}", 497.611, 0);

// more exotic strange hadrons (life is pain):
// unknown masses are taken from pythia
Hadron Dsplus = Hadron("Ds+", 431, "D_{s}^{+}", 1968.30, -1);
Hadron Dsminus = Hadron("Ds-", -431, "D_{s}^{-}", 1968.30, 1);

Hadron Bszero = Hadron("Bs0", 531, "B_{s}^{0}", 5366.77, 1);
Hadron Bszerobar = Hadron("Bs0bar", -531, "#bar{B_{s}^{0}}", 5366.77, -1);

// these are not strange you dipshit
// Hadron Sigmacplusplus = Hadron("Sigmac++", "Sigmac++bar", 4222, "#Sigma_{c}^{++}", "#bar{#Sigma_{c}^{++}}", 2453.97);
// Hadron Sigmacplus = Hadron("Sigmac+", "Sigmac+bar", 4212, "#Sigma_{c}^{+}", "#bar{#Sigma_{c}^{+}}", 2452.9);
// Hadron Sigmaczero = Hadron("Sigmac0", "Sigmac0bar", 4112, "#Sigma_{c}^{0}", "#bar{#Sigma_{c}^{0}}", 2453.75);
// Hadron Sigmabplus = Hadron("Sigmab+", "Sigmab+bar", 5222, "#Sigma_{b}^{+}", "#bar{#Sigma_{b}^{+}}", 5810.56);
// Hadron Sigmabzero = Hadron("Sigmab0", "Sigmab0bar", 5212, "#Sigma_{b}^{0}", "#bar{#Sigma_{b}^{0}}", 5800.00); // pythia mass
// Hadron Sigmabminus = Hadron("Sigmab-", "Sigmab-bar", 5112, "#Sigma_{b}^{-}", "#bar{#Sigma_{b}^{-}}", 5815.64);

Hadron Xiczero = Hadron("Xic0", 4132, "#Xi_{c}^{0}", 2470.90, 1);
Hadron Xiczerobar = Hadron("Xic0bar", -4132, "#bar{#Xi_{c}^{0}}", 2470.90, -1);

Hadron Xicplus = Hadron("Xic+", 4232, "#Xi_{c}^{+}", 2467.94, 1);
Hadron Xicminus = Hadron("Xic-", -4232, "#Xi_{c}^{-}", 2467.94, -1);

Hadron Xibzero = Hadron("Xib0", 5232, "#Xi_{b}^{0}", 5791.9, 1);
Hadron Xibzerobar = Hadron("Xib0bar", -5232, "#bar{#Xi_{b}^{0}}", 5791.9, -1);

Hadron Xibminus = Hadron("Xib-", 5132, "#Xi_{b}^{-}", 5797.0, 1);
Hadron Xibplus = Hadron("Xib+", -5132, "#Xi_{b}^{+}", 5797.0, -1);

Hadron Omegac = Hadron("Omegac0", 4332, "#Omega_{c}^{0}", 2695.2, 2);
Hadron Omegacbar = Hadron("Omegac0bar", -4332, "#bar{#Omega_{c}^{0}}", 2695.2, -2);

Hadron Omegaccplus = Hadron("Omegacc+", 4432, "#Omega_{cc}^{+}", 3786.63, 1); // pythia mass
Hadron Omegaccminus = Hadron("Omegacc-", -4432, "#Omega_{cc}^{-}", 3786.63, -1); // pythia mass

Hadron Omegabczero = Hadron("Omegabc0", 5342, "#Omega_{bc}^{0}", 7190.99, 1); // pythia mass
Hadron Omegabczerobar = Hadron("Omegabc0bar", -5342, "#bar{#Omega_{bc}^{0}}", 7190.99, -1); // pythia mass

Hadron Omegabminus = Hadron("Omegab-", 5332, "#Omega_{b}^{-}", 6046.1, 2);
Hadron Omegabplus = Hadron("Omegab+", -5332, "#Omega_{b}^{+}", 6046.1, -2);

Hadron Omegabbminus = Hadron("Omegabb-", 5532, "#Omega_{bb}^{-}", 10602.09, 1); // pythia mass
Hadron Omegabbplus = Hadron("Omegabb+", -5532, "#Omega_{bb}^{+}", 10602.09, -1); // pythia mass

// store the hadrons in an iteratable container, with the antiparticle always following the particle for quick linking
const std::array<Hadron*, 40> StrangeHadrons = {&Kzeroshort, &Kzerolong, &Kminus, &Kplus, &Lambda, &Lambdabar, &Sigmaminus, &Sigmaminusbar, &Sigmazero, &Sigmazerobar, &Sigmaplus, &Sigmaplusbar, &Ximinus, &Xiplus, &Xizero, &Xizerobar, &Omegaminus, &Omegaplus, &Dsplus, &Dsminus, &Bszero, &Bszerobar, &Xicplus, &Xicminus, &Xiczero, &Xiczerobar, &Xibminus, &Xibplus, &Xibzero, &Xibzerobar, &Omegac, &Omegacbar, &Omegabminus, &Omegabplus, &Omegaccplus, &Omegaccminus, &Omegabbminus, &Omegabbplus, &Omegabczero, &Omegabczerobar};

// set the antiparticles
int setAntiParticles(){
    for (int i = 0; i < StrangeHadrons.size()/2; i++){
        StrangeHadrons[2*i]->setAntiParticle(StrangeHadrons[2*i+1]);
        StrangeHadrons[2*i+1]->setAntiParticle(StrangeHadrons[2*i]);
    }
    return 0;
}
int dummy = setAntiParticles(); // TODO: find a nicer way to do this

// create a map that maps the pdg to the hadron
std::unordered_map<int, Hadron*> createMap(){
    std::unordered_map<int, Hadron*> map;
    for (Hadron *hadron : StrangeHadrons){
        map[hadron->getPDG()] = hadron;
    }
    return map;
}
const std::unordered_map<int, Hadron*> StrangeHadronPDGMap = createMap();

#endif