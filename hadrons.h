#pragma once
#ifndef HADRONS_
#define HADRONS_
// container for all info regarding a specific hadron, including it the info on its antiparticle
class Hadron {
    private:
        //// PROPERTIES ////
        TString name;
        TString antiname;
        Int_t pdg;
        Int_t antipdg;
        TString latex;
        TString antilatex;
        Double_t mass; // in MeV

    public:
        //// CONSTRUCTORS ////
        // default constructor, doesn't do anything but allows for declaring instances to be "filled in" later
        Hadron() = default;

        // constructor for hadrons with a well defined anti-particle, such as the charged Kaon
        template <class T1, class T2, class T3, class T4> // use template so that we can use any combination of 'const char *' and 'TString'
        Hadron(T1 _name, T2 _antiname, const Int_t _pdg, T3 _latex, T4 _antilatex, const Double_t _mass)
        {
            name = (TString) _name;
            antiname = (TString) _antiname;
            pdg = _pdg;
            antipdg = -1*pdg;
            latex = (TString) _latex;
            antilatex = (TString) _antilatex;
            mass = _mass;
        }

        // constructor for hadrons without a well defined anti-particle, such as the neutral K0_short/long
        template <class T1, class T2>
        Hadron(T1 _name, Int_t _pdg, T2 _latex, const Double_t _mass)
        {
            name = (TString) _name;
            pdg = _pdg;
            latex = (TString) _latex;
            mass = _mass;
        }

        //// GETTERS ////
        TString getName() {return name;}
        TString getAntiName() {return antiname;}
        Int_t getPDG() {return pdg;}
        Int_t getAntiPDG() {return antipdg;}
        TString getLatex() {return latex;}
        TString getAntiLatex() {return antilatex;}
        Double_t getMass() {return mass;}
};


// declare the hadrons needed for analysis
// Hadron Kzero = Hadron("K0", "K0bar", 311, "K^{0}", "#bar{K^{0}}", 497.611); // this is an antistrange hadron!
Hadron Kzerobar = Hadron("K0bar", "K0", -311, "#bar{K^{0}}", "K^{0}", 497.611);
Hadron Kminus = Hadron("K-", "K+", -321, "K^{-}", "K^{+}", 493.677);
Hadron Lambda = Hadron("Lambda", "Lambdabar", 3122, "#Lambda", "#bar{#Lambda}", 1115.683);
Hadron Sigmazero = Hadron("Sigma0", "Sigma0bar", 3212, "#Sigma^{0}", "#bar{#Sigma^{0}}", 1192.642);
Hadron Sigmaminus = Hadron("Sigma-", "Sigma-bar", 3112, "#Sigma^{-}", "#bar{#Sigma^{-}}", 1197.449);
Hadron Sigmaplus = Hadron("Sigma+", "Sigma+bar", 3222, "#Sigma^{+}", "#bar{#Sigma^{+}}", 1189.37);
Hadron Xizero = Hadron("Xi0", "Xi0bar", 3322, "#Xi^{0}", "#bar{#Xi^{0}}", 1314.86);
Hadron Ximinus = Hadron("Xi-", "Xi+", 3312, "#Xi^{-}", "#Xi^{+}", 1321.71);
Hadron Omegaminus = Hadron("Omega-", "Omega+", 3334, "#Omega^{-}", "#Omega^{+}", 1672.45);
Hadron Kzeroshort = Hadron("K0_S", 310, "K^{0}_{S}", 497.611);
Hadron Kzerolong = Hadron("K0_L", 130, "K^{0}_{L}", 497.611);

// more exotic strange hadrons (life is pain):
// unknown masses are taken from pythia
Hadron Dsubs = Hadron("Ds-", "Ds+", -431, "D_{s}^{-}", "D_{s}^{+}", 1968.30);
Hadron Bsubs = Hadron("Bs0", "Bs0bar", 531, "B_{s}^{0}", "#bar{B_{s}^{0}}", 5366.77);

Hadron Sigmacplusplus = Hadron("Sigmac++", "Sigmac++bar", 4222, "#Sigma_{c}^{++}", "#bar{#Sigma_{c}^{++}}", 2453.97);
Hadron Sigmacplus = Hadron("Sigmac+", "Sigmac+bar", 4212, "#Sigma_{c}^{+}", "#bar{#Sigma_{c}^{+}}", 2452.9);
Hadron Sigmaczero = Hadron("Sigmac0", "Sigmac0bar", 4112, "#Sigma_{c}^{0}", "#bar{#Sigma_{c}^{0}}", 2453.75);
Hadron Sigmabplus = Hadron("Sigmab+", "Sigmab+bar", 5222, "#Sigma_{b}^{+}", "#bar{#Sigma_{b}^{+}}", 5810.56);
Hadron Sigmabzero = Hadron("Sigmab0", "Sigmab0bar", 5212, "#Sigma_{b}^{0}", "#bar{#Sigma_{b}^{0}}", 5800.00); // pythia mass
Hadron Sigmabminus = Hadron("Sigmab-", "Sigmab-bar", 5112, "#Sigma_{b}^{-}", "#bar{#Sigma_{b}^{-}}", 5815.64);

Hadron Xiczero = Hadron("Xic0", "Xic0bar", 4132, "#Xi_{c}^{0}", "#bar{#Xi_{c}^{0}}", 2470.90);
Hadron Xicplus = Hadron("Xic+", "Xic-", 4232, "#Xi_{c}^{+}", "#Xi_{c}^{-}", 2467.94);
Hadron Xibzero = Hadron("Xib0", "Xib0bar", 5232, "#Xi_{b}^{0}", "#bar{#Xi_{b}^{0}}", 5791.9);
Hadron Xibmin = Hadron("Xib-", "Xib+", 5132, "#Xi_{b}^{-}", "#Xi_{b}^{+}", 5797.0);

Hadron Omegac = Hadron("Omegac0", "Omegac0bar", 4332, "#Omega_{c}^{0}", "#bar{#Omega_{c}^{0}}", 2695.2);
Hadron Omegacc = Hadron("Omegacc+", "Omegacc-", 4432, "#Omega_{cc}^{+}", "#Omega_{cc}^{-}", 3786.63); // pythia mass

Hadron Omegabc = Hadron("Omegabc0", "Omegabc0bar", 5342, "#Omega_{bc}^{0}", "#bar{#Omega_{bc}^{0}", 7190.99); // pythia mass

Hadron Omegab = Hadron("Omegab-", "Omegab+", 5332, "#Omega_{b}^{-}", "#Omega_{b}^{+}", 6046.1);
Hadron Omegabb = Hadron("Omegabb-", "Omegabb+", 5532, "#Omega_{bb}^{-}", "#Omega_{bb}^{+}", 10602.09); // pythia mass

#endif