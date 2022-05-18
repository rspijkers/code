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
Hadron Kzero = Hadron("K0", "K0bar", 311, "K^{0}", "#bar{K^{0}}", 497.611);
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
#endif