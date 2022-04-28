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

    public:
        //// CONSTRUCTORS ////
        // constructor for hadrons with a well defined anti-particle, such as the charged Kaon
        template <class T1, class T2, class T3, class T4> // use template so that we can use any combination of 'const char *' and 'TString'
        Hadron(T1 _name, T2 _antiname, Int_t _pdg, T3 _latex, T4 _antilatex)
        {
            name = (TString) _name;
            antiname = (TString) _antiname;
            pdg = _pdg;
            antipdg = -1*pdg;
            latex = (TString) _latex;
            antilatex = (TString) _antilatex;
        }

        // constructor for hadrons without a well defined anti-particle, such as the neutral K0_short/long
        template <class T1, class T2>
        Hadron(T1 _name, Int_t _pdg, T2 _latex)
        {
            name = (TString) _name;
            pdg = _pdg;
            latex = (TString) _latex;
        }

        //// GETTERS ////
        TString getName() {return name;}
        TString getAntiName() {return antiname;}
        Int_t getPDG() {return pdg;}
        Int_t getAntiPDG() {return antipdg;}
        TString getLatex() {return latex;}
        TString getAntiLatex() {return antilatex;}
};


// declare the hadrons needed for analysis, use pointers so we only create one instance per hadron -- no need for copying
static Hadron Kzero = Hadron("K0", "K0bar", 311, "K^{0}", "#bar{K^{0}}");
static Hadron Kminus = Hadron("K-", "K+", -321, "K^{-}", "K^{+}");
static Hadron Lambda = Hadron("Lambda", "Lambdabar", 3122, "#Lambda", "#bar{#Lambda}");
static Hadron Sigmazero = Hadron("Sigma0", "Sigma0bar", 3212, "#Sigma^{0}", "#bar{#Sigma^{0}}");
static Hadron Sigmaminus = Hadron("Sigma-", "Sigma-bar", 3112, "#Sigma^{-}", "#bar{#Sigma^{-}}");
static Hadron Sigmaplus = Hadron("Sigma+", "Sigma+bar", 3222, "#Sigma^{+}", "#bar{#Sigma^{+}}");
static Hadron Xizero = Hadron("Xi0", "Xi0bar", 3322, "#Xi^{0}", "#bar{#Xi^{0}}");
static Hadron Ximinus = Hadron("Xi-", "Xi+", 3312, "#Xi^{-}", "#Xi^{+}");
static Hadron Omegaminus = Hadron("Omega-", "Omega+", 3334, "#Omega^{-}", "#Omega^{+}");
static Hadron Kzeroshort = Hadron("K0_S", 310, "K^{0}_{S}");
static Hadron Kzerolong = Hadron("K0_L", 130, "K^{0}_{L}");
#endif