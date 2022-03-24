#!/usr/bin/env python3

# NEVER IMPORT ROOT BEFORE NUMPY!!!
# incompatibilities may lead to segmentation violation
from hadrons import * # import hadron classes to easily access well defined properties such as pdg codes.
import numpy as np
from joblib import Parallel, delayed
from ROOT import TFile, TH1, TH1D, gStyle
# import os#, sys

# function that iterates over axes in THn, prints names and indices,
def produce_axes_dict(THn):
    # TODO: add function description
    """
    """
    axes_array = THn.GetListOfAxes()
    print(f"Number of axes: {axes_array.GetEntries()}")
    axesdict = {}; axisindex = 0
    for axis in axes_array:
        axisname = axis.GetName()
        print(f"axis {axisindex} is named {axisname}")
        axesdict[axisname] = axisindex
        axisindex += 1
    return axesdict

# function that creates the projections from the THn
def project(THn, *axes, **kwargs):
    # TODO: add function description
    """
    """

    assert (THn.ClassName() == 'THnSparseT<TArrayD>'), "Error in project(): input histogram does not match the type 'THnSparseD' (THnSparseT<TArrayD>)!"

    n = len(axes) # more than 3D is complicated and not necessary
    assert (0 < n < 4), "Error in project(): number of axes is not 1, 2, or 3! Aborting function..."

    THn = THn.Clone()

    # if trigger and/or associated pdg's are specified in kwargs, set the ranges BEFORE projection
    trigger = kwargs.get('trigger') 
    if trigger is not None:
        THn.GetAxis(0).SetRangeUser(trigger.pdg, trigger.pdg + 1) # FIXME: assumes pdgTrigger is axis 0
    assoc = kwargs.get('associated')
    if assoc is not None:
        THn.GetAxis(1).SetRangeUser(assoc.pdg, assoc.pdg + 1) # FIXME: assumes pdgAssoc is axis 1

    # set range of continuous axis (i.e. not pdg, so eta or pT or whatever)
    setrange = kwargs.get('setrange')
    if setrange is not None:
        assert len(setrange[0]) == 3, "the first element in setrange is not a list of length 3! You are suspected to give a list of elements containing [axis, min, max]. If you only want to specify the range of one axis, you need to nest it like [[axis, min, max]]"
        for [axis,amin,amax] in setrange:
            THn.GetAxis(axis).SetRangeUser(amin, amax)

    # make the projection now
    if n == 2: # ROOT's 2D projection switches x and y axes for some fucking reason
        projection = THn.Projection(axes[1], axes[0])
    else:
        projection = THn.Projection(*axes)

    # pass optional kwargs to output histo (projection)
    title = kwargs.get('title')
    if title is None: # default title
        title = "Projection of "
        for i in axes: title += THn.GetAxis(i).GetName() + ", "
        title += "with "
        if trigger is not None: title += "trigger " + trigger.latex + ", "
        if assoc is not None: title += "associated " + assoc.latex
    projection.SetTitle(title)

    name = kwargs.get('name')
    if name is None:
        name = f"h{n}"
        if assoc is not None: name += assoc.name
        for i in axes: name += "_" + THn.GetAxis(i).GetName()
    projection.SetName(name)

    range_user = kwargs.get('plotrange')
    if range_user is not None:
        # check for correct dimensionality of range:
        range_user = np.array(range_user)
        rangedim = np.shape(range_user)
        if n == 1 and rangedim == (2,): # exception for the 1D case
            projection.GetXaxis().SetRangeUser(range_user[0], range_user[1])
        else:
            assert (rangedim == (n, 2)), f"Warning in project(): could not apply the range specified by user. Expected array of shape ({n}, 2), found {rangedim}"
            # manually set the ranges, since ROOT is too dumb to understand GetAxis(n)
            # don't ask me why...
            projection.GetXaxis().SetRangeUser(range_user[0, 0], range_user[0, 1])
            if n >= 2:
                projection.GetYaxis().SetRangeUser(range_user[1, 0], range_user[1, 1])
                if n >= 3:
                    projection.GetZaxis().SetRangeUser(range_user[2, 0], range_user[2, 1])

    # other options? draw options (colz)?

    # TODO: add Legend

    # cleanup to prevent RAM crash
    THn.Delete()

    return projection

# make all relevant deltaPhi projections
def make_projections(inputfilepath, parallel = False):
    print("Started make_projections()...")

    inputfile = TFile(inputfilepath, "READ") # READ should be default
    if inputfile.IsZombie():
        print(f"Error: Failed to open file {inputfilepath}")
        raise Exception('IsZombie')
    hTHn = inputfile.Get("hSS")
        
    # get the trigger PDG and pT vs eta histo's to propagate to projections outputfile
    hPDG = inputfile.Get("hPDG").Clone()
    hEtaPt = inputfile.Get("hEtaPt").Clone()

    axes = produce_axes_dict(hTHn)

    # TODO: set global plot settings? smt like set_plotstyle()
    TH1.AddDirectory(False);

    outputfile = TFile("output/Projections.root", "RECREATE")
    # propagate plots to projections output, and close inputfile
    hPDG.SetDirectory(outputfile); hPDG.Write()
    hEtaPt.SetDirectory(outputfile); hEtaPt.Write()
    inputfile.Close()

    etarange = [axes['etaTrigger'], -3., 3.] # eta cut on the trigger, this way we don't have "half jets"
    # make the deltaPhi correlation histo's
    if parallel:
        def write_corr(trig, i, assoc, j):
            trigbar = antistrangelist[i]
            assocbar = strangelist[j]
            histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = assoc, name = assoc.name, setrange = [etarange])
            histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trigbar, associated = assocbar, setrange = [etarange]))
            return histo.Write()
        def write_ss_corr(trig, i, assoc, j):
            trigbar = antistrangelist[i]
            assocbar = antistrangelist[j]
            histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = assoc, name = assoc.name, setrange = [etarange])
            histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trigbar, associated = assocbar, setrange = [etarange]))
            return histo.Write()

        for i,trig in enumerate(strangelist):
            # opposite sign correlations
            x = outputfile.mkdir(trig.name, f"strangeness correlations with {trig.name} trigger")
            x.cd()
            Parallel(n_jobs = len(strangelist), require='sharedmem')(delayed(write_corr)(trig, i, assoc, j) for j,assoc in enumerate(antistrangelist))

            # K0_S/L exception - 0 (net) strangeness so opposite/same sign doesn't make sense
            # We do have to count them though, since they are strange hadrons
            trigbar = antistrangelist[i]
            histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = Kzeroshort, name = Kzeroshort.name, setrange = [etarange])
            histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trigbar, associated = Kzeroshort, setrange = [etarange]))
            histo.Write()
            histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = Kzerolong, name = Kzerolong.name, setrange = [etarange])
            histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trigbar, associated = Kzerolong, setrange = [etarange]))
            histo.Write()

            # same sign correlations
            x = outputfile.mkdir(f"{trig.name}ss", f"samesign strangeness correlations with {trig.name} trigger")
            x.cd()
            Parallel(n_jobs = len(strangelist), require='sharedmem')(delayed(write_ss_corr)(trig, i, assoc, j) for j,assoc in enumerate(strangelist))
    elif not parallel: # TODO: add K0_S/L assoc correlations here
        for i,trig in enumerate(strangelist):
            trigbar = antistrangelist[i]
            # opposite sign correlations
            x = outputfile.mkdir(trig.name, f"strangeness correlations with {trig.name} trigger")
            x.cd()
            for j,assoc in enumerate(antistrangelist): 
                assocbar = strangelist[j]
                histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = assoc, name = assoc.name, setrange = [etarange])
                histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trigbar, associated = assocbar, setrange = [etarange]))
                histo.Write()
                
            # same sign correlations
            x = outputfile.mkdir(f"{trig.name}ss", f"samesign strangeness correlations with {trig.name} trigger")
            x.cd()
            for j,assoc in enumerate(strangelist): 
                assocbar = antistrangelist[j]
                histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = assoc, name = assoc.name, setrange = [etarange])
                histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trigbar, associated = assocbar, setrange = [etarange]))
                histo.Write()

    # outputfile.cd()    # return to parent
    # outputfile.Write()

    # histo = project(hTHn, axes['etaTrigger'], name = 'etaproject')
    # total = histo.Integral()
    # histo.GetXaxis().SetRangeUser(-3.5, 3.5)
    # etapart = histo.Integral()
    # print(float(etapart)/float(total))
    # histo = project(hTHn, axes['etaTrigger'], trigger = Ximinus, name = 'etaprojectXi')
    # total = histo.Integral()
    # histo.GetXaxis().SetRangeUser(-3.5, 3.5)
    # etapart = histo.Integral()
    # print(float(etapart)/float(total))

    outputfile.Close()
    print("Finished make_projections()")
    return 

def analysis(inputpath):
    print("Started analysis()...")
    inputfile = TFile(inputpath, "READ")
    TH1.AddDirectory(False);

    hPDG = inputfile.Get("hPDG")
    hEtaPt = inputfile.Get("hEtaPt").Clone()
    hEtaPt.GetYaxis().SetRangeUser(4, 20)
    total = hEtaPt.Integral()
    hEtaPt.GetXaxis().SetRangeUser(-3., 3.)
    etapart = hEtaPt.Integral()
    factor = float(etapart)/float(total)
    # factor = .75

    outputfile = TFile("output/AnalysisResults.root", "RECREATE")

    ndims = len(strangelist)
    # primitively substract samesign correlations from opposite sign
    for i,trig in enumerate(strangelist):
        trigbar = antistrangelist[i]
        # get number of triggers
        trigbin = hPDG.FindBin(trig.pdg); trigbarbin = hPDG.FindBin(trigbar.pdg)
        N_trig = hPDG.GetBinContent(trigbin) + hPDG.GetBinContent(trigbarbin)
        N_trig *= factor

        x = outputfile.mkdir(f"{trig.name}", f"strangeness correlations with {trig.name} trigger, background subtracted")
        x.cd()
        hratio = TH1D("hratio", "title", ndims, 0, ndims)
        hratio.SetStats(0)
        hratio.SetOption("HIST")
        hratio.SetFillColor(38)
        hratio.SetMarkerStyle(0)
        hratio.SetLineColor(1)
        hratio.SetLineStyle(1)
        
        strangeness_sum = 0.

        # K0_S/L together so we can use K+ for bkg subtraction
        histo = inputfile.Get(f"{trig.name}/{Kzeroshort.name}").Clone() 
        histo.Add(inputfile.Get(f"{trig.name}/{Kzerolong.name}"))
        histo = histo.ProjectionX()
        bkghisto = inputfile.Get(f"{trig.name}ss/{Kminus.name}")
        bkghisto = bkghisto.ProjectionX()
        bkgnormhisto = inputfile.Get(f"{trig.name}/{Kplus.name}")
        bkgnormhisto = bkgnormhisto.ProjectionX()
        bkgscaling = 0.5*float(histo.Integral())/float(bkgnormhisto.Integral())
        # bkgscaling = 1.
        print(bkgscaling)
        histo.Add(bkghisto, -2*bkgscaling) # TODO: check that the weights/errors are properly propagated
        histo.SetName("K0SL_bkgsub")
        histo.SetDirectory(outputfile)
        histo.Write()
        N_assoc = histo.Integral()/2 # divide by 2 because we have strange + antistrange here (i.e. K0 short + long)
        print("K0sl", N_assoc)
        hratio.Fill("K0SL_bkgsub", N_assoc)
        strangeness_sum += N_assoc
        
        # loop for other associated hadrons
        for j,assoc in enumerate(antistrangelist):
            assocbar = strangelist[j]
            histo = inputfile.Get(f"{trig.name}/{assoc.name}").Clone()
            histo.Add(inputfile.Get(f"{trig.name}ss/{assocbar.name}"), -1) # TODO: check that the weights/errors are properly propagated
            histo.SetName(f"{assoc.name}_bkgsub")
            histo.SetDirectory(outputfile)
            histo.Write()
            N_assoc = histo.Integral()
            # test = histo.Integral()
            print(assoc.name, N_assoc)
            hratio.Fill(assoc.name, N_assoc)
            if assoc == Xizerobar or assoc == Xiplus: # double strange hadrons obviously count double
                strangeness_sum += 2*N_assoc
            elif assoc == Omegaplus:
                strangeness_sum += 3*N_assoc
            else:
                strangeness_sum += N_assoc
        
        hratio.Write()
        print(f"total strangeness associated with {trig.name} trigger is {strangeness_sum/N_trig}, found with {N_trig} triggers")

    inputfile.Close()
    outputfile.Close()
    print(factor)

    print("Finished analysis()")
    return


antistrangelist = [Kplus, Lambdabar, Sigmaminusbar, Sigmaplusbar, Sigmazerobar, Xizerobar, Xiplus, Omegaplus]
strangelist = [Kminus, Lambda, Sigmaminus, Sigmaplus, Sigmazero, Xizero, Ximinus, Omegaminus]

# make_projections("output/ssbarv4_500M_14TeV_Monash.root", parallel = True)

analysis("output/Projections.root")

# inputfile = TFile("output/ssbarv4_500M_14TeV_Monash.root", "READ") # READ should be default
# hTHn = inputfile.Get("hSS")
# inputfile.Close()
# outputfile = TFile("output/AnalysisResults.root", "READ")
# test = project(hTHn, 4, trigger = Ximinus)
# test.Draw()
