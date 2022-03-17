#!/usr/bin/env python3

# NEVER IMPORT ROOT BEFORE NUMPY!!!
# incompatibilities may lead to segmentation violation
from hadrons import * # import hadron classes to easily access well defined properties such as pdg codes.
import numpy as np
from joblib import Parallel, delayed
from ROOT import TFile, TH1
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

    range_user = kwargs.get('range')
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

    inputfile = TFile(inputfilepath, "READ") # READ should be default
    if inputfile.IsZombie():
        print(f"Error: Failed to open file {inputfile}")
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
    hPDG.SetDirectory(outputfile)
    hEtaPt.SetDirectory(outputfile)
    hPDG.Write()
    hEtaPt.Write()
    inputfile.Close()

    # # make the 2D pT,eta histogram for trigger
    # histo = project(hTHn, axes['pTTrigger'], axes['etaTrigger'], trigger = Xizero, name = "h2_Xi0_Trigger_pT_eta")
    # histo.Add(project(hTHn, axes['pTTrigger'], axes['etaTrigger'], trigger = Xizerobar))
    # histo.Write()
    # histo = project(hTHn, axes['pTTrigger'], axes['etaTrigger'], trigger = Ximinus, name = "h2_Xi-_Trigger_pT_eta")
    # histo.Add(project(hTHn, axes['pTTrigger'], axes['etaTrigger'], trigger = Xiplus))
    # histo.Write()

    # make the deltaPhi correlation histo's
    if parallel:
        def write_corr(trig, i, assoc, j):
            histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = assoc, name = assoc.name)
            histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = antistrangelist[i], associated = strangelist[j]))
            return histo.Write()
            # antihisto.Delete()
        def write_ss_corr(trig, i, assoc, j):
            histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = assoc, name = assoc.name)
            histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = antistrangelist[i], associated = antistrangelist[j]))
            return histo.Write()
            # antihisto.Delete()

        for i,trig in enumerate(strangelist):
            # opposite sign correlations
            x = outputfile.mkdir(trig.name, f"strangeness correlations with {trig.name} trigger")
            x.cd()
            Parallel(n_jobs = len(strangelist), require='sharedmem')(delayed(write_corr)(trig, i, assoc, j) for j,assoc in enumerate(antistrangelist))
            
            # same sign correlations
            x = outputfile.mkdir(f"{trig.name}ss", f"samesign strangeness correlations with {trig.name} trigger")
            x.cd()
            Parallel(n_jobs = len(strangelist), require='sharedmem')(delayed(write_ss_corr)(trig, i, assoc, j) for j,assoc in enumerate(strangelist))
    elif not parallel:
        for i,trig in enumerate(strangelist):
            # opposite sign correlations
            x = outputfile.mkdir(trig.name, f"strangeness correlations with {trig.name} trigger")
            x.cd()
            for j,assoc in enumerate(antistrangelist): 
                histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = assoc, name = assoc.name)
                histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = antistrangelist[i], associated = strangelist[j]))
                histo.Write()
                
            # same sign correlations
            x = outputfile.mkdir(f"{trig.name}ss", f"samesign strangeness correlations with {trig.name} trigger")
            x.cd()
            for j,assoc in enumerate(strangelist): 
                histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = trig, associated = assoc, name = assoc.name)
                histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = antistrangelist[i], associated = antistrangelist[j]))
                histo.Write()

    # outputfile.cd()    # return to parent
    # outputfile.Write()
    outputfile.Close()
    return "Finished make_projections()"

def analysis(inputpath):
    inputfile = TFile(inputpath, "READ")
    TH1.AddDirectory(False);

    # primitively substract samesign correlations from opposite sign
    hPDG = inputfile.Get("hPDG")
    # assert that x-axis range of the pdg histo is symmetrical around 0 (same amount of bins on each side)
    nbins = hPDG.GetNbinsX()
    assert(-1 * int(hPDG.GetBinLowEdge(1)) == int(hPDG.GetBinLowEdge(nbins) + 1)), f"Error: PDG histogram range is not symmetrical around 0! Please do something about it (xmin = {int(hPDG.GetBinLowEdge(1))}, xmax = {int(hPDG.GetBinLowEdge(nbins) + 1)})"

    # get number of Xi0/Xi0bar
    Xi0bin = hPDG.FindBin(Xizero.pdg); Xi0barbin = hPDG.FindBin(Xizerobar.pdg)
    N_Xi0 = hPDG.GetBinContent(Xi0bin) + hPDG.GetBinContent(Xi0barbin)

    outputfile = TFile("output/AnalysisResults.root", "RECREATE")
    # do first only Xi0 trigger, later generalize
    strangeness_sum = 0.
    for i,assoc in enumerate(antistrangelist):
        histo = inputfile.Get(f"Xi0/{assoc.name}").Clone()
        histo.Add(inputfile.Get(f"Xi0ss/{strangelist[i].name}"), -1) # TODO: check that the weights/errors are properly propagated
        histo.SetName(f"Xi0{assoc.name}_bkgsub")
        histo.SetDirectory(outputfile)
        histo.Write()
        if assoc == Xizerobar or assoc == Xiplus: # double strange hadrons obviously count double
            strangeness_sum += 2*histo.GetEntries()
        elif assoc == Omegaplus:
            strangeness_sum += 3*histo.GetEntries()
        else:
            strangeness_sum += histo.GetEntries()
    
    print(strangeness_sum/N_Xi0)

    Ximinusbin = hPDG.FindBin(Ximinus.pdg); Xiplusbin = hPDG.FindBin(Xiplus.pdg)
    N_Xipm = hPDG.GetBinContent(Ximinusbin) + hPDG.GetBinContent(Xiplusbin)

    # do first only Xi- trigger, later generalize
    strangeness_sum = 0.
    for i,assoc in enumerate(antistrangelist):
        histo = inputfile.Get(f"Xi-/{assoc.name}").Clone()
        histo.Add(inputfile.Get(f"Xi-ss/{strangelist[i].name}"), -1) # TODO: check that the weights/errors are properly propagated
        histo.SetName(f"Xi-{assoc.name}_bkgsub")
        histo.SetDirectory(outputfile)
        histo.Write()
        if assoc == Xizerobar or assoc == Xiplus: # double strange hadrons obviously count double
            strangeness_sum += 2*histo.GetEntries()
        elif assoc == Omegaplus:
            strangeness_sum += 3*histo.GetEntries()
        else:
            strangeness_sum += histo.GetEntries()
    
    print(strangeness_sum/N_Xipm)

    inputfile.Close()
    outputfile.Close()
    return "Finished analysis()"

antistrangelist = [Kplus, Kzero, Lambdabar, Xizerobar, Xiplus, Sigmaminusbar, Sigmaplusbar, Sigmazerobar, Omegaplus]
strangelist = [Kminus, Kzerobar, Lambda, Xizero, Ximinus, Sigmaminus, Sigmaplus, Sigmazero, Omegaminus]

# make_projections("output/ssbarv3_500M_14TeV_Monash.root", parallel = True)

analysis("output/Projections.root")