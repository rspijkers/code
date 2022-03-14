#!/usr/bin/env python3

# NEVER IMPORT ROOT BEFORE NUMPY!!!
# incompatibilities may lead to segmentation violation
from hadrons import * # import hadron classes to easily access well defined properties such as pdg codes.
import numpy as np
from ROOT import TFile, TDirectory
import os#, sys

# function that loads THn from given inputfile
def load_THn(filename, THn_name):
    # TODO: add function description
    """
    """
    inputfile = TFile(filename, "READ") # READ should be default
    if inputfile.IsZombie():
        print(f"Error: Failed to open file {inputfile}")
        raise Exception('IsZombie')
    THn = inputfile.Get(THn_name)
    inputfile.Close()
    return THn
    # maybe implement that the function searches for THn in file, no need to specify THn name

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

    return projection

# main function, called at the end of the script
def main():
    inputfilepath = "output/ssbarv2_100M_14TeV_Monash.root"
    hTHn = load_THn(inputfilepath, "hSS")
    hNtrigger = TFile(inputfilepath, "READ").Get("hPDG")

    axes = produce_axes_dict(hTHn)

    # NOTE: histo's created after this point are automatically added to outputfile
    outputfile = TFile("output/PyRootOutput.root", "RECREATE")

    # TODO: set global plot settings? smt like set_plotstyle()

    # create all histograms here

    # create directories in root output file
    Xi0dir = outputfile.mkdir('Xi0', 'strangeness correlations with Xi0 trigger');
    Xi0ssdir = outputfile.mkdir('Xi0ss', 'same-sign strangeness correlations with Xi0 trigger');
    Ximinusdir = outputfile.mkdir('Xi-', 'strangeness correlations with Ximinus trigger');
    Ximinusssdir = outputfile.mkdir('Xi-ss', 'same-sign strangeness correlations with Ximinus trigger');
    Lambdadir = outputfile.mkdir('Lambda', 'strangeness correlations with Lambda trigger');
    Lambdassdir = outputfile.mkdir('Lambdass', 'same-sign strangeness correlations with Lambda trigger');
    K0dir = outputfile.mkdir('K0', 'strangeness correlations with K0 trigger');
    K0ssdir = outputfile.mkdir('K0ss', 'same-sign strangeness correlations with K0 trigger');
    Kminusdir = outputfile.mkdir('K-', 'strangeness correlations with Kminus trigger');
    Kminusssdir = outputfile.mkdir('K-ss', 'same-sign strangeness correlations with Kminus trigger');

    # make the 2D pT,eta histogram for trigger
    histo = project(hTHn, axes['pTTrigger'], axes['etaTrigger'], trigger = Xizero, name = "h2_Xi0_Trigger_pT_eta")
    histo.Add(project(hTHn, axes['pTTrigger'], axes['etaTrigger'], trigger = Xizerobar))
    histo.Write()
    histo = project(hTHn, axes['pTTrigger'], axes['etaTrigger'], trigger = Ximinus, name = "h2_Xi-_Trigger_pT_eta")
    histo.Add(project(hTHn, axes['pTTrigger'], axes['etaTrigger'], trigger = Xiplus))
    histo.Write()

    # make the deltaPhi correlation histo's
    assoclist = [Kplus, Kzero, Lambdabar, Xizerobar, Xiplus]
    antiassoclist = [Kminus, Kzerobar, Lambda, Xizero, Ximinus]
    for i in range(5):
        # Xi0
        Xi0dir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Xizero, associated = assoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Xizerobar, associated = antiassoclist[i]))
        histo.Write()
        # Xi0 samesign
        Xi0ssdir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Xizero, associated = antiassoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Xizerobar, associated = assoclist[i]))
        histo.Write()
        # Ximinus
        Ximinusdir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Ximinus, associated = assoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Xiplus, associated = antiassoclist[i]))
        histo.Write()
        # Ximinus
        Ximinusssdir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Ximinus, associated = antiassoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Xiplus, associated = assoclist[i]))
        histo.Write()

        # Lambda
        Lambdadir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Lambda, associated = assoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Lambdabar, associated = antiassoclist[i]))
        histo.Write()
        # Lambda samesign
        Lambdassdir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Lambda, associated = antiassoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Lambdabar, associated = assoclist[i]))
        histo.Write()

        # TODO: Sigma

        # K0
        K0dir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Kzero, associated = assoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Kzerobar, associated = antiassoclist[i]))
        histo.Write()
        # K0 samesign
        K0ssdir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Kzero, associated = antiassoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Kzerobar, associated = assoclist[i]))
        histo.Write()
        # Kminus
        Kminusdir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Kminus, associated = assoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Kplus, associated = antiassoclist[i]))
        histo.Write()
        # Kminus
        Kminusssdir.cd()
        histo = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Kminus, associated = antiassoclist[i])
        histo.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Kplus, associated = assoclist[i]))
        histo.Write()

    # try if we can access a histo we just wrote to file:
    # test = outputfile.Get("Xi0/h2Xi0bar_deltaPhi_pTAssoc").ProjectionX()
    # for i in range(4):
    #     print(f"mean of axis {i} is {test.GetMean(i)}")

    outputfile.cd()    # return to parent

    # quick test with lambda-lambda(bar) correlation
    hll = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Lambda, associated = Lambda)
    hll.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Lambdabar, associated = Lambdabar))
    hll.Write()
    hllbar = project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Lambda, associated = Lambdabar)
    hllbar.Add(project(hTHn, axes['deltaPhi'], axes['pTAssoc'], trigger = Lambdabar, associated = Lambda))
    hllbar.Write()

    # outputfile.Write()
    outputfile.Close()

main()
