#!/usr/bin/python3

# beam-gas properties

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import RDataFrame

import sys
sys.path.append("../../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 0

    func = {}
    func[0] = theta_on_z

    func[iplot]()

#main

#_____________________________________________________________________________
def theta_on_z():

    pass

#theta_on_z

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetFrameLineWidth(2)

    main()







