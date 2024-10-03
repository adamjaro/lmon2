#!/usr/bin/python3

# generator true kinematic quantities

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import TH2D

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

  lQ2()

#main

#_____________________________________________________________________________
def lQ2():

    #GeV^2
    qbin = 0.1
    qmin = -8
    qmax = -1

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ax1/acc.root"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    hxy = ut.prepare_TH1D("hx", qbin, qmin, qmax)
    tree.Draw("TMath::Log10(true_Q2)", "s1_allhit==1")

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#lQ2

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()






























