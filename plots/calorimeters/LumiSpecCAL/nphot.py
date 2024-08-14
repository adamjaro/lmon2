#!/usr/bin/python3

# number of optical photons in SiPMs

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import RDataFrame

import sys
sys.path.append("../../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    inp = "/home/jaroslav/sim/lmon2/macro/calorimeters/LumiSpecCAL/cal_opt_0p5GeV_8evt.root"

    xbin = 10
    xmin = 0
    xmax = 150000

    df = RDataFrame("DetectorTree", inp)

    #sum of optical photons from individual SiPMs
    df = df.Define("nphot_all", "return LumiSpecCAL_opdet_time.size();")

    can = ut.box_canvas()
    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )
    hx = df.Histo1D(hx, "nphot_all").GetValue()
    ut.set_H1D(hx)

    hx.Draw("e1")

    gPad.SetGrid()

    #gPad.SetLogy()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#main


#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetFrameLineWidth(2)

    main()

























