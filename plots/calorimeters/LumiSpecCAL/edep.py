#!/usr/bin/python3

# deposited energy in fibers

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import RDataFrame

import sys
sys.path.append("../../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    #inp = "/home/jaroslav/sim/lmon2/macro/calorimeters/LumiSpecCAL/cal_1kevt.root"
    inp = "/home/jaroslav/sim/lmon2/macro/calorimeters/LumiSpecCAL/cal_10GeV_1kevt.root"

    #MeV
    xbin = 1
    xmin = 0
    #xmax = 60
    xmax = 600

    df = RDataFrame("DetectorTree", inp)

    #sum of deposited energy over all hits, conversion to MeV
    df = df.Define("edep_all", "Double_t e=0; for(auto& i: LumiSpecCAL_en) {e+=i;} return 1e3*e;")

    can = ut.box_canvas()
    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )
    hx = df.Histo1D(hx, "edep_all").GetValue()
    ut.set_H1D(hx)

    hx.Draw("e1")

    gPad.SetGrid()

    #gPad.SetLogy()

    print(hx.GetMean())
    print(hx.GetStdDev())
    print(hx.GetStdDev()/hx.GetMean())

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

























