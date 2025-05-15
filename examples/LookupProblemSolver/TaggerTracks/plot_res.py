#!/usr/bin/python3

import ROOT as rt
from ROOT import gPad, gROOT, gStyle
from ROOT import TFile

import sys
sys.path.append("../../../plots")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    infile = TFile.Open("rec_tracks.root")
    tree = infile.Get("rec_trk")

    xbin = 0.5
    xmin = 5
    xmax = 19

    can = ut.box_canvas()

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, xbin, xmin, xmax)

    tree.Draw("rec_en:true_en >> hxy")

    ytit = "Reconstructed energy #it{E_{e}} (GeV)"
    #xtit = "Generated MC particle energy #it{E_{e,mc}} (GeV)"
    xtit = "Generated true electron energy #it{E_{e,true}} (GeV)"
    ut.put_yx_tit(hxy, ytit, xtit, 1.5, 1.4)

    ut.set_margin_lbtr(gPad, 0.12, 0.11, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetLogz()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetFrameLineWidth(2)

    main()

