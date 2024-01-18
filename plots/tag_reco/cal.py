#!/usr/bin/python3

# calorimeter basic parameters

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 0

    func = {}
    func[0] = xy

    func[iplot]()

#main

#_____________________________________________________________________________
def xy():

    #centroid position in x and y

    #mm
    xbin = 1
    xmin = -75
    xmax = 75

    #mm
    ybin = 0.5
    ymin = -30
    ymax = 30

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_v3.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_v3.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, ybin, ymin, ymax)

    tree.Draw(det+"_cal_y:"+det+"_cal_x >> hxy")

    ut.put_yx_tit(hxy, "#it{y} (mm)", "#it{x} (mm)", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetGrid()

    gPad.SetLogz()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#xy

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()





















