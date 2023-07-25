#!/usr/bin/python3

# track basic parameters

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 3

    func = {}
    func[0] = theta_x
    func[1] = theta_y
    func[2] = xy
    func[3] = chi2

    func[iplot]()

#main

#_____________________________________________________________________________
def theta_x():

    #mrad
    xbin = 1
    xmin = -10
    xmax = 40

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6a/trk_v1.root"

    #det = "s1_tracks"
    det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw(det+"_theta_x*1e3 >> hx")

    ut.put_yx_tit(hx, "Counts", "theta_x (mrad)", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#theta_x

#_____________________________________________________________________________
def theta_y():

    #mrad
    xbin = 1
    xmin = -10
    xmax = 40

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6a/trk_v1.root"

    #det = "s1_tracks"
    det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw(det+"_theta_y*1e3 >> hx")

    ut.put_yx_tit(hx, "Counts", "theta_y (mrad)", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#theta_y

#_____________________________________________________________________________
def xy():

    #track position in x and y

    #mrad
    xbin = 1
    xmin = -80
    xmax = 80

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6a/trk_v1.root"

    #det = "s1_tracks"
    det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, xbin, xmin, xmax)

    tree.Draw(det+"_y:"+det+"_x >> hxy")

    ut.put_yx_tit(hxy, "y (mm)", "x (mm)", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#xy

#_____________________________________________________________________________
def chi2():

    xbin = 1e-3
    xmin = 0
    xmax = 0.06

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6a/trk_v1.root"

    #det = "s1_tracks"
    det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw(det+"_chi2_xy >> hx")

    ut.put_yx_tit(hx, "Counts", "chi2", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#chi2

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()








