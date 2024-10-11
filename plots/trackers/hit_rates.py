#!/usr/bin/python3

# hit rates for individual pixels

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import TH2D

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 2

    func = {}
    func[0] = rate_xy
    func[1] = rate_1d
    func[2] = rate_s1_s2

    func[iplot]()

#main

#_____________________________________________________________________________
def rate_xy(show_plot=True):

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9a/hits_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax1/hits_s2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax2/hits_s2.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ax2/hits.root"

    #number of simulated events
    #nev = 1e6
    #nev = 9980000 # 9ax2
    nev = 12000000

    #bunch crossing frequency
    freq = 22676 # kHz

    infile = TFile.Open(inp)
    #hxy = infile.Get("h_counts")
    hxy = infile.Get("lowQ2_s2_4_h_counts")

    scale = freq/nev
    print("scale:", scale)

    hxy.Scale(scale)

    print("min:", hxy.GetBinContent(hxy.GetMinimumBin()))
    print("max:", hxy.GetBinContent(hxy.GetMaximumBin()))

    if not show_plot:
        return infile, hxy

    can = ut.box_canvas(3000, 2000)

    hxy.Draw("colz")

    hxy.SetMinimum(scale*0.98)
    hxy.SetMaximum(1.1*hxy.GetBinContent(hxy.GetMaximumBin()))

    hxy.SetContour(300)

    hxy.SetTitle("")

    ut.put_yx_tit(hxy, "Pixel number in #it{y}", "Pixel number in #it{x}", 1.6, 1.4)

    ut.set_margin_lbtr(gPad, 0.11, 0.12, 0.03, 0.14)

    hxy.SetZTitle("Hit rate per pixel (kHz)")
    hxy.SetTitleOffset(1.2, "Z")

    gPad.SetLogz()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.png")

#rate_xy

#_____________________________________________________________________________
def rate_1d():

    xbin = 0.04
    #xmax = 1.5
    #xbin = 0.07
    xmax = 2.5

    hx = ut.prepare_TH1D("hx", xbin, 0, xmax)

    infile, hxy = rate_xy(False)

    for ix in range(1, hxy.GetNbinsX()+1):
        for iy in range(1, hxy.GetNbinsY()+1):
            hx.Fill( hxy.GetBinContent(ix, iy) )

    can = ut.box_canvas()

    ut.set_H1D_col(hx, rt.kRed)

    hx.Draw()

    ut.put_yx_tit(hx, "Number of pixels", "Hit rate per pixel (kHz)", 1.9, 1.4)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.02, 0.02)

    gPad.SetLogy()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#rate_1d

#_____________________________________________________________________________
def rate_1d_inp(inp, nev, freq, hx, plane_nam="h_counts"):

    infile = TFile.Open(inp)
    hxy = infile.Get(plane_nam)

    scale = freq/nev

    hxy.Scale(scale)

    for ix in range(1, hxy.GetNbinsX()+1):
        for iy in range(1, hxy.GetNbinsY()+1):
            hx.Fill( hxy.GetBinContent(ix, iy) )

#rate_1d_inp

#_____________________________________________________________________________
def rate_s1_s2():

    #inp_s1 = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax1/hits_s1.root"
    #inp_s2 = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax1/hits_s2.root"
    #inp_s1 = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax2/hits_s1.root"
    #inp_s2 = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax2/hits_s2.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ax2/hits.root"

    #xbin = 0.07
    #xmax = 2.5
    xbin = 0.04
    xmax = 1.4

    #number of simulated events
    #nev = 1e6
    #nev = 9980000 # tag9ax2
    nev = 12000000

    #bunch crossing frequency
    freq = 22676 # kHz

    hs1 = ut.prepare_TH1D("hs1", xbin, 0, xmax)
    hs2 = ut.prepare_TH1D("hs2", xbin, 0, xmax)

    #rate_1d_inp(inp_s1, nev, freq, hs1)
    #rate_1d_inp(inp_s2, nev, freq, hs2)
    rate_1d_inp(inp, nev, freq, hs1, "lowQ2_s1_4_h_counts")
    rate_1d_inp(inp, nev, freq, hs2, "lowQ2_s2_4_h_counts")

    can = ut.box_canvas()

    ut.set_H1D_col(hs1, rt.kBlue)
    ut.set_H1D_col(hs2, rt.kRed)

    hs1.Draw()
    hs2.Draw("e1same")

    ut.put_yx_tit(hs1, "Number of pixels", "Hit rate at single pixel (kHz)", 1.9, 1.4)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.02, 0.03)

    leg = ut.prepare_leg(0.7, 0.81, 0.24, 0.1, 0.035) # x, y, dx, dy, tsiz
    leg.AddEntry(hs1, "Tagger 1", "lp")
    leg.AddEntry(hs2, "Tagger 2", "lp")
    leg.Draw("same")

    gPad.SetLogy()

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#rate_s1_s2

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()












