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

    iplot = 0

    func = {}
    func[0] = rate_xy
    func[1] = rate_1d

    func[iplot]()

#main

#_____________________________________________________________________________
def rate_xy(show_plot=True):

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9a/hits_v2.root"

    #number of simulated events
    nev = 1e6

    #bunch crossing frequency
    freq = 22676 # kHz

    infile = TFile.Open(inp)
    hxy = infile.Get("h_counts")

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

    hxy.SetZTitle("Hit rate (kHz)")
    hxy.SetTitleOffset(1.2, "Z")

    gPad.SetLogz()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.png")

#rate_xy

#_____________________________________________________________________________
def rate_1d():

    hx = ut.prepare_TH1D("hx", 0.05, 0, 1.5)

    infile, hxy = rate_xy(False)

    for ix in range(1, hxy.GetNbinsX()+1):
        for iy in range(1, hxy.GetNbinsY()+1):
            hx.Fill( hxy.GetBinContent(ix, iy) )

    can = ut.box_canvas()

    ut.set_H1D_col(hx, rt.kRed)

    hx.Draw()

    ut.put_yx_tit(hx, "Number of pixels", "Hit rate (kHz)", 1.9, 1.4)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.02, 0.02)

    gPad.SetLogy()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#rate_1d

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()

