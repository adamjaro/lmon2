#!/usr/bin/python3

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TGaxis
from ROOT import TFile

import sys
sys.path.append("../../../plots")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    #gStyle.SetPalette(1)

    inp = TFile("gauss_solve.root", "read")

    can = ut.box_canvas()

    g2 = inp.Get("g2")

    #print(type(g2))

    g2.SetNameTitle("", "")

    g2.GetXaxis().SetTitle("dddddd")
    #g2.GetYaxis().SetTitle("dddddd")
    #g2.GetXaxis().SetTitleOffset(1)

    g2.GetXaxis().CenterTitle()

    g2.GetXaxis().SetLabelColor(rt.kRed)

    g2.GetYaxis().SetTitle("fff")

    g2.SetLineColor(rt.kRed)

    g2.Draw("tri") # surf1 atri triw

    #g2.Draw("p") # axes titles seen with p
    #g2.Draw("triwsame")
    #g2.Draw("p")

    #g2.GetXaxis().Draw()

    #gPad.SetFrameFillColor(rt.kCyan+2)
    gPad.SetFrameFillColor(rt.kGreen+2)
    #gPad.SetFrameFillColor(rt.kBlue)
    #gPad.SetLineColor(rt.kBlue)

    #gPad.SetFrameLineColor(rt.kGreen+2)

    #gPad.GetFrame().SetLineColor(rt.kBlue)
    #gPad.GetFrame().SetFillColor(rt.kYellow)

    #axisE = TGaxis(-3, -3, 0, -3, 3, 0)

    #axisE.SetTitle("ggg")

    #axisE.Draw("same")

    #gPad.SetGrid(0, 0)
    gPad.SetFillColor(rt.kBlue)

    #gPad.Update()

    #g2.GetZaxis().Draw("same")

    #gStyle.SetAxisColor(rt.kYellow, "X")

    gPad.SetTopMargin(0.03)
    gPad.SetRightMargin(0.03)
    gPad.SetBottomMargin(0.05)
    gPad.SetLeftMargin(0.08)

    #gStyle.SetPadColor(rt.kGreen)

    ut.print_pad(gPad)

    ut.invert_col(gPad)
    can.SaveAs("01fig.pdf")

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetFrameLineWidth(2)

    main()

