#!/usr/bin/env python3

#PDE as a function of wavelength

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, gSystem, gInterpreter, RDataFrame, RDF
from ROOT import TGraph

import sys
sys.path.append("../..")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    df = RDF.FromCSV("../../../../lmon2-data/pde_wavelength_30035.csv")

    #df.Display().Print()

    df = df.Define("PDE_5V", "0.01*DetectionEfficiency_5V_percent")
    df = df.Define("PDE_2_5V", "0.01*DetectionEfficiency_2_5V_percent")

    #g1 = df.Graph("Wavelength_nm", "DetectionEfficiency_5V_percent")
    #g2 = df.Graph("Wavelength_nm", "DetectionEfficiency_2_5V_percent")
    g1 = df.Graph("Wavelength_nm", "PDE_5V")
    g2 = df.Graph("Wavelength_nm", "PDE_2_5V")

    xmin = 300
    xmax = 950

    npx = 2000
    dx = (xmax-xmin)/npx

    g1n = TGraph()
    g2n = TGraph()

    for i in range(npx):
        wi = xmin + i*dx
        g1n.AddPoint(wi, g1.Eval(wi))
        g2n.AddPoint(wi, g2.Eval(wi))


    ut.set_graph(g1, rt.kRed)
    ut.set_graph(g2, rt.kRed)
    #g2.SetLineStyle(rt.kDashed)

    ut.set_graph(g1n, rt.kGreen)
    ut.set_graph(g2n, rt.kGreen)
    g1n.SetLineStyle(rt.kDashed)
    g2n.SetLineStyle(rt.kDashed)

    can = ut.box_canvas(1024, 768)

    #hx = gPad.DrawFrame(300, 0, 950, 45)
    hx = gPad.DrawFrame(xmin, 0, xmax, 0.45)

    g1.Draw("l")
    g2.Draw("lsame")
    g1n.Draw("lsame")
    g2n.Draw("lsame")

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#main

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    #rt.EnableImplicitMT()

    main()








