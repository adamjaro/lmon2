#!/usr/bin/env python3

from ctypes import c_double
import code

import ROOT as rt
from ROOT import gROOT, gStyle, gPad, TGraph2D

import sys
sys.path.append("../../../../plots/")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    gROOT.ProcessLine(".include ../../../../calorimeters/include")
    gROOT.ProcessLine(".L ../../../../calorimeters/src/FiberYZ.cxx")

    fib = rt.FiberYZ(10, 8, 1.1)

    #return

    #can = ut.box_canvas()

    g = TGraph2D()

    for i in range(fib.GetNSlice()):
        slc = fib.GetSlice(i)

        for ip in range(slc.points.size()):

            g.AddPoint(slc.points.at(ip).at(2), slc.points.at(ip).at(1), slc.points.at(ip).at(0))



    g.SetLineColor(rt.kViolet)
    g.SetMarkerColor(rt.kViolet)
    g.SetMarkerStyle(rt.kStar)

    g.Draw("ap")



    variables = globals().copy()
    variables.update(locals())
    shell = code.InteractiveConsole(variables)
    shell.interact()

#main


#_____________________________________________________________________________
if __name__ == "__main__":

    #gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()



