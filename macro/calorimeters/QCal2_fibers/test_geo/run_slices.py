#!/usr/bin/env python3

from ctypes import c_double

import ROOT as rt
from ROOT import gROOT, gStyle, gPad, TGraph, TMath

import sys
sys.path.append("../../../../plots/")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    gROOT.ProcessLine(".include ../../../../calorimeters/include")
    gROOT.ProcessLine(".L ../../../../calorimeters/src/FiberYZ.cxx")

    fib = rt.FiberYZ(10, 8, 1.1, 0.9)

    #print_facets(fib)
    return

    can = ut.box_canvas()

    #frame = can.DrawFrame(-1, -0.1, 11, 1.1)
    frame = can.DrawFrame(-1.6, -1.6, 11, 11)
    frame.Draw()

    axis = make_axis(fib)
    axis.Draw("lsame")

    points = make_slice_points(fib)
    points.Draw("psame")

    lin = make_tangent_all(fib)
    for i in lin: i.Draw("lsame")

    #lin = make_tangent_form(fib, fib.GetSlice(7).z, 1)
    #lin.SetLineColor(rt.kCyan)
    #lin.SetLineStyle(rt.kDashed)
    #lin.Draw("lsame")

    gPad.SetGrid()

    ut.invert_col(gPad)
    can.SaveAs("01fig.pdf")

#_____________________________________________________________________________
def make_axis(fib):

    g = TGraph(fib.GetNSlice())

    for i in range(g.GetN()):
        slc = fib.GetSlice(i)
        #print(i, z)
        g.SetPoint(i, slc.z, slc.y)

    g.SetLineColor(rt.kRed)

    return g

#make_axis

#_____________________________________________________________________________
def make_slice_points(fib):

    #g = TGraph(fib.GetNSlice()*2)
    g = TGraph(fib.GetNSlice())
    ip = 0;

    for i in range(fib.GetNSlice()):
        slc = fib.GetSlice(i)

        g.SetPoint(ip, slc.points.at(0).at(2), slc.points.at(0).at(1))
        ip += 1
        #g.SetPoint(ip, slc.points.at(1).at(2), slc.points.at(1).at(1))
        #ip += 1

    g.SetLineColor(rt.kViolet)
    g.SetMarkerColor(rt.kViolet)
    g.SetMarkerStyle(rt.kStar)

    return g

#make_slice_points

#_____________________________________________________________________________
def make_tangent_form(fib, z1, r):

    #tangent line by direct formula
    a1 = c_double(0)
    b1 = c_double(0)
    fib.MakeNormVect(z1, a1, b1)

    a1 = a1.value
    b1 = b1.value

    npoints = 100
    tmin = -r
    dt = 2*r/(npoints-1)

    y1 = fib.f(z1)

    g = TGraph(npoints)

    for i in range(npoints):
        t = tmin + i*dt
        z = z1 + a1*t
        y = y1 + b1*t

        g.SetPoint(i, z, y)

    g.SetLineColor(rt.kCyan)

    return g

#make_tangent_form

#_____________________________________________________________________________
def make_tangent_all(fib):

    lin = []

    for i in range(fib.GetNSlice()):

        lx = make_tangent_form(fib, fib.GetSlice(i).z, fib.GetR())
        lx.SetLineColor(rt.kCyan)
    #lin.SetLineStyle(rt.kDashed)
    #lx.Draw("lsame")

        lin.append(lx)

    return lin

#make_tangent_all

#_____________________________________________________________________________
def print_facets(fib):

    for i in range(fib.GetNFacets()):

        fct = fib.GetFacet(i);
        print("Facet", i)
        print("p0:", fct.p0.at(0), fct.p0.at(1), fct.p0.at(2))
        print("p1:", fct.p1.at(0), fct.p1.at(1), fct.p1.at(2))
        print("p2:", fct.p2.at(0), fct.p2.at(1), fct.p2.at(2))

#print_facets

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()














