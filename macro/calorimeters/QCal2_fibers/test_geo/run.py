#!/usr/bin/env python3

import ROOT as rt
from ROOT import gROOT, gStyle, gPad, TGraph, TMath

import sys
sys.path.append("../../../../plots/")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    gROOT.ProcessLine(".include ../../../../calorimeters/include")
    gROOT.ProcessLine(".L ../../../../calorimeters/src/FiberYZ.cxx")

    fib = rt.FiberYZ(10, 8, 1)

    print(fib.f(0))
    print(fib.f(10))

    #can = ut.box_canvas(2000, 1000)
    can = ut.box_canvas()

    #frame = can.DrawFrame(-1, -0.1, 11, 1.1)
    frame = can.DrawFrame(-1, -1, 11, 11)
    frame.Draw()

    axis = make_axis(fib)
    axis.Draw("lsame")

    #lin = make_line(fib, 5)
    #lin.Draw("lsame")
    #lin.Draw("*same")

    z1 = 7

    lin = make_tangent_angle(fib, z1, 4)
    lin.Draw("lsame")

    lin2 = make_tangent_angle(fib, z1, 4, TMath.Pi()/2)
    lin2.SetLineColor(rt.kYellow)
    lin2.Draw("lsame")

    lin3 = make_tangent_form(fib, z1, 3)
    lin3.SetLineColor(rt.kRed)
    lin3.SetLineStyle(rt.kDashed)
    lin3.Draw("lsame")

    lin4 = make_tangent_form(fib, z1, 3, True)
    lin4.SetLineColor(rt.kViolet)
    lin4.SetLineStyle(rt.kDashed)
    lin4.Draw("lsame")

    gPad.SetGrid()

    ut.invert_col(gPad)
    can.SaveAs("01fig.pdf")

#_____________________________________________________________________________
def make_axis(fib):

    npoints = 100

    g = TGraph(npoints)

    dz = 10./(npoints-1)

    for i in range(npoints):
        z = i*dz
        y = fib.f(z)
        g.SetPoint(i, z, y)

    g.SetLineColor(rt.kRed)

    return g

#make_axis

#_____________________________________________________________________________
def make_line(fib, z1):

    #perpendicular line

    #g = TGraph(1)
    #g.SetPoint(0, z1, fib.f(z1))
    #g.SetMarkerColor(rt.kCyan)


    fP = fib.fP(z1)
    b1 = -TMath.Sqrt(1./(1.+fP*fP))
    a1 = -b1*fP
    y1 = fib.f(z1)

    print(fP, a1, b1, z1, y1)

    #return g

    g = TGraph(100)

    for i in range(100):
        t = i*0.1/99
        z = z1 + a1*t
        y = y1 + b1*t

        g.SetPoint(i, z, y)

    g.SetLineColor(rt.kCyan)

    return g

#_____________________________________________________________________________
def make_tangent_angle(fib, z1, r, theta_add=0):

    #tangent line with use of direction angle

    theta = TMath.ATan(fib.fP(z1)) + theta_add
    a1 = TMath.Cos(theta)
    b1 = TMath.Sin(theta)
    y1 = fib.f(z1)

    print(theta, a1, b1)

    npoints = 100
    tmin = -r
    dt = 2*r/(npoints-1)

    g = TGraph(npoints)

    for i in range(npoints):
        t = tmin + i*dt
        z = z1 + a1*t
        y = y1 + b1*t

        g.SetPoint(i, z, y)

    g.SetLineColor(rt.kCyan)

    return g

#make_tangent_angle

#_____________________________________________________________________________
def make_tangent_form(fib, z1, r, rot=False):

    #tangent line by direct formula
    fP = fib.fP(z1)
    a1 = TMath.Sqrt(1./(1+fP*fP))
    b1 = a1*fP
    y1 = fib.f(z1)

    if rot:
        a1old = a1
        a1 = b1
        b1 = -a1old

    print(a1, b1)

    npoints = 100
    tmin = -r
    dt = 2*r/(npoints-1)

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
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()




