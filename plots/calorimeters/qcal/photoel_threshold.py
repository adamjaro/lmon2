#!/usr/bin/env python3

#number of sensors having photoelectron count in event above a threshold

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, gSystem, gInterpreter
from ROOT import TGraph, RDataFrame, RDF

from photoelectrons import photoelectrons

import sys
sys.path.append("../..")
import plot_utils as ut


#_____________________________________________________________________________
def main():

    iplot = 1

    func = {}
    func[0] = run_single
    func[1] = run_all

    func[iplot]()

#main

#_____________________________________________________________________________
def run_all():

    #inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal3cx2/en_","/lmon.root"]
    inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal3cx3/en_","/lmon.root"]

    #energy = [1, 5, 9, 14, 18]
    energy = [1, 9, 18]

    #threshold in fired microcells
    thres = 10

    hx = {}
    for i in energy:
        hx[i] = run_single(inp[0]+str(i)+inp[1], thres)

    hx[1].SetLineColor(rt.kViolet)
    #hx[5].SetLineColor(rt.kCyan+1)
    hx[9].SetLineColor(rt.kRed)
    #hx[14].SetLineColor(rt.kBlue)
    hx[18].SetLineColor(rt.kGreen+1)

    can = ut.box_canvas()

    h1 = hx[1]
    ut.put_yx_tit(h1, "Counts / event", "Sensors with at least "+str(thres)+" fired microcells in event", 1.6, 1.3)

    h1.Draw()

    for i in energy[1:]:
        hx[i].Draw("same")

    ut.set_margin_lbtr(gPad, 0.11, 0.1, 0.02, 0.02)

    leg = ut.prepare_leg(0.67, 0.7, 0.25, 0.25, 0.04)
    leg.SetHeader("Generated #it{E}_{#gamma}:")
    for i in energy:
        leg.AddEntry(hx[i], str(i)+" GeV", "l")

    leg.Draw("same")

    gPad.SetGrid()

    #gPad.SetLogy()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#_____________________________________________________________________________
def run_single(inp=None, thres=10):

    draw = False
    if inp is None:
        inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal3cx2/en_9/lmon.root"
        draw = True

    xbin = 1
    xmin = 0
    xmax = 24

    df = RDataFrame("DetectorTree", inp)

    #detected flag
    df = df.Define("qcal_opdet_is_detected", "sens(qcal_opdet_phot_en)")

    #number of photoelectrons in sensor in event
    df = df.Define("qcal_nphotoel_sens", "nphotoel_sens(qcal_opdet_cell_id, qcal_opdet_is_detected)")

    #number of sensors above photoelectron threshold
    df = df.Define("qcal_nsens_photoel", "nsens_photoel(qcal_nphotoel_sens, "+str(thres)+")")

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )
    hx = df.Histo1D(hx, "qcal_nsens_photoel").GetValue()

    #number of events processed
    nev = df.Count().GetValue()

    ut.line_h1(hx, rt.kBlue, 3)

    hx.Scale(1./nev)
    print(hx.Integral())

    if not draw: return hx

    can = ut.box_canvas()

    hx.Draw("")

    gPad.SetGrid()

    gPad.SetLogy()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#run_single

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    #rt.EnableImplicitMT()

    #photoelectron functions
    photoelectrons()

    #call to main
    main()

































