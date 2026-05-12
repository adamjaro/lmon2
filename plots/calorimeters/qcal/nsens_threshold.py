#!/usr/bin/env python3

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, RDataFrame

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

    inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal2bx1/en_","/lmon.root"]

    #energy = [1, 5, 9, 14, 18]
    #energy = [5, 9]
    energy = [1, 9, 18]

    hx = {}
    for i in energy:
        hx[i] = run_single(inp[0]+str(i)+inp[1])

    hx[1].SetLineColor(rt.kViolet)
    #hx[5].SetLineColor(rt.kCyan+1)
    hx[9].SetLineColor(rt.kRed)
    #hx[14].SetLineColor(rt.kBlue)
    hx[18].SetLineColor(rt.kGreen+1)

    #hx[5].SetLineStyle(rt.kDashed)
    #hx[18].SetLineStyle(rt.kDashed)

    can = ut.box_canvas()

    #h1 = hx[energy[len(energy)-1]]
    h1 = hx[energy[0]]
    ut.put_yx_tit(h1, "Counts", "Sensors with at least 10 incident photons in event", 1.8, 1.3)

    h1.Draw()

    #for i in energy[:-1]:
    for i in energy[1:]:
        hx[i].Draw("same")

    ut.set_margin_lbtr(gPad, 0.12, 0.1, 0.02, 0.02)

    leg = ut.prepare_leg(0.67, 0.72, 0.25, 0.2, 0.04)
    leg.SetHeader("Generated #it{E}_{#gamma} =")
    for i in energy:
        leg.AddEntry(hx[i], str(i)+" GeV", "l")

    leg.Draw("same")

    gPad.SetGrid()

    #gPad.SetLogy()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#run_all

#_____________________________________________________________________________
def run_single(inp=None):

    draw = False
    if inp is None:
        inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal2bx1/en_9/lmon.root"
        draw = True

    xbin = 1
    xmin = 0
    xmax = 19

    df = RDataFrame("DetectorTree", inp)

    #number of sensors above a threshold in event
    df = df.Define("qcal_nsens_th0", "map<Int_t, Int_t> sens; for(Int_t i: qcal_opdet_cell_id) {\
        if(sens.find(i) == sens.end()) {sens.emplace(i, 1);} else {sens[i]++;}}\
        Int_t cnt=0; for(pair<Int_t, Int_t> p: sens) {if(p.second >= 10) {cnt++;}} return cnt;")

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )
    hx = df.Histo1D(hx, "qcal_nsens_th0").GetValue()

    ut.line_h1(hx, rt.kBlue, 2)

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

    rt.EnableImplicitMT()

    main()

