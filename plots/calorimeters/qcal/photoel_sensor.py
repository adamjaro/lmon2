#!/usr/bin/env python3

#number of photoelectrons per sensor

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

    #inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal2cx3/en_","/lmon.root"]
    #inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal2cx4/en_","/lmon.root"]
    #inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal3cx1/en_","/lmon.root"]
    inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal3cx2/en_","/lmon.root"]

    energy = [1, 5, 9, 14, 18] # , 14

    hx = {}
    for i in energy:
        hx[i] = run_single(inp[0]+str(i)+inp[1])

    hx[1].SetLineColor(rt.kViolet)
    hx[5].SetLineColor(rt.kCyan+1)
    hx[9].SetLineColor(rt.kRed)
    hx[14].SetLineColor(rt.kBlue)
    hx[18].SetLineColor(rt.kGreen+1)

    can = ut.box_canvas()

    h1 = hx[energy[len(energy)-1]]
    ut.put_yx_tit(h1, "Counts", "Number of photoelectrons in sensor in event", 1.4, 1.3)

    h1.Draw()

    for i in energy[:-1]:
        hx[i].Draw("same")

    ut.set_margin_lbtr(gPad, 0.1, 0.1, 0.02, 0.04)

    leg = ut.prepare_leg(0.67, 0.66, 0.25, 0.3, 0.04)
    leg.SetHeader("Generated #it{E}_{#gamma}:")
    for i in energy:
        leg.AddEntry(hx[i], str(i)+" GeV", "l")

    leg.Draw("same")

    gPad.SetGrid()

    gPad.SetLogy()
    #gPad.SetLogx()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#_____________________________________________________________________________
def run_single(inp=None):

    draw = False
    if inp is None:
        #inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal2cx2/en_1/lmon.root"
        inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal2cx3/en_18/lmon.root"
        draw = True

    xbin = 8
    #xbin = 2
    xmin = 0
    #xmax = 100
    xmax = 650

    df = RDataFrame("DetectorTree", inp)
    #df = df.Range(12)

    #detected flag
    df = df.Define("qcal_opdet_is_detected", "sens(qcal_opdet_phot_en)")

    #number of photons in sensor in event
    df = df.Define("qcal_nphot_sens", "map<Int_t, Int_t> sens; for(Int_t i: qcal_opdet_cell_id) {\
        if(sens.find(i) == sens.end()) {sens.emplace(i, 1);} else {sens[i] = sens[i]+1;}}\
        vector<Int_t> v; for(pair<Int_t, Int_t> p: sens) {v.push_back(p.second);} return v;")

    #number of photoelectrons in sensor in event
    df = df.Define("qcal_nphotoel_sens", "nphotoel_sens(qcal_opdet_cell_id, qcal_opdet_is_detected)")

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )
    hx1 = rt.RDF.TH1DModel( ut.prepare_TH1D("hx1", xbin, xmin, xmax) )
    #hx = df.Histo1D(hx, "qcal_nphot_sens").GetValue()
    hx1 = df.Histo1D(hx1, "qcal_nphotoel_sens").GetValue()

    #print("nphot:", hx.GetEntries())
    #print("nphotel:", hx1.GetEntries())

    #ut.line_h1(hx, rt.kBlue, 2)
    #ut.line_h1(hx1, rt.kRed, 2)
    ut.line_h1(hx1, rt.kBlue, 2)

    if not draw: return hx1

    can = ut.box_canvas()

    #hx.Draw("")
    #hx1.Draw("same")
    hx1.Draw("")

    gPad.SetGrid()

    gPad.SetLogy()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")



    #hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )
    #hx = df.Histo1D(hx, "qcal_opdet_is_detected").GetValue()

    #can = ut.box_canvas()

    #ut.line_h1(hx, rt.kBlue, 2)

    #hx.Draw("")

    #gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    #can.SaveAs("01fig.pdf")



#main

#_____________________________________________________________________________
def make_pde():

    df = RDF.FromCSV("../../../../lmon2-data/pde_wavelength_30035.csv")

    df = df.Define("PDE_5V", "0.01*DetectionEfficiency_5V_percent")
    df = df.Define("PDE_2_5V", "0.01*DetectionEfficiency_2_5V_percent")

    g = df.Graph("Wavelength_nm", "PDE_2_5V")

    g.SetBit(TGraph.kIsSortedX)

    return g

#make_pde

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
























