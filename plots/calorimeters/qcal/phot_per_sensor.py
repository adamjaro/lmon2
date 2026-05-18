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

    energy = [1, 5, 9, 14, 18]
    #energy = [5, 9]

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
    ut.put_yx_tit(h1, "Counts", "Number of optical photons incident on sensor in event", 1.4, 1.3)

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

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#run_all

#_____________________________________________________________________________
def run_single(inp=None):

    draw = False
    if inp is None:
        #inp = "/home/jaroslav/sim/lmon2/macro/calorimeters/QCal2_fibers/lmon.root"
        inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal2bx1/en_1/lmon.root"
        draw = True

    xbin = 4
    xmin = 0
    xmax = 300

    df = RDataFrame("DetectorTree", inp)

    #number of photons in sensor in event
    df = df.Define("qcal_nphot_sens", "map<Int_t, Int_t> sens; for(Int_t i: qcal_opdet_cell_id) {\
        if(sens.find(i) == sens.end()) {sens.emplace(i, 1);} else {sens[i]++;}}\
        vector<Int_t> v; for(pair<Int_t, Int_t> p: sens) {v.push_back(p.second);} return v;")

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )
    hx = df.Histo1D(hx, "qcal_nphot_sens").GetValue()

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
def event_averages():

    #unused attempt for average counts over all events

    finp = TFile.Open(inp)
    nev = finp.Get("DetectorTree").GetEntries() # number of events
    finp.Close()

    #total number of sensors
    nsens = 16*11 # from xml

    can = ut.box_canvas()
    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_n("hx", nsens, 0, nsens) )

    #horizontal: sensor ID, vertical: number of photons in given sensor ID
    hx = df.Histo1D(hx, "qcal_opdet_cell_id").GetValue()
    ut.set_H1D(hx)

    #horizontal: number of photons in event, vertical: number of sensors
    hphot = ut.prepare_TH1D("hphot", xbin, xmin, xmax)

    cntsum = 0

    for ibin in range(1, hphot.GetNbinsX()+1):

        lower = hphot.GetBinLowEdge(ibin)
        upper = lower + hphot.GetBinWidth(ibin)

        cnt = 0

        #sensor loop
        for i in range(1, hx.GetNbinsX()+1):

            #num of photons in a given cell ID
            nid = hx.GetBinContent(i)

            if nid <= 0: continue

            nid = nid/nev

            #print(i, hx.GetBinLowEdge(i), nid)

            if nid >= lower and nid < upper: cnt += 1

        print(lower, upper, cnt)

        cntsum += lower*cnt

    print("cntsum:", cntsum)

    #hx.Draw("e1")

    #gPad.SetGrid()

    #gPad.SetLogy()

    #ut.invert_col(rt.gPad)
    #can.SaveAs("01fig.pdf")

#event_averages

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    #rt.EnableImplicitMT()

    main()







