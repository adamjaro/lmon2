#!/usr/bin/python3

# acceptance calculation

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import RDataFrame

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 1

    func = {}
    func[0] = energy_1d
    func[1] = energy_pitheta
    func[2] = logx_logQ2
    func[3] = energy_mlt
    func[4] = energy_1d_s12

    func[iplot]()

#main

#_____________________________________________________________________________
def energy_1d():

    #GeV
    xmin = 0
    xmax = 19

    amin = 0.01
    amax = 2

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6a/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/trk.root"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    acc = rt.acalc(tree, "true_el_E", "s1_ntrk", "s2_ntrk")
    acc.prec = 0.6
    acc.bmin = 0.1
    #as2.nev = int(1e4)
    gacc = acc.get()

    can = ut.box_canvas()
    frame = gPad.DrawFrame(xmin, amin, xmax, amax)

    ut.put_yx_tit(frame, "Tagger acceptance", "Electron energy #it{E} (GeV)", 1.6, 1.3)

    frame.Draw()

    ut.set_margin_lbtr(gPad, 0.11, 0.1, 0.03, 0.02)

    ut.set_graph(gacc, rt.kRed)
    gacc.Draw("psame")

    gPad.SetGrid()

    gPad.SetLogy()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#energy_1d

#_____________________________________________________________________________
def energy_pitheta():

    #reconstruction efficiency in energy (GeV) and pi - theta (mrad)

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag11ax1/tracks_lps_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ax3/tracks_v3_mcp.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10bx1/tracks_v0.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10dx1/tracks_v0.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10cx1/tracks_v0.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10cx2/tracks_v0.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ex1/tracks_v0.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ex3/tracks_v0.root"

    #tagger 1 or 2
    tag = 2

    #bins in energy, GeV
    xbin = 0.3
    xmin = 3.5
    xmax = 20

    #bins in theta, mrad
    ybin = 0.2
    ymin = 0
    ymax = 11

    if tag == 1:
        sel = "s1_ntrk>0"
        #sel = "s1_nrec>0"
        #sel = "s1_ntrk>0 && s1_nrec==0"
        lab_sel = "Tagger 1"
    elif tag == 2:
        sel = "s2_ntrk>0"
        #sel = "s2_nrec>0"
        lab_sel = "Tagger 2"
    elif tag == 3:
        sel = "s1_ntrk>0 || s2_ntrk>0"
        lab_sel = ""

    df = RDataFrame("event", inp)
    df = df.Define("pitheta", "(TMath::Pi()-true_el_theta)*1e3")
    df = df.Define("s1_nrec", "int n=0; for(auto& i:s1_tracks_is_rec) { if(i) n++; } return n;")
    df = df.Define("s2_nrec", "int n=0; for(auto& i:s2_tracks_is_rec) { if(i) n++; } return n;")
    df = df.Define("s1_ntrk", "s1_tracks_x.size()")
    df = df.Define("s2_ntrk", "s2_tracks_x.size()")

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, ybin, ymin, ymax) )

    hAll = df.Histo2D(hx, "true_el_E", "pitheta")
    #hAll = df.Filter("s1_ntrk>0").Histo2D(hx, "true_el_E", "pitheta")
    hTag = df.Filter(sel).Histo2D(hx, "true_el_E", "pitheta")

    hAll = hAll.GetValue()
    hTag = hTag.GetValue()

    print("all: ", hAll.GetEntries())
    print("tag: ", hTag.GetEntries())
    print("rto: ", hTag.GetEntries()/hAll.GetEntries())

    hTag.Divide(hAll)

    xtit = "Electron energy #it{E} (GeV)"
    ytit = "Electron polar angle #it{#pi}-#it{#theta} (mrad)"
    ut.put_yx_tit(hTag, ytit, xtit, 1.4, 1.3)

    hTag.SetTitleOffset(1.4, "Z")
    hTag.SetZTitle("Acceptance")

    ut.set_margin_lbtr(gPad, 0.1, 0.1, 0.015, 0.15)

    #gPad.SetLogz()

    gPad.SetGrid()

    hTag.SetMinimum(0.)
    hTag.SetMaximum(1.)
    hTag.SetContour(300)

    hTag.Draw("colz")

    leg = ut.prepare_leg(0.1, 0.85, 0.24, 0.1, 0.035) # x, y, dx, dy, tsiz
    leg.AddEntry("", "#bf{"+lab_sel+"}", "")
    leg.AddEntry("", "Dist. from beam axis: 5 cm", "")
    #leg.AddEntry("", "Dist. from beam axis: 4 cm", "")
    #leg.Draw("same")

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#energy_pitheta

#_____________________________________________________________________________
def logx_logQ2():

    #acceptance x efficiency in Bjorken x (horizontal) and Q^2 (vertical), both as log_10

    #log_10(x)
    xbin = 0.05
    xmin = -12
    xmax = -1

    #log_10(GeV^2)
    ybin = 0.05
    ymin = -8
    ymax = 0

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ax1/acc.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ax3/tracks_v2_tracks_only.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ax3/tracks_v3_mcp.root"

    df = RDataFrame("event", inp)
    df = df.Define("s1_ntrk", "s1_tracks_x.size()")
    df = df.Define("s2_ntrk", "s2_tracks_x.size()")
    df = df.Define("logQ2", "TMath::Log10(true_Q2)")
    df = df.Define("logx", "TMath::Log10(true_x)")

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, ybin, ymin, ymax) )

    hxy_all = df.Histo2D(hx, "logx", "logQ2").GetValue()
    hxy_sel = df.Filter("s1_ntrk>0").Histo2D(hx, "logx", "logQ2").GetValue()
    #hxy_sel = df.Filter("s2_ntrk>0").Histo2D(hx, "logx", "logQ2").GetValue()

    print("Selected: ", hxy_sel.GetEntries())

    hxy_sel.Divide(hxy_all)

    ytit = "Virtuality log_{10}(#it{Q}^{2}) (GeV^{2})"
    xtit = "Bjorken log_{10}(#it{x})"
    ut.put_yx_tit(hxy_sel, ytit, xtit, 1.35, 1.3)

    #hxy_sel.SetZTitle("Acceptance #times Efficiency #it{A}#times#it{E}")
    hxy_sel.SetZTitle("Acceptance")
    hxy_sel.SetTitleOffset(1.4, "Z")

    ut.set_margin_lbtr(gPad, 0.11, 0.11, 0.03, 0.16)

    gPad.SetGrid()

    hxy_sel.SetMinimum(0)
    hxy_sel.SetMaximum(1)
    hxy_sel.SetContour(300)

    hxy_sel.Draw("colz")

    leg = ut.prepare_leg(0.12, 0.85, 0.2, 0.1, 0.035) # 0.035
    #leg.AddEntry("", "Tagger 1, V6.3", "")
    #leg.AddEntry("", "Tagger 2, V6.3", "")
    #leg.Draw("same")

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#logx_logQ2

#_____________________________________________________________________________
def energy_mlt():

    #reconstruction efficiency in energy (GeV) and mlt defined as -log_10(pi - theta) (rad)

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx3/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag10ax3/tracks_v3_mcp.root"

    #bins in energy, GeV
    xbin = 0.3
    xmin = 1
    xmax = 20

    #bins in mlt (theta in rad)
    ybin = 0.05
    ymin = 1.1
    ymax = 6

    df = RDataFrame("event", inp)
    df = df.Define("mlt", "-TMath::Log10(TMath::Pi()-true_el_theta)")
    df = df.Define("s1_ntrk", "s1_tracks_x.size()")
    df = df.Define("s2_ntrk", "s2_tracks_x.size()")
    df = df.Define("s1_nrec", "int n=0; for(auto& i:s1_tracks_is_rec) { if(i) n++; } return n;")

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, ybin, ymin, ymax) )

    hAll = df.Histo2D(hx, "true_el_E", "mlt")
    #hTag = df.Filter("s1_ntrk>0").Histo2D(hx, "true_el_E", "mlt")
    hTag = df.Filter("s2_ntrk>0").Histo2D(hx, "true_el_E", "mlt")

    hTag = hTag.GetValue()

    hTag.Divide(hAll.GetValue())

    xtit = "Electron energy #it{E} (GeV)"
    ytit = "mlt"
    ut.put_yx_tit(hTag, ytit, xtit, 1.4, 1.3)

    hTag.SetTitleOffset(1.4, "Z")
    hTag.SetZTitle("Acceptance")

    ut.set_margin_lbtr(gPad, 0.1, 0.1, 0.015, 0.15)

    #gPad.SetLogz()

    gPad.SetGrid()

    hTag.SetMinimum(0.)
    hTag.SetMaximum(1.)
    hTag.SetContour(300)

    hTag.Draw("colz")

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#energy_mlt

#_____________________________________________________________________________
def energy_1d_s12():

    #acceptance on energy for both taggers in one plot

    #GeV
    xmin = 0
    xmax = 18.5

    amin = 1e-6
    amax = 2

    #inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/trk.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax3/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax4/trk_v0.root"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    as1 = rt.acalc(tree, "true_el_E", "s1_ntrk")
    #as1.prec = 0.6
    #as1.bmin = 0.1
    as1.prec = 0.3
    as1.bmin = 0.01
    #as1.nev = int(1e4)
    ga1 = as1.get()

    as2 = rt.acalc(tree, "true_el_E", "s2_ntrk")
    #as2.prec = 0.4
    #as2.bmin = 0.01
    as2.prec = 0.3
    as2.bmin = 0.01
    #as2.nev = int(1e4)
    ga2 = as2.get()

    can = ut.box_canvas()
    frame = gPad.DrawFrame(xmin, amin, xmax, amax)

    ut.put_yx_tit(frame, "Tagger acceptance", "Electron energy #it{E} (GeV)", 1.5, 1.3)

    frame.Draw()

    ut.set_margin_lbtr(gPad, 0.11, 0.1, 0.03, 0.02)

    ut.set_graph(ga1, rt.kBlue)
    ga1.Draw("psame")

    ut.set_graph(ga2, rt.kRed)
    ga2.Draw("psame")

    #leg = ut.prepare_leg(0.48, 0.68, 0.2, 0.1, 0.035) # 0.035
    leg = ut.prepare_leg(0.14, 0.83, 0.2, 0.1, 0.035) # 0.035
    leg.AddEntry(ga1, "Tagger 1", "lp")
    leg.AddEntry(ga2, "Tagger 2", "lp")
    leg.Draw("same")

    gPad.SetGrid()

    gPad.SetLogy()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#energy_1d

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    #build and load the acceptance calculator
    gROOT.ProcessLine(".L ../acalc.h+")

    main()




















