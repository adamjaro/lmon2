#!/usr/bin/python3

# track resolution in kinematic quantities

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import RDataFrame

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 5

    func = {}
    func[0] = energy
    func[1] = pitheta
    func[2] = phi
    func[3] = logQ2
    func[4] = energy_err
    func[5] = theta_err

    func[iplot]()

#main

#_____________________________________________________________________________
def energy():

    #GeV
    ebin = 0.1
    emin = 3
    emax = 19

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_pass1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_pass1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    sel = det+"_is_rec==1"
    #sel = det+"_is_rec==1 && "+det+"_itrk==1"
    #sel = det+"_is_rec==1 && "+det+"_prim_id==1"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    #counts for all tracks
    hcnt = ut.prepare_TH1D("hcnt", 1, 0, 1)
    tree.Draw(det+"_x >> hcnt")
    print("All tracks:  ", hcnt.GetEntries())

    hxy = ut.prepare_TH2D("hxy", ebin, emin, emax, ebin, emin, emax)

    #tree.Draw(det+"_rec_en:"+det+"_true_en >> hxy", sel)
    tree.Draw(det+"_rec_en:"+det+"_mcp_en >> hxy", sel)

    print("On the plot:  ", hxy.Integral())

    ytit = "Reconstructed energy #it{E_{e}} (GeV)"
    xtit = "Generated MC particle energy #it{E_{e,mc}} (GeV)"
    ut.put_yx_tit(hxy, ytit, xtit, 1.4, 1.3)

    ut.set_margin_lbtr(gPad, 0.11, 0.11, 0.02, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetLogz()

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#energy

#_____________________________________________________________________________
def pitheta():

    #mrad
    xbin = 0.1
    xmin = 0
    xmax = 11

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_pass1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_pass1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    sel = det+"_is_rec==1"
    #sel = det+"_is_rec==1 && "+det+"_itrk==1"
    #sel = det+"_is_rec==1 && "+det+"_prim_id==1"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    #counts for all tracks
    hcnt = ut.prepare_TH1D("hcnt", 1, 0, 1)
    tree.Draw(det+"_x >> hcnt")
    print("All tracks:  ", hcnt.GetEntries())

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, xbin, xmin, xmax)

    #tree.Draw("(TMath::Pi()-"+det+"_rec_theta)*1e3:(TMath::Pi()-"+det+"_true_theta)*1e3 >> hxy", sel)
    tree.Draw("(TMath::Pi()-"+det+"_rec_theta)*1e3:(TMath::Pi()-"+det+"_mcp_theta)*1e3 >> hxy", sel)

    print("On the plot:  ", hxy.Integral())

    ytit = "Reconstructed #it{#pi-#theta_{e}} (mrad)"
    #xtit = "Generated true #it{#pi-#theta_{e,gen}} (mrad)"
    xtit = "Generated MC particle #it{#pi-#theta_{e,gen}} (mrad)"
    ut.put_yx_tit(hxy, ytit, xtit, 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetLogz()

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#pitheta

#_____________________________________________________________________________
def phi():

    #rad
    xbin = 0.1
    xmin = -TMath.Pi()-0.1
    xmax = TMath.Pi()+0.1

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_pass1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_pass1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    sel = det+"_is_rec==1"
    #sel = det+"_is_rec==1 && "+det+"_itrk==1"
    #sel = det+"_is_rec==1 && "+det+"_prim_id==1"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    #counts for all tracks
    hcnt = ut.prepare_TH1D("hcnt", 1, 0, 1)
    tree.Draw(det+"_x >> hcnt")
    print("All tracks:  ", hcnt.GetEntries())

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, xbin, xmin, xmax)

    tree.Draw(det+"_rec_phi:"+det+"_true_phi >> hxy", sel)
    #tree.Draw(det+"_rec_phi:"+det+"_mcp_phi >> hxy", sel)

    print("On the plot:  ", hxy.Integral())

    ytit = "Reconstructed #it{#phi_{e}} (rad)"
    xtit = "Generated true #it{#phi_{e,gen}} (rad)"
    ut.put_yx_tit(hxy, ytit, xtit, 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.1, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetLogz()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#phi

#_____________________________________________________________________________
def logQ2():

    #GeV^2
    #xbin = 0.1
    xbin = 0.05
    xmin = -8
    xmax = -1

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_pass1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_pass1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    sel = det+"_is_rec==1"
    #sel = det+"_is_rec==1 && "+det+"_itrk==1"
    #sel = det+"_is_rec==1 && "+det+"_prim_id==1"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    #counts for all tracks
    hcnt = ut.prepare_TH1D("hcnt", 1, 0, 1)
    tree.Draw(det+"_x >> hcnt")
    print("All tracks:  ", hcnt.GetEntries())

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, xbin, xmin, xmax)

    #tree.Draw("TMath::Log10("+det+"_rec_Q2):TMath::Log10("+det+"_true_Q2) >> hxy", sel)
    tree.Draw("TMath::Log10("+det+"_rec_Q2):TMath::Log10("+det+"_mcp_Q2) >> hxy", sel)

    print("On the plot:  ", hxy.Integral())

    ytit = "Reconstructed electron log_{10}(#it{Q}^{2}) (GeV^{2})"
    #xtit = "Generated true #it{Q}^{2} (GeV^{2})"
    xtit = "Generated MC particle log_{10}(#it{Q}^{2}) (GeV^{2})"
    ut.put_yx_tit(hxy, ytit, xtit, 1.9, 1.4)

    ut.set_margin_lbtr(gPad, 0.14, 0.11, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetLogz()

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#logQ2

#_____________________________________________________________________________
def energy_err():

    #relative error in energy

    #a.u.
    xbin = 0.02
    xmin = 0
    xmax = 0.2

    #inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/trk_ax2.root"
    inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/trk_bx2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_pass1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    #reconstructed value and its error
    val = det+"_rec_en"
    err = det+"_rec_en_err"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw(err+"/"+val+" >> hx", err+">=0")

    ytit = "Counts"
    xtit = "sigma_E/E"
    ut.put_yx_tit(hx, ytit, xtit, 1.4, 1.3)

    ut.set_margin_lbtr(gPad, 0.11, 0.11, 0.02, 0.11)

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#energy_err

#_____________________________________________________________________________
def theta_err():

    #relative error in theta

    #a.u.
    #xbin = 1e-5
    xmin = 0
    #xmax = 1e-3
    xbin = 0.001
    xmax = 0.01

    #inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/trk_ax2.root"
    inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/trk_bx2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_pass1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    #reconstructed value and its error
    val = "TMath::Pi()-"+det+"_rec_theta"
    err = "TMath::Pi()-"+det+"_rec_theta_err"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw(err+"/"+val+" >> hx", det+"_rec_theta_err>=0")

    ytit = "Counts"
    xtit = "sigma_theta/theta"
    ut.put_yx_tit(hx, ytit, xtit, 1.4, 1.3)

    ut.set_margin_lbtr(gPad, 0.11, 0.11, 0.02, 0.11)

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#theta_err

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()












