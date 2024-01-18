#!/usr/bin/python3

# track basic parameters

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 8

    func = {}
    func[0] = theta_x
    func[1] = theta_y
    func[2] = xy
    func[3] = chi2
    func[4] = lQ2_resolution
    func[5] = theta_x_x
    func[6] = theta_x_true_en
    func[7] = track_en_cal_en
    func[8] = ntrk

    func[iplot]()

#main

#_____________________________________________________________________________
def theta_x():

    #mrad
    xbin = 1
    xmin = -10
    xmax = 40

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx1/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx2/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx3/trk_v1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw(det+"_theta_x*1e3 >> hx")

    ut.put_yx_tit(hx, "Counts", "Tracks #it{#theta}_{x} angle (mrad)", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#theta_x

#_____________________________________________________________________________
def theta_y():

    #mrad
    xbin = 0.4
    xmin = -12
    xmax = 12

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx1/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx2/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx3/trk_v1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw(det+"_theta_y*1e3 >> hx")

    ut.put_yx_tit(hx, "Counts", "Tracks #it{#theta_{y}} angle (mrad)", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#theta_y

#_____________________________________________________________________________
def xy():

    #track position in x and y

    #mrad
    xbin = 1
    xmin = -80
    xmax = 80

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx1/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx2/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx3/trk_v1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, xbin, xmin, xmax)

    tree.Draw(det+"_y:"+det+"_x >> hxy")

    ut.put_yx_tit(hxy, "#it{y} (mm)", "#it{x} (mm)", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#xy

#_____________________________________________________________________________
def chi2():

    xbin = 0.003
    xmin = 0
    xmax = 0.06

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_v2.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx1/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx3/trk_v1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw(det+"_chi2_xy >> hx")

    ut.put_yx_tit(hx, "Counts", "Tracks #it{#chi}^{2}", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    gPad.SetLogy()

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#chi2

#_____________________________________________________________________________
def lQ2_resolution():

    #track resolution in log_10(Q^2)

    #GeV^2
    xbin = 0.1
    #xmin = -8
    #xmax = 1
    xmin = -10
    xmax = 0

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax3/trk_v1.root" # quasi-real single particle
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax4/trk_v1.root" # bremsstrahlung single particle

    #det = "s1_tracks"
    det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, xbin, xmin, xmax)

    #tree.Draw("TMath::Log10("+det+"_rec_Q2):TMath::Log10("+det+"_true_Q2) >> hxy")
    tree.Draw("TMath::Log10("+det+"_rec_Q2):TMath::Log10("+det+"_mcp_Q2) >> hxy")

    xtit = "MC particle log_{10}(#it{Q}^{2}) (GeV^{2})"
    ytit = "Reconstructed track log_{10}(#it{Q}^{2}) (GeV^{2})"
    ut.put_yx_tit(hxy, ytit, xtit, 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetLogz()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#lQ2_resolution

#_____________________________________________________________________________
def theta_x_x():

    #track theta_x angle and position in x

    #mrad on y axis
    ybin = 0.5
    ymin = -6
    ymax = 30

    #mm on x axis
    xbin = 2
    xmin = -80
    xmax = 80

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax3/trk_v4.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, ybin, ymin, ymax)

    tree.Draw("("+det+"_theta_x)*1e3:"+det+"_x >> hxy")

    xtit = "Track horizontal position #it{x} (mm)"
    ytit = "Track horizontal angle #it{#theta}_{x} (mrad)"
    ut.put_yx_tit(hxy, ytit, xtit, 1.5, 1.3)

    ut.set_margin_lbtr(gPad, 0.12, 0.12, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetGrid()

    gPad.SetLogz()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#theta_x_x

#_____________________________________________________________________________
def theta_x_true_en():

    #track theta_x angle and true generated energy

    #mrad on y axis
    ybin = 0.5
    ymin = -10
    ymax = 30

    #GeV on x axis
    xbin = 0.3
    xmin = 2
    xmax = 19

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax3/trk_v4.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, ybin, ymin, ymax)

    tree.Draw("("+det+"_theta_x)*1e3:"+det+"_true_en >> hxy")

    xtit = "Generated true electron energy (GeV)"
    ytit = "Track horizontal angle #it{#theta}_{x} (mrad)"
    ut.put_yx_tit(hxy, ytit, xtit, 1.5, 1.3)

    ut.set_margin_lbtr(gPad, 0.12, 0.12, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetGrid()

    gPad.SetLogz()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#theta_x_true_en

#_____________________________________________________________________________
def track_en_cal_en():

    #track reconstructed energy (x) and calorimeter energy (y)

    #GeV on y axis
    #ybin = 0.01
    #ymin = 0
    #ymax = 0.7
    ybin = 0.2
    ymin = 2
    ymax = 16

    #GeV on x axis
    xbin = 0.2
    xmin = 5
    xmax = 15

    #inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/trk.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_v3.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag8ax1/trk_v1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hxy = ut.prepare_TH2D("hxy", xbin, xmin, xmax, ybin, ymin, ymax)

    tree.Draw(det+"_cal_en:"+det+"_rec_en >> hxy", det+"_cal_x<64") # 65

    xtit = "Track reconstructed energy (GeV)"
    ytit = "Deposited energy in calorimeter (GeV)"
    ut.put_yx_tit(hxy, ytit, xtit, 1.5, 1.3)

    ut.set_margin_lbtr(gPad, 0.12, 0.1, 0.03, 0.11)

    hxy.SetMinimum(0.98)
    hxy.SetContour(300)

    gPad.SetGrid()

    gPad.SetLogz()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#theta_x_true_en

#_____________________________________________________________________________
def ntrk():

    #mrad
    xbin = 1
    xmin = 0
    xmax = 20

    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax5/trk_v1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw("s1_ntrk >> hx") #, "s1_ntrk==14"

    ut.line_h1(hx)

    ut.put_yx_tit(hx, "Counts", "Number of tracks per event", 1.9, 1.3)

    ut.set_margin_lbtr(gPad, 0.14, 0.12, 0.03, 0.03)

    gPad.SetGrid()

    gPad.SetLogy()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#ntrk

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()








