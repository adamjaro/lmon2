#!/usr/bin/python3

#relative resolutions in reconstructed kinematics

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, RDataFrame, gInterpreter, TMath

from ROOT import RooFit as rf
from ROOT import RooRealVar, RooCBShape, RooDataHist, RooArgList

import sys
sys.path.append("../")
import plot_utils as ut
from parameter_descriptor import parameter_descriptor as pdesc

#_____________________________________________________________________________
def main():

    iplot = 2

    func = {}
    func[0] = pitheta_2d
    func[1] = pitheta_1d
    func[2] = lQ2_2d

    func[iplot]()

#main

#_____________________________________________________________________________
def pitheta_2d():

    #mrad
    xbin = 0.1
    xmin = 0
    xmax = 11

    #A.U.
    ybin = 0.1
    ymin = -2
    ymax = 20

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax3/trk_v4.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax4/trk_v2.root"

    #sel = "s1_tracks_is_rec==1 && s1_tracks_has_mcp==1"
    sel = "s1_tracks_is_rec==1 && s1_tracks_has_mcp==1 && (TMath::Pi()-s1_tracks_mcp_theta)*1e3<1"

    df = RDataFrame("event", inp)
    df = df.Define("s1_tracks_rec_pitheta", "(TMath::Pi()-s1_tracks_rec_theta["+sel+"])*1e3")
    df = df.Define("s1_tracks_mcp_pitheta", "(TMath::Pi()-s1_tracks_mcp_theta["+sel+"])*1e3")
    df = df.Define("s1_tracks_rel_mcp_pitheta", "(s1_tracks_rec_pitheta-s1_tracks_mcp_pitheta)/s1_tracks_mcp_pitheta")

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, ybin, ymin, ymax) )

    hx = df.Histo2D(hx, "s1_tracks_mcp_pitheta", "s1_tracks_rel_mcp_pitheta").GetValue()

    hx.Draw("colz")

    print("Entries: ", hx.GetEntries())

    ytit = "(rec-gen)/gen"
    xtit = "mcp pitheta (mrad)"
    ut.put_yx_tit(hx, ytit, xtit, 1.9, 1.4)

    ut.set_margin_lbtr(gPad, 0.14, 0.11, 0.03, 0.11)

    hx.SetMinimum(0.98)
    hx.SetContour(300)

    gPad.SetLogz()
    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#pitheta_2d

#_____________________________________________________________________________
def pitheta_1d():

    #A.U.
    xbin = 0.1
    xmin = -4
    xmax = 50

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax3/trk_v4.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax4/trk_v2.root"

    #sel = "s1_tracks_is_rec==1 && s1_tracks_has_mcp==1"
    sel = "s1_tracks_is_rec==1 && s1_tracks_has_mcp==1 && (TMath::Pi()-s1_tracks_mcp_theta)*1e3<1"
    sel2 = "s2_tracks_is_rec==1 && s2_tracks_has_mcp==1 && (TMath::Pi()-s2_tracks_mcp_theta)*1e3<1"

    df = RDataFrame("event", inp)
    df = df.Define("s1_tracks_rec_pitheta", "(TMath::Pi()-s1_tracks_rec_theta["+sel+"])*1e3")
    df = df.Define("s1_tracks_mcp_pitheta", "(TMath::Pi()-s1_tracks_mcp_theta["+sel+"])*1e3")
    df = df.Define("s1_tracks_rel_mcp_pitheta", "(s1_tracks_rec_pitheta-s1_tracks_mcp_pitheta)/s1_tracks_mcp_pitheta")

    df = df.Define("s2_tracks_rec_pitheta", "(TMath::Pi()-s2_tracks_rec_theta["+sel2+"])*1e3")
    df = df.Define("s2_tracks_mcp_pitheta", "(TMath::Pi()-s2_tracks_mcp_theta["+sel2+"])*1e3")
    df = df.Define("s2_tracks_rel_mcp_pitheta", "(s2_tracks_rec_pitheta-s2_tracks_mcp_pitheta)/s2_tracks_mcp_pitheta")

    can = ut.box_canvas()
    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )

    hx = df.Histo1D(hx, "s1_tracks_rel_mcp_pitheta").GetValue()
    #hx = df.Histo1D(hx, "s2_tracks_rel_mcp_pitheta").GetValue()

    ut.set_H1D(hx)

    #hx.Draw("e1")

    #reversed Crystal Ball
    res = RooRealVar("res", "res", xmin, xmax)
    mean = RooRealVar("mean", "mean", 0, -5, 5)
    n = RooRealVar("n", "n", 2, 0., 20.)
    sigma = RooRealVar("sigma", "sigma", 1, 0.001, 10)
    alpha = RooRealVar("alpha", "alpha", -0.3, -4.1, -0.01)
    cb = RooCBShape("cb", "reversed Crystal Ball", res, mean, sigma, alpha, n)

    #cb.fitTo(RooDataHist("dh", "dh", RooArgList(res), hx))
    hx = RooDataHist("dh", "dh", RooArgList(res), hx)
    cb.fitTo(hx)

    frame = res.frame()
    frame.SetTitle("")

    hx.plotOn(frame, rf.Name("data"))
    cb.plotOn(frame, rf.Precision(1e-6), rf.Name("CrystalBall"))

    frame.Draw()

    #fit parameters on the plot
    desc = pdesc(frame, 0.6, 0.7, 0.057); #x, y, sep
    desc.set_text_size(0.03)

    desc.itemD("#chi^{2}/ndf", frame.chiSquare("CrystalBall", "data", 4), -1)
    desc.itemR("mean", mean)
    desc.itemR("sigma", sigma)
    #desc.prec = 3
    desc.itemR("alpha", alpha)
    desc.itemR("n", n)
    desc.draw()

    #frame.SetMinimum(1.1)
    #frame.SetMaximum(2e5)

    ytit = "Counts"
    xtit = "(rec-gen)/gen"
    #ut.put_yx_tit(hx, ytit, xtit, 1.6, 1.4)

    ut.set_margin_lbtr(gPad, 0.11, 0.11, 0.03, 0.02)

    #gPad.SetLogy()
    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#pitheta_1d

#_____________________________________________________________________________
def lQ2_2d():

    #log_10(Q^2)
    xbin = 0.1
    xmin = -10
    xmax = 0

    #A.U.
    ybin = 0.02
    ymin = -1
    ymax = 1

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax2/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7ax3/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx1/trk_pass1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx2/trk_pass1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7bx3/trk_pass1.root"

    sel = "s1_tracks_is_rec==1 && s1_tracks_has_mcp==1"
    #sel = "s2_tracks_is_rec==1 && s2_tracks_has_mcp==1"

    df = RDataFrame("event", inp)
    df = df.Define("s1_tracks_mcp_good_Q2", "s1_tracks_mcp_Q2["+sel+"]")
    df = df.Define("s1_tracks_rec_good_Q2", "s1_tracks_rec_Q2["+sel+"]")
    #df = df.Define("s1_tracks_mcp_good_Q2", "s2_tracks_mcp_Q2["+sel+"]")
    #df = df.Define("s1_tracks_rec_good_Q2", "s2_tracks_rec_Q2["+sel+"]")
    df = df.Define("s1_tracks_mcp_lQ2",\
        "std::vector<Double_t> v; for(auto& i: s1_tracks_mcp_good_Q2) {v.push_back(TMath::Log10(i));} return v;")
    df = df.Define("s1_tracks_rec_lQ2",\
        "std::vector<Double_t> v; for(auto& i: s1_tracks_rec_good_Q2) {v.push_back(TMath::Log10(i));} return v;")

    df = df.Define("s1_tracks_rel_mcp_lQ2",\
        "std::vector<Double_t> v; for(size_t i=0; i<s1_tracks_rec_lQ2.size(); i++)\
        {v.push_back( (s1_tracks_rec_lQ2[i]-s1_tracks_mcp_lQ2[i])/s1_tracks_mcp_lQ2[i] );} return v;")

    gInterpreter.Declare("int calc_nout(std::vector<Double_t>& v) {\
        int n=0; for(auto& i: v) {if(i>0.3 || i<-0.3) n++;} return n;}")

    df = df.Define("s1_tracks_rel_mcp_lQ2_nout", "calc_nout(s1_tracks_rel_mcp_lQ2)")

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, ybin, ymin, ymax) )

    hx = df.Histo2D(hx, "s1_tracks_mcp_lQ2", "s1_tracks_rel_mcp_lQ2")
    #nOut = df.Filter("for(auto i: s1_tracks_rel_mcp_lQ2 ) return true;").Count()
    nout = df.Sum("s1_tracks_rel_mcp_lQ2_nout")
    ntrk_all = df.Sum("s1_ntrk")

    hx = hx.GetValue()
    hx.Draw("colz")

    #fraction outside 30%
    nall = hx.GetEntries()
    nsel = nout.GetValue()
    f30 = nsel/nall
    f30_sigma = f30*TMath.Sqrt( (nall-nsel)/(nall*nsel) )

    print("All tracks:", ntrk_all.GetValue())
    print("Reconstruction efficiency:", nall/ntrk_all.GetValue())

    #print("Entries: ", hx.GetEntries())
    #print("nout:", nout.GetValue())
    #print("Fraction outside 30%:", nout.GetValue()/hx.GetEntries())
    print("Entries: ", nall, ", nout:", nsel, ", fraction outside 30%:", f30, "+/-", f30_sigma)

    ytit = "(rec-gen)/gen"
    xtit = "mcp lQ2"
    ut.put_yx_tit(hx, ytit, xtit, 1.9, 1.4)

    ut.set_margin_lbtr(gPad, 0.14, 0.11, 0.03, 0.11)

    hx.SetMinimum(0.98)
    hx.SetContour(300)

    gPad.SetLogz()
    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#lQ2_2d


#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    rt.EnableImplicitMT()

    main()










