#!/usr/bin/python3

# signal fraction for reconstructed tracks

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TGraphAsymmErrors, TH1D
from ROOT import TGaxis

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    #log_10(GeV^2)
    qbin = 0.2
    qmin = -10
    qmax = -0.5

    #kHz
    rmin = 1e-3
    rmax = 1e5

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag5dx12/trk_v1.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag5dx13/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag6ax2/trk_v1.root"

    det = "s1_tracks"
    #det = "s2_tracks"

    #18x275 GeV
    sigma_qr = 0.053266 # mb, quasi-real cross section, qr_bx_18x275_T3p3_10Mevt.log
    lumi_cmsec = 1.54e33 # cm^-2 sec^-1, instantaneous luminosity, 
    #rate_bkg = 22.6760075 # MHz, bunch frequency, qr_bx_18x275_T3p3_10Mevt.log
    rate_bkg = 22676.0075 # kHz, bunch frequency, qr_bx_18x275_T3p3_10Mevt.log

    #quasi-real production rate
    #rate_qr = sigma_qr*lumi_cmsec*1e-27*1e-3 # kHz
    rate_qr = 82.02964 # kHz, result of calculation above put directly here

    print("Quasi-real production rate (kHz):", rate_qr)

    infile = TFile.Open(inp)
    tree = infile.Get("event")

    #number of simulated events
    nsim = infile.Get("event").GetEntries()

    print("Num of simulated events:", nsim)

    #selection for signal and background tracks
    #sig_sel = "itrk==1"
    #bkg_sel = "itrk!=1"
    sig_sel = det+"_prim_id==1"
    bkg_sel = det+"_prim_id!=1"
    #bkg_sel = "prim_id!=1 && prim_id==itrk"

    #background
    hbkg = ut.prepare_TH1D("hbkg", qbin, qmin, qmax)
    tree.Draw("TMath::Log10("+det+"_rec_Q2) >> hbkg", bkg_sel)
    print("Background tracks per event:", hbkg.Integral()/nsim)
    #ut.norm_to_integral(hbkg, rate_bkg, rt.kBlue, True) # changed to manual example below

    #scale for background
    hbkg.Scale( rate_bkg/hbkg.Integral("width") )

    print("hbkg_int_w:", hbkg.Integral("width"))

    #quasi-real signal
    hsig = ut.prepare_TH1D("hsig", qbin, qmin, qmax)
    tree.Draw("TMath::Log10("+det+"_rec_Q2) >> hsig", sig_sel)
    print("Signal tracks per event:", hsig.Integral()/nsim)
    #ut.norm_to_integral(hsig, (hsig.GetEntries()/nsim)*rate_qr, rt.kBlue, True) # rate is in kHz

    #scale for signal
    hsig.Scale( ( rate_qr*hsig.GetEntries()/infile.Get("event").GetEntries() )/hsig.Integral("width") )

    #all events as signal + background
    hall = ut.prepare_TH1D("hall", qbin, qmin, qmax)
    hall.Add(hbkg)
    hall.Add(hsig)

    #signal fraction in all events
    sig_frac = TGraphAsymmErrors(hsig, hall);
    #sig_frac = ut.prepare_TH1D("sig_frac", qbin, qmin, qmax) # placeholder
    ut.set_graph(sig_frac)
    ut.set_H1D_col(sig_frac, rt.kGreen+1)

    #canvas frame for plot
    gStyle.SetPadTickY(1)
    can = ut.box_canvas(1086, 543) # square area is still 768^2
    can.SetMargin(0, 0, 0, 0)
    can.Divide(2, 1, 0, 0)
    gStyle.SetLineWidth(1)

    #plot on event rates
    can.cd(1)
    gPad.SetLogy()
    ut.set_margin_lbtr(gPad, 0.11, 0.11, 0.03, 0)

    frame = gPad.DrawFrame(qmin, rmin, qmax, rmax)
    ut.put_yx_tit(frame, "Event rate (kHz)", "Reconstructed #it{Q}^{2} (GeV^{2})", 1.5, 1.4)

    frame.Draw()

    #all events
    ut.line_h1(hall, rt.kBlue)
    hall.Draw("][ same")

    #signal
    ut.line_h1(hsig, rt.kRed)
    hsig.SetLineStyle(rt.kDashed)
    hsig.Draw("][ same")

    #legend on event rates
    leg = ut.prepare_leg(0.14, 0.77, 0.28, 0.16, 0.035) # 0.035
    tlab = {"s1_tracks": "Tagger 1", "s2_tracks": "Tagger 2"}
    leg.AddEntry("", tlab[det], "")
    leg.AddEntry(hall, "All tracks (sig + bkg)", "l")
    leg.AddEntry(hsig, "Signal only", "l")
    leg.Draw("same")

    gPad.SetGrid()

    ut.frame_pow10_labels(frame, -10, -1, "x", 1, 0.015)

    #plot on signal fraction
    can.cd(2)
    ut.set_margin_lbtr(gPad, 0, 0.11, 0.03, 0.12)

    fraction_frame = gPad.DrawFrame(-4.8, 0, -1.1, 1.1)
    fraction_frame.SetXTitle("Zoom on reconstructed #it{Q}^{2} (GeV^{2})")
    fraction_frame.GetXaxis().SetTitleOffset(1.4)

    fraction_frame.Draw()
    sig_frac.Draw("psame")

    #legend on signal fraction
    frac_leg = ut.prepare_leg(0.05, 0.9, 0.28, 0.05, 0.035)
    frac_leg.AddEntry(sig_frac, "Signal fraction in all tracks", "lp")
    frac_leg.Draw("same")

    gPad.SetGrid()

    ut.frame_pow10_labels_float(fraction_frame, -4.5, -1.5, "x", 0.5, 0.015)

    #vertical axis for fraction plot
    frac_xpos = fraction_frame.GetXaxis().GetXmax()
    frac_ypos = fraction_frame.GetMaximum()
    frac_ymin = fraction_frame.GetMinimum()

    frac_axis = TGaxis(frac_xpos, 0, frac_xpos, frac_ypos, frac_ymin, frac_ypos, 510, "+L")
    ut.set_axis(frac_axis)

    frac_axis.SetTitle("Signal fraction")
    frac_axis.SetTitleOffset(1.6)

    frac_axis.Draw()

    ut.invert_col_can(can)
    can.SaveAs("01fig.pdf")

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()

