#!/usr/bin/python3

# kinematics, true generator and after beam effects

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import RDataFrame

import sys
sys.path.append("../../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 10

    func = {}
    func[0] = phot_en
    func[1] = phot_pitheta
    func[2] = phot_phi
    func[3] = el_en
    func[4] = el_pitheta
    func[5] = el_phi

    func[6] = phot_px
    func[7] = phot_py
    func[8] = phot_pz
    func[9] = el_px
    func[10] = el_py
    func[11] = el_pz

    func[iplot]()

#main

#_____________________________________________________________________________
def phot_en():

    #GeV
    xbin = 0.08
    xmin = 0.1
    xmax = 18.1

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_T3p3_6Mevt.root"

    df = RDataFrame("ltree", inp)

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )

    gtrue = ut.h1_to_graph( df.Histo1D(hx, "true_phot_en") )
    gbeff = ut.h1_to_graph( df.Histo1D(hx, "beff_phot_en") )

    ut.set_line_graph(gtrue, rt.kBlue, rt.kSolid, 3)
    ut.set_line_graph(gbeff, rt.kRed, rt.kDashed)

    can = ut.box_canvas()

    gtrue.Draw("al")
    gbeff.Draw("lsame")

    gPad.SetLogy()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#phot_en

#_____________________________________________________________________________
def phot_pitheta():

    #mrad
    xbin = 0.08
    xmin = 0
    xmax = 4.5
    #xmax = 3
    bin2 = 0.18
    xmid = 1.4

    #mbarn
    sigma = 274.7891

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    df = RDataFrame("ltree", inp)
    df = df.Define("true_phot_pitheta", "(TMath::Pi()-true_phot_theta)*1e3")
    df = df.Define("beff_phot_pitheta", "(TMath::Pi()-beff_phot_theta)*1e3")

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_mid("hx", xbin, bin2, xmin, xmax, xmid) )

    htrue = df.Histo1D(hx, "true_phot_pitheta")
    hbeff = df.Histo1D(hx, "beff_phot_pitheta")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    gtrue = ut.h1_to_graph(htrue, rt.kBlue, rt.kSolid, 3)
    gbeff = ut.h1_to_graph(hbeff, rt.kRed, rt.kDashed)

    can = ut.box_canvas()

    gtrue.Draw("al")
    gbeff.Draw("lsame")

    gPad.SetLogy()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#phot_pitheta

#_____________________________________________________________________________
def phot_phi():

    #rad
    nbins = 50
    xmin = -TMath.Pi()
    xmax = TMath.Pi()

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_n("hx", nbins, xmin, xmax) )

    htrue = df.Histo1D(hx, "true_phot_phi")
    hbeff = df.Histo1D(hx, "beff_phot_phi")

    #ut.norm_to_integral_list((htrue, hbeff), sigma)

    gtrue = ut.h1_to_graph(htrue, rt.kBlue, rt.kSolid, 3)
    gbeff = ut.h1_to_graph(hbeff, rt.kRed, rt.kDashed)

    can = ut.box_canvas()

    gbeff.Draw("al")
    gtrue.Draw("lsame")

    gPad.SetLogy()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#phot_phi

#_____________________________________________________________________________
def el_en():

    #GeV
    xbin = 0.08
    xmin = 0
    xmax = 18.1

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    df = RDataFrame("ltree", inp)

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )

    gtrue = ut.h1_to_graph( df.Histo1D(hx, "true_el_en") )
    gbeff = ut.h1_to_graph( df.Histo1D(hx, "beff_el_en") )

    ut.set_line_graph(gtrue, rt.kBlue, rt.kSolid, 3)
    ut.set_line_graph(gbeff, rt.kRed, rt.kDashed)

    can = ut.box_canvas()

    gtrue.Draw("al")
    gbeff.Draw("lsame")

    gPad.SetLogy()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#el_en

#_____________________________________________________________________________
def el_pitheta():

    #mrad
    xbin = 0.08
    xmin = 0
    xmax = 4.5
    bin2 = 0.18
    xmid = 1.4

    #mbarn
    sigma = 274.7891

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    df = RDataFrame("ltree", inp)
    df = df.Define("true_el_pitheta", "(TMath::Pi()-true_el_theta)*1e3")
    df = df.Define("beff_el_pitheta", "(TMath::Pi()-beff_el_theta)*1e3")

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_mid("hx", xbin, bin2, xmin, xmax, xmid) )

    htrue = df.Histo1D(hx, "true_el_pitheta")
    hbeff = df.Histo1D(hx, "beff_el_pitheta")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    gtrue = ut.h1_to_graph(htrue, rt.kBlue, rt.kSolid, 3)
    gbeff = ut.h1_to_graph(hbeff, rt.kRed, rt.kDashed)

    can = ut.box_canvas()

    gtrue.Draw("al")
    gbeff.Draw("lsame")

    gPad.SetLogy()

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#el_pitheta

#_____________________________________________________________________________
def el_phi():

    #rad
    nbins = 50
    xmin = -TMath.Pi()
    xmax = TMath.Pi()

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_n("hx", nbins, xmin, xmax) )

    htrue = df.Histo1D(hx, "true_el_phi")
    hbeff = df.Histo1D(hx, "beff_el_phi")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    gtrue = ut.h1_to_graph(htrue, rt.kBlue, rt.kSolid, 3)
    gbeff = ut.h1_to_graph(hbeff, rt.kRed, rt.kDashed)

    can = ut.box_canvas()
    frame = gPad.DrawFrame(xmin, 30, xmax, 60)
    frame.Draw()

    ut.put_yx_tit(frame, "d#it{#sigma}/d#it{#phi} (mb/rad)", "Electron azimuthal angle #it{#phi} (rad)", 1.5, 1.3)
    ut.set_margin_lbtr(gPad, 0.11, 0.1, 0.03, 0.02)

    gbeff.Draw("lsame")
    gtrue.Draw("lsame")

    leg = ut.prepare_leg(0.14, 0.83, 0.2, 0.1, 0.035) # 0.035
    leg.AddEntry(gtrue, "True generator", "l")
    leg.AddEntry(gbeff, "After divergence", "l")
    leg.Draw("same")

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#el_phi

#_____________________________________________________________________________
def phot_px():

    #MeV
    xmin = -60
    xmax = 60
    xbin = 0.5
    xbin1 = 6
    xbin3 = 6
    xmid1 = -15
    xmid2 = 15

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)
    df = df.Define("true_phot_px_MeV", "true_phot_px*1e3") # to MeV
    df = df.Define("beff_phot_px_MeV", "beff_phot_px*1e3") # to MeV

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_3pt("hx", xbin1, xbin, xbin3, xmin, xmax, xmid1, xmid2) )
    #hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )

    htrue = df.Histo1D(hx, "true_phot_px_MeV")
    hbeff = df.Histo1D(hx, "beff_phot_px_MeV")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    gtrue = ut.h1_to_arrays_centers( htrue )
    gbeff = ut.h1_to_arrays_centers( hbeff )

    #plt.plot(gtrue[0], gtrue[1], "-", label="gtrue", color="blue")
    plt.semilogy(gtrue[0], gtrue[1], "-", label="gtrue", color="blue")
    plt.semilogy(gbeff[0], gbeff[1], "--", label="gbeff", color="red")

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("Photon momentum $p_x$ (MeV)")
    ax.set_ylabel("d$\sigma$/d$p_x$ (mb/MeV)")

    leg = legend()
    leg.add_entry(leg_lin("blue"), "True generator")
    leg.add_entry(leg_lin("red", "--"), "After divergence")
    leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#phot_px

#_____________________________________________________________________________
def phot_py():

    #MeV
    xmin = -60
    xmax = 60
    xbin = 0.5
    xbin1 = 6
    xbin3 = 6
    xmid1 = -15
    xmid2 = 15

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)
    df = df.Define("true_phot_py_MeV", "true_phot_py*1e3") # to MeV
    df = df.Define("beff_phot_py_MeV", "beff_phot_py*1e3") # to MeV

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_3pt("hx", xbin1, xbin, xbin3, xmin, xmax, xmid1, xmid2) )

    htrue = df.Histo1D(hx, "true_phot_py_MeV")
    hbeff = df.Histo1D(hx, "beff_phot_py_MeV")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    gtrue = ut.h1_to_arrays_centers( htrue )
    gbeff = ut.h1_to_arrays_centers( hbeff )

    plt.semilogy(gtrue[0], gtrue[1], "-", label="gtrue", color="blue")
    plt.semilogy(gbeff[0], gbeff[1], "--", label="gbeff", color="red")

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("Photon momentum $p_y$ (MeV)")
    ax.set_ylabel("d$\sigma$/d$p_y$ (mb/MeV)")

    leg = legend()
    leg.add_entry(leg_lin("blue"), "True generator")
    leg.add_entry(leg_lin("red", "--"), "After divergence")
    leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#phot_py

#_____________________________________________________________________________
def phot_pz():

    #GeV
    xmin = -18
    xmax = 0
    xbin = 0.05

    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )

    htrue = df.Histo1D(hx, "true_phot_pz")
    hbeff = df.Histo1D(hx, "beff_phot_pz")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    gtrue = ut.h1_to_arrays_centers( htrue )
    gbeff = ut.h1_to_arrays_centers( hbeff )

    plt.semilogy(gtrue[0], gtrue[1], "-", label="gtrue", color="blue")
    plt.semilogy(gbeff[0], gbeff[1], "--", label="gbeff", color="red")

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("Photon momentum $p_z$ (GeV)")
    ax.set_ylabel("d$\sigma$/d$p_z$ (mb/GeV)")

    leg = legend()
    leg.add_entry(leg_lin("blue"), "True generator")
    leg.add_entry(leg_lin("red", "--"), "After divergence")
    leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#phot_pz

#_____________________________________________________________________________
def el_px():

    #MeV
    xmin = -62
    xmax = 62
    xbin = 0.5
    xbin1 = 8
    xbin3 = 8
    xmid1 = -15
    xmid2 = 15

    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)
    df = df.Define("true_el_px_MeV", "true_el_px*1e3") # to MeV
    df = df.Define("beff_el_px_MeV", "beff_el_px*1e3") # to MeV

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_3pt("hx", xbin1, xbin, xbin3, xmin, xmax, xmid1, xmid2) )

    htrue = df.Histo1D(hx, "true_el_px_MeV")
    hbeff = df.Histo1D(hx, "beff_el_px_MeV")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    gtrue = ut.h1_to_arrays_centers( htrue )
    gbeff = ut.h1_to_arrays_centers( hbeff )

    plt.semilogy(gtrue[0], gtrue[1], "-", label="gtrue", color="magenta") # blue
    plt.semilogy(gbeff[0], gbeff[1], "--", label="gbeff", color="red")

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("Electron momentum $p_x$ (MeV)")
    ax.set_ylabel("d$\sigma$/d$p_y$ (mb/MeV)")

    leg = legend()
    leg.add_entry(leg_lin("magenta"), "True generator") # blue
    leg.add_entry(leg_lin("red", "--"), "After divergence")
    leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#el_px

#_____________________________________________________________________________
def el_pz():

    #MeV
    xmin = -18
    xmax = 0
    xbin = 0.05

    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )

    htrue = df.Histo1D(hx, "true_phot_pz")
    hbeff = df.Histo1D(hx, "beff_phot_pz")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    gtrue = ut.h1_to_arrays_centers( htrue )
    gbeff = ut.h1_to_arrays_centers( hbeff )

    plt.semilogy(gtrue[0], gtrue[1], "-", label="gtrue", color="blue")
    plt.semilogy(gbeff[0], gbeff[1], "--", label="gbeff", color="red")

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("Photon momentum $p_z$ (GeV)")
    ax.set_ylabel("d$\sigma$/d$p_z$ (mb/GeV)")

    leg = legend()
    leg.add_entry(leg_lin("blue"), "True generator")
    leg.add_entry(leg_lin("red", "--"), "After divergence")
    leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#el_pz

#_____________________________________________________________________________
def el_py():

    #MeV
    xmin = -62
    xmax = 62
    xbin = 0.5
    xbin1 = 8
    xbin3 = 8
    xmid1 = -15
    xmid2 = 15

    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)
    df = df.Define("true_el_py_MeV", "true_el_py*1e3") # to MeV
    df = df.Define("beff_el_py_MeV", "beff_el_py*1e3") # to MeV

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_3pt("hx", xbin1, xbin, xbin3, xmin, xmax, xmid1, xmid2) )

    htrue = df.Histo1D(hx, "true_el_py_MeV")
    hbeff = df.Histo1D(hx, "beff_el_py_MeV")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    gtrue = ut.h1_to_arrays_centers( htrue )
    gbeff = ut.h1_to_arrays_centers( hbeff )

    plt.semilogy(gtrue[0], gtrue[1], "-", label="gtrue", color="blue")
    plt.semilogy(gbeff[0], gbeff[1], "--", label="gbeff", color="red")

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("Electron momentum $p_y$ (MeV)")
    ax.set_ylabel("d$\sigma$/d$p_y$ (mb/MeV)")

    leg = legend()
    leg.add_entry(leg_lin("blue"), "True generator")
    leg.add_entry(leg_lin("red", "--"), "After divergence")
    leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#el_py

#_____________________________________________________________________________
def el_pz():

    #GeV
    xmin = -18
    xmax = 0
    xbin = 0.05

    inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_lif_18x275_T3p3_12Mevt.root"
    #inp = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_July02/bh_zeus_18x275_4p5mrad_T3p3_12Mevt.root"

    #mbarn
    sigma = 274.7891

    df = RDataFrame("ltree", inp)

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )

    htrue = df.Histo1D(hx, "true_el_pz")
    hbeff = df.Histo1D(hx, "beff_el_pz")

    ut.norm_to_integral_list((htrue, hbeff), sigma)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    gtrue = ut.h1_to_arrays_centers( htrue )
    gbeff = ut.h1_to_arrays_centers( hbeff )

    plt.semilogy(gtrue[0], gtrue[1], "-", label="gtrue", color="blue")
    plt.semilogy(gbeff[0], gbeff[1], "--", label="gbeff", color="red")

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("Electron momentum $p_z$ (GeV)")
    ax.set_ylabel("d$\sigma$/d$p_z$ (mb/GeV)")

    leg = legend()
    leg.add_entry(leg_lin("blue"), "True generator")
    leg.add_entry(leg_lin("red", "--"), "After divergence")
    leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#el_pz

#_____________________________________________________________________________
def set_axes_color(ax, col):

    ax.xaxis.label.set_color(col)
    ax.yaxis.label.set_color(col)
    ax.tick_params(which = "both", colors = col)
    ax.spines["bottom"].set_color(col)
    ax.spines["left"].set_color(col)
    ax.spines["top"].set_color(col)
    ax.spines["right"].set_color(col)

#set_axes_color

#_____________________________________________________________________________
class legend:
    def __init__(self):
        self.items = []
        self.data = []
    def add_entry(self, i, d):
        self.items.append(i)
        self.data.append(d)
    def draw(self, px, col=None, **kw):
        leg = px.legend(self.items, self.data, **kw)
        if col is not None:
            px.setp(leg.get_texts(), color=col)
            if col != "black":
                leg.get_frame().set_edgecolor("orange")
        return leg

#_____________________________________________________________________________
def leg_lin(col, sty="-"):
    return Line2D([0], [0], lw=2, ls=sty, color=col)

#_____________________________________________________________________________
def leg_txt():
    return Line2D([0], [0], lw=0)

#_____________________________________________________________________________
def leg_dot(fig, col, siz=8):
    return Line2D([0], [0], marker="o", color=fig.get_facecolor(), markerfacecolor=col, markersize=siz)

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetFrameLineWidth(2)

    main()


















