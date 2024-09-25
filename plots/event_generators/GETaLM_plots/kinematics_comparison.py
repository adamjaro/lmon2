#!/usr/bin/python3

# comparison in kinematic variables for several generator configuration

import os
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

    #phot_px()
    #phot_py()
    phot_phi()

#main

#_____________________________________________________________________________
def phot_px():

    #MeV
    xbin1 = 0.5
    #xbin2 = 0.07
    #xbin3 = 0.3
    xmin = -60
    xmax = 60
    #xmid1 = 0.2
    #xmid2 = 1.1

    nmax = 100000
    #nmax = 0


    sigma_lifbx = 214.879 # mb, 10x100 GeV

    #nominal divergence
    inp_nom = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_T3p3_20Mevt.hepmc3.tree.root"

    hmod = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin1, xmin, xmax) )
    #hmod = rt.RDF.TH1DModel( ut.prepare_TH1D_3pt("hx", xbin1, xbin2, xbin3, 0, xmax, xmid1, xmid2) )

    #nominal
    df = load_hepmc3_attrib(inp_nom, ["true_phot_px", "beff_phot_px"], 0, nmax)
    df = df.Define("true_phot_px_MeV", "true_phot_px*1e3") # to MeV
    df = df.Define("beff_phot_px_MeV", "beff_phot_px*1e3") # to MeV

    gnom_true = ut.h1_to_arrays_norm( df.Histo1D(hmod, "true_phot_px_MeV"), sigma_lifbx )
    gnom_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_px_MeV"), sigma_lifbx )

    #zeus
    #df = load_hepmc3_attrib(inp_zeus, ["true_phot_theta", "beff_phot_theta"], 1, nmax)
    #df = df.Define("true_phot_pitheta", "(TMath::Pi()-true_phot_theta)*1e3")
    #df = df.Define("beff_phot_pitheta", "(TMath::Pi()-beff_phot_theta)*1e3")

    #gzeus_true = ut.h1_to_arrays_norm( df.Histo1D(hmod, "true_phot_pitheta"), sigma_zeus )
    #gzeus_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_pitheta"), sigma_zeus )

    os.system("rm tmp*")


    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    #ax.set_xlabel(r"Photon polar angle $\pi-\theta_\gamma$ (mrad)")
    #ax.set_ylabel(r"Cross section $\frac{\mathrm{d}\sigma}{\mathrm{d}\theta_\gamma}$ (mb/mrad)")

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    plt.semilogy(gnom_true[0], gnom_true[1], "-", label="nom_true", color="blue")
    plt.semilogy(gnom_beff[0], gnom_beff[1], "--", label="nom_beff", color="red")

    #plt.semilogy(glifbx_true[0], glifbx_true[1], "--", label="lifbx_true", color="red")
    #plt.semilogy(glifbx_beff[0], glifbx_beff[1], ":", label="lifbx_beff", color="red")

    leg = legend()
    leg.add_entry(leg_txt(), "18x275 GeV")
    leg.add_entry(leg_lin("blue", "-"), "ZEUS true")
    leg.add_entry(leg_lin("blue", "-."), "ZEUS with divergence")
    leg.add_entry(leg_lin("red", "--"), "Lifshitz true")
    leg.add_entry(leg_lin("red", ":"), "Lifshitz with divergence")
    #leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#phot_px

#_____________________________________________________________________________
def phot_py():


    #MeV
    xbin1 = 3
    xbin2 = 0.3
    xbin3 = 3
    xmin = -22
    xmax = 22
    xmid1 = -5
    xmid2 = 5

    #inner -5, 5, outer bins 3

    nmax = 1000000
    #nmax = 0

    sigma_lifbx = 214.879 # mb, 10x100 GeV

    #nominal divergence
    inp_nom = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_T3p3_20Mevt.hepmc3.tree.root"

    #Div4
    inp_div4 = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_Div4_20Mevt.hepmc3.tree.root"

    #Div2
    inp_div2 = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_Div2_20Mevt.hepmc3.tree.root"

    #hmod = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin1, xmin, xmax) )
    hmod = rt.RDF.TH1DModel( ut.prepare_TH1D_3pt("hx", xbin1, xbin2, xbin3, xmin, xmax, xmid1, xmid2) )

    #nominal
    df = load_hepmc3_attrib(inp_nom, ["true_phot_py", "beff_phot_py"], 0, nmax)
    df = df.Define("true_phot_py_MeV", "true_phot_py*1e3") # to MeV
    df = df.Define("beff_phot_py_MeV", "beff_phot_py*1e3") # to MeV

    gnom_true = ut.h1_to_arrays_norm( df.Histo1D(hmod, "true_phot_py_MeV"), sigma_lifbx )
    gnom_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_py_MeV"), sigma_lifbx )

    #div4
    df = load_hepmc3_attrib(inp_div4, ["beff_phot_py"], 1, nmax)
    df = df.Define("beff_phot_py_MeV", "beff_phot_py*1e3") # to MeV

    gdiv4_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_py_MeV"), sigma_lifbx )

    #div2
    df = load_hepmc3_attrib(inp_div2, ["beff_phot_py"], 2, nmax)
    df = df.Define("beff_phot_py_MeV", "beff_phot_py*1e3") # to MeV

    gdiv2_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_py_MeV"), sigma_lifbx )

    os.system("rm tmp*")

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    #ax.set_xlabel(r"Photon polar angle $\pi-\theta_\gamma$ (mrad)")
    #ax.set_ylabel(r"Cross section $\frac{\mathrm{d}\sigma}{\mathrm{d}\theta_\gamma}$ (mb/mrad)")

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    plt.semilogy(gnom_true[0], gnom_true[1], "-", label="nom_true", color="blue")
    plt.semilogy(gnom_beff[0], gnom_beff[1], "--", label="nom_beff", color="red")
    plt.semilogy(gdiv4_beff[0], gdiv4_beff[1], "-.", label="div4_beff", color="goldenrod")
    plt.semilogy(gdiv2_beff[0], gdiv2_beff[1], ":", label="div2_beff", color="orange")

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#phot_py

#_____________________________________________________________________________
def phot_phi():

    #rad
    nbins = 50
    xmin = -TMath.Pi()
    xmax = TMath.Pi()

    #50 bins, tried: 100, 60

    nmax = 1000000
    #nmax = 0

    sigma_lifbx = 214.879 # mb, 10x100 GeV

    #nominal divergence
    inp_nom = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_T3p3_20Mevt.hepmc3.tree.root"

    #Div4
    inp_div4 = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_Div4_20Mevt.hepmc3.tree.root"

    #Div2
    inp_div2 = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_Div2_20Mevt.hepmc3.tree.root"

    #Div1
    inp_div1 = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_Div1_20Mevt.hepmc3.tree.root"

    #Div3
    inp_div3 = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/10x100/bh_lifbx_10x100_Div3_20Mevt.hepmc3.tree.root"

    hmod = rt.RDF.TH1DModel( ut.prepare_TH1D_n("hx", nbins, xmin, xmax) )

    #nominal
    df = load_hepmc3_attrib(inp_nom, ["true_phot_phi", "beff_phot_phi"], 0, nmax)
    gnom_true = ut.h1_to_arrays_norm( df.Histo1D(hmod, "true_phot_phi"), sigma_lifbx )
    gnom_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_phi"), sigma_lifbx )

    #div4
    df = load_hepmc3_attrib(inp_div4, ["beff_phot_phi"], 1, nmax)
    gdiv4_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_phi"), sigma_lifbx )

    #div2
    df = load_hepmc3_attrib(inp_div2, ["beff_phot_phi"], 2, nmax)
    gdiv2_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_phi"), sigma_lifbx )

    #div1
    df = load_hepmc3_attrib(inp_div1, ["beff_phot_phi"], 3, nmax)
    gdiv1_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_phi"), sigma_lifbx )

    #div3
    df = load_hepmc3_attrib(inp_div3, ["beff_phot_phi"], 4, nmax)
    gdiv3_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_phi"), sigma_lifbx )

    os.system("rm tmp*")

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    plt.semilogy(gnom_true[0], gnom_true[1], "-", label="nom_true", color="blue")
    plt.semilogy(gnom_beff[0], gnom_beff[1], "--", label="nom_beff", color="red")
    plt.semilogy(gdiv4_beff[0], gdiv4_beff[1], "-.", label="div4_beff", color="goldenrod")
    plt.semilogy(gdiv2_beff[0], gdiv2_beff[1], ":", label="div2_beff", color="orange")
    plt.semilogy(gdiv1_beff[0], gdiv1_beff[1], "-.", label="div1_beff", color="darkgreen")
    plt.semilogy(gdiv3_beff[0], gdiv3_beff[1], ":", label="div3_beff", color="magenta")

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#phot_phi

#_____________________________________________________________________________
def load_hepmc3_attrib(inp, attrib, icall=0, nmax=0):

    gSystem.AddIncludePath(" -I/home/jaroslav/sim/hepmc/HepMC3-3.3.0/install/include")

    #chr(97) = 'a', chr(122) = 'z'

    #attribute character identifier and name in hepmc
    iattr = 0
    attr = {}
    for i in attrib:
        attr[chr(97+iattr)] = i
        iattr += 1

    #hepmc input and tree output
    tmp = open("tmp"+str(icall)+".C", "w")
    tmp.write('#include <iostream>\n #include "HepMC3/ReaderRootTree.h"\n \
        using namespace std; using namespace HepMC3;\
        void tmp'+str(icall)+'() {\
            shared_ptr<ReaderRootTree> read = make_shared<ReaderRootTree>("'+inp+'");\
            TFile out("tmp.root", "recreate");\
            TTree otree("tmp", "tmp");')

    #attribute branches by character identifiers and names
    for i in attr:
        tmp.write('\
            Double_t '+i+';\
            otree.Branch("'+attr[i]+'", &'+i+', "'+attr[i]+'/D");\
        ')

    #loop over hepmc events
    tmp.write('\
            Long64_t iev=0;\
            while( !read->failed() ) {\
                if('+str(nmax)+' > 0 and iev >= '+str(nmax)+') break;\
                GenEvent mc(Units::GEV,Units::MM);\
                read->read_event(mc);\
                if( read->failed() ) break;')

    #load the attributes by names to variables by identifiers
    for i in attr:
        tmp.write('\
                '+i+' = mc.attribute<DoubleAttribute>("'+attr[i]+'")->value();\
        ')

    #fill the output tree and finish the loop
    tmp.write('\
                otree.Fill();\
                iev++;\
            }\
            otree.Write();\
            out.Close();\
        }')
    tmp.close()

    gROOT.ProcessLine(".x tmp"+str(icall)+".C+g")

    df = RDataFrame("tmp", "tmp.root")

    return df

#load_hepmc3_attrib

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

















