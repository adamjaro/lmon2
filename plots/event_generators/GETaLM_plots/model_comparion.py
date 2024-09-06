#!/usr/bin/python3

# kinematics comparison among different models

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

    pitheta()

#main

#_____________________________________________________________________________
def pitheta():

    #mrad
    xbin1 = 0.01
    xbin2 = 0.05
    xbin3 = 0.3
    xmax = 4.5
    xmid1 = 0.2
    xmid2 = 1.1

    inp_lifbx = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/18x275/bh_lifbx_18x275_T3p3_20Mevt.hepmc3.tree.root"

    inp_zeus = "/media/jaroslav/SanDisk2T/datafiles/GETaLM_data/BH/div_test_2024_Sept06/18x275/bh_zeus_18x275_T3p3_20Mevt.hepmc3.tree.root"

    sigma_lifbx = 274.7913 # mb

    sigma_zeus = 276.3467 # mb

    #nmax = 1000000
    nmax = 0

    #hmod = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, 0, xmax) )
    hmod = rt.RDF.TH1DModel( ut.prepare_TH1D_3pt("hx", xbin1, xbin2, xbin3, 0, xmax, xmid1, xmid2) )

    #lifbx
    df = load_hepmc3_attrib(inp_lifbx, ["true_phot_theta", "beff_phot_theta"], 0, nmax)
    df = df.Define("true_phot_pitheta", "(TMath::Pi()-true_phot_theta)*1e3")
    df = df.Define("beff_phot_pitheta", "(TMath::Pi()-beff_phot_theta)*1e3")

    glifbx_true = ut.h1_to_arrays_norm( df.Histo1D(hmod, "true_phot_pitheta"), sigma_lifbx )
    glifbx_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_pitheta"), sigma_lifbx )

    #zeus
    df = load_hepmc3_attrib(inp_zeus, ["true_phot_theta", "beff_phot_theta"], 1, nmax)
    df = df.Define("true_phot_pitheta", "(TMath::Pi()-true_phot_theta)*1e3")
    df = df.Define("beff_phot_pitheta", "(TMath::Pi()-beff_phot_theta)*1e3")

    gzeus_true = ut.h1_to_arrays_norm( df.Histo1D(hmod, "true_phot_pitheta"), sigma_zeus )
    gzeus_beff = ut.h1_to_arrays_norm( df.Histo1D(hmod, "beff_phot_pitheta"), sigma_zeus )

    os.system("rm tmp*")


    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    plt.semilogy(gzeus_true[0], gzeus_true[1], "-", label="zeus_true", color="blue")
    plt.semilogy(gzeus_beff[0], gzeus_beff[1], "-.", label="zeus_beff", color="blue")

    plt.semilogy(glifbx_true[0], glifbx_true[1], "--", label="zeus_true", color="red")
    plt.semilogy(glifbx_beff[0], glifbx_beff[1], "-.", label="zeus_beff", color="red")

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#pitheta

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

















