#!/usr/bin/env python3

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, RDataFrame

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import mplhep as hep

from photoelectrons import photoelectrons

import sys
sys.path.append("../..")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    #inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal3cx3/en_9/lmon.root"
    inp = "/home/jaroslav/sim/lmon2/macro/calorimeters/QCal2_fibers/fib_1mm_sens_3mm/lmon.root"

    xbin = 10
    xmin = 200
    xmax = 1400

    df = RDataFrame("DetectorTree", inp)

    #wavelength from energy, E (eV) = 1239.841 / lambda (nm), PDG 2024 Tab 1.1
    df = df.Define("qcal_opdet_phot_lambda", "1239.841/qcal_opdet_phot_en")

    #detected flag
    df = df.Define("qcal_opdet_is_detected", "sens(qcal_opdet_phot_en)")

    #wavelength for detected photons
    df = df.Define("qcal_opdet_phot_lambda_det", "vector<double> v;\
        for(size_t i=0; i<qcal_opdet_phot_lambda.size(); i++) {\
            if( !qcal_opdet_is_detected[i] ) continue;\
            v.push_back( qcal_opdet_phot_lambda[i] );\
        }\
        return v;")


    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )
    hx1 = rt.RDF.TH1DModel( ut.prepare_TH1D("hx1", xbin, xmin, xmax) )
    hx = df.Histo1D(hx, "qcal_opdet_phot_lambda")
    hx1 = df.Histo1D(hx1, "qcal_opdet_phot_lambda_det")

    hx = hx.GetValue()
    hx1 = hx1.GetValue()

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig, ax = plt.subplots()
    fig.set_size_inches(5, 5)

    set_axes_color(ax, col)

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Counts")

    leg = legend()
    leg.add_entry(leg_lin("blue"), "Photons incident on sensors")
    leg.add_entry(leg_lin("red"), "Detected photons")
    leg.draw(plt, col)

    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    hx = hep.histplot(hx, color="blue")
    hx1 = hep.histplot(hx1, color="red")

    ax.set_yscale("log")

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#main

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
if __name__ == "__main__":

    gROOT.SetBatch()

    #rt.EnableImplicitMT()

    #photoelectron functions
    photoelectrons()

    main()



