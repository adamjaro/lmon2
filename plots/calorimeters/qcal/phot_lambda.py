#!/usr/bin/env python3

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, RDataFrame

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import mplhep as hep

import sys
sys.path.append("../..")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal2bx1/en_9/lmon.root"

    df = RDataFrame("DetectorTree", inp)

    #wavelength from energy, E (eV) = 1240.637 / lambda (nm)
    df = df.Define("qcal_opdet_phot_lambda", "1240.637/qcal_opdet_phot_en")

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", 10, 150, 850) )
    hx = df.Histo1D(hx, "qcal_opdet_phot_lambda")

    hx = hx.GetValue()

    plt.style.use("dark_background")
    col = "lime"

    fig, ax = plt.subplots()
    fig.set_size_inches(5, 5)

    set_axes_color(ax, col)

    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    hx = hep.histplot(hx, color="red")

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
if __name__ == "__main__":

    gROOT.SetBatch()

    rt.EnableImplicitMT()

    main()



