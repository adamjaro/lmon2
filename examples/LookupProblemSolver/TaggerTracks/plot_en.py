#!/usr/bin/python3

import matplotlib.pyplot as plt

import ROOT as rt
from ROOT import gPad, gROOT, gStyle
from ROOT import TFile

import sys
sys.path.append("../../../plots")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    infile = TFile.Open("rec_tracks.root")
    tree = infile.Get("rec_trk")

    xbin = 0.3
    xmin = 5
    xmax = 19

    can = ut.box_canvas()

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    tree.Draw("rec_en >> hx")

    gx = ut.h1_to_arrays(hx)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel(r"Reconstructed energy $E_{e}$ (GeV)")
    ax.set_ylabel(r"Counts")

    set_axes_color(ax, col)
    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    #plt.semilogy(gx[0], gx[1], "-", label="gx", color="red")
    plt.plot(gx[0], gx[1], "-", label="gx", color="red")



    #ytit = "Reconstructed energy #it{E_{e}} (GeV)"
    #xtit = "Generated MC particle energy #it{E_{e,mc}} (GeV)"
    #xtit = "Generated true electron energy #it{E_{e,true}} (GeV)"
    #ut.put_yx_tit(hxy, ytit, xtit, 1.5, 1.4)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

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
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetFrameLineWidth(2)

    main()

