#!/usr/bin/python3

# beam-gas properties

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pandas import DataFrame
from scipy.interpolate import interp1d
import numpy as np

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import RDataFrame

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 1

    func = {}
    func[0] = beta_xy
    func[1] = sigma_xy
    func[2] = divergence_xy

    func[iplot]()

#main

#_____________________________________________________________________________
def beta_xy():

    df = load_lattice()

    #print(lat["S"])
    print(df)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    set_axes_color(ax, col)
    set_grid(plt, col)

    plt.plot(df["S"], df["BETX"], "o", markersize=3, color="blue", lw=1)
    plt.plot(df["S"], df["BETY"], "o", markersize=3, color="red", lw=1)

    beta_ix = interp1d(df["S"], df["BETX"], kind="linear")
    beta_iy = interp1d(df["S"], df["BETY"], kind="linear")

    xs = np.linspace(df["S"][ df["S"].index[0] ], df["S"][ df["S"].index[-1] ], 300)

    plt.plot(xs, [beta_ix(i) for i in xs], "-", color="blue", lw=1)
    plt.plot(xs, [beta_iy(i) for i in xs], "--", color="red", lw=1)

    fig.savefig("01fig.pdf", bbox_inches = "tight")
    plt.close()

#beta_xy

#_____________________________________________________________________________
def sigma_xy():

    df = load_lattice()

    df = df.query("S>-44 and S<32")

    print(df)

    #emittance, meters
    eps_x = 24e-9 # m
    eps_y = 2e-9 # m

    #number of sigmas to show
    nsig_x = 15
    nsig_y = 23

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(6, 5)
    ax = fig.add_subplot(1, 1, 1)
    set_axes_color(ax, col)
    set_grid(plt, col)

    plt.plot(df["S"], beam_sigma(eps_x, df["BETX"], nsig_x), "o", markersize=3, color="blue", lw=1)
    plt.plot(df["S"], beam_sigma(eps_y, df["BETY"], nsig_y), "o", markersize=3, color="red", lw=1)

    sig_ix = interp1d(df["S"], beam_sigma(eps_x, df["BETX"], nsig_x), kind="linear")
    sig_iy = interp1d(df["S"], beam_sigma(eps_y, df["BETY"], nsig_y), kind="linear")

    xs = np.linspace(df["S"][ df["S"].index[0] ], df["S"][ df["S"].index[-1] ], 300)

    plt.plot(xs, [sig_ix(i) for i in xs], "-", color="blue", lw=1)
    plt.plot(xs, [sig_iy(i) for i in xs], "--", color="red", lw=1)

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("$s$ (meters)")
    ax.set_ylabel(r"Beam $\sigma$ (mm)")

    leg = legend()
    leg.add_entry(leg_lin("blue"), "15$\sigma_x$")
    leg.add_entry(leg_lin("red", "--"), "23$\sigma_y$")
    leg.add_entry(leg_txt(), "$\epsilon_{x/y}$ = 24/2 nm")
    leg.draw(plt, col)

    #print beam size at IP with BETX and Y at IP6 from the dataframe
    print("IP6 sigma_x (mm):", beam_sigma(eps_x, df.loc[ df["NAME"]=="IP6" ]["BETX"].values[0]))
    print("IP6 sigma_y (mm):", beam_sigma(eps_y, df.loc[ df["NAME"]=="IP6" ]["BETY"].values[0]))

    #print(df.loc[ df["NAME"]=="IP6" ]["BETY"].values[0])

    #front of Q3ER
    df_oww = df.loc[ df["NAME"]=="OWW_SH" ]
    #df_oww = df.loc[ df["NAME"]=="Q3ER_6" ]
    print("Q3 front sigma_x (mm):", beam_sigma(eps_x, df_oww["BETX"].values[0]))
    print("Q3 front sigma_y (mm):", beam_sigma(eps_y, df_oww["BETY"].values[0]))

    fig.savefig("01fig.pdf", bbox_inches = "tight")
    plt.close()

#sigma_xy

#_____________________________________________________________________________
def divergence_xy():

    df = load_lattice()

    df = df.query("S>-44 and S<32")

    print(df)

    #emittance, meters
    eps_x = 24e-9 # m
    eps_y = 2e-9 # m

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    fig = plt.figure()
    fig.set_size_inches(6, 5)
    ax = fig.add_subplot(1, 1, 1)
    set_axes_color(ax, col)
    set_grid(plt, col)

    plt.plot(df["S"], beam_divergence(eps_x, df["ALFX"], df["BETX"]), "o", markersize=3, color="blue", lw=1)
    plt.plot(df["S"], beam_divergence(eps_y, df["ALFY"], df["BETY"]), "o", markersize=3, color="red", lw=1)

    div_ix = interp1d(df["S"], beam_divergence(eps_x, df["ALFX"], df["BETX"]), kind="linear")
    div_iy = interp1d(df["S"], beam_divergence(eps_y, df["ALFY"], df["BETY"]), kind="linear")

    xs = np.linspace(df["S"][ df["S"].index[0] ], df["S"][ df["S"].index[-1] ], 300)

    plt.plot(xs, [div_ix(i) for i in xs], "-", color="blue", lw=1)
    plt.plot(xs, [div_iy(i) for i in xs], "--", color="red", lw=1)

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("$s$ (meters)")
    ax.set_ylabel(r"Angular divergence $\sigma_\theta$ ($\mathrm{\mu}$rad)")

    leg = legend()
    leg.add_entry(leg_lin("blue"), r"$\sigma_{\theta x}$")
    leg.add_entry(leg_lin("red", "--"), r"$\sigma_{\theta y}$")
    leg.add_entry(leg_txt(), "$\epsilon_{x/y}$ = 24/2 nm")
    leg.draw(plt, col)

    #print beam divergence at IP with ALFX and Y BETX and Y at IP6 from the dataframe
    df_ip6 = df.loc[ df["NAME"]=="IP6" ]
    print("IP6 div_x (urad):", beam_divergence(eps_x, df_ip6["ALFX"].values[0], df_ip6["BETX"].values[0]))
    print("IP6 div_y (urad):", beam_divergence(eps_y, df_ip6["ALFY"].values[0], df_ip6["BETY"].values[0]))

    #front of Q3ER
    df_oww = df.loc[ df["NAME"]=="OWW_SH" ]
    print("Q3ER front div_x (urad):", beam_divergence(eps_x, df_oww["ALFX"].values[0], df_oww["BETX"].values[0]))
    print("Q3ER front div_y (urad):", beam_divergence(eps_y, df_oww["ALFY"].values[0], df_oww["BETY"].values[0]))

    #print("IP6 sigma_y (mm):", beam_sigma(eps_y, df.loc[ df["NAME"]=="IP6" ]["BETY"].values[0]))

    #print(df_ip6)

    #print(df.loc[ df["NAME"]=="IP6" ]["BETY"].values[0])

    fig.savefig("01fig.pdf", bbox_inches = "tight")
    plt.close()

#divergence_xy

#_____________________________________________________________________________
def beam_sigma(eps, beta, mult=1):

    #beam sigma in mm

    return mult*1e3*np.sqrt(eps*beta)

#beam_sigma

#_____________________________________________________________________________
def beam_divergence(eps, alpha, beta):

    #angular divergence in micro rad

    return 1e6*np.sqrt( eps*( (1. + alpha**2)/beta ) )

#beam_divergence

#_____________________________________________________________________________
def load_lattice():

    infile = "/home/jaroslav/Knihy/Clanky/EIC/Magnets_IR/Lattice/Andrii_20240606/ir6twiss.tfs"

    inp = open(infile, "r")

    #skip the first 50 lines
    for i in range(50): inp.readline()

    #header line, skip the first '*' character
    head = inp.readline().split()[1:]

    #line with data types
    types = inp.readline().split()[1:]

    #read the data
    dat = []

    for i in inp:
        ii = i.split()
        lin = []

        for j in ii:

            #remove beginning and trailing '"' character
            if j[0] == '"': j = j[1:]
            if j[-1] == '"': j = j[:-1]

            lin.append( j )

        #string and float types
        for j in range(len(lin)):
            if types[j] == "%s":
                lin[j] = str(lin[j])
            else:
                lin[j] = float(lin[j])


        dat.append( lin )

    df = DataFrame(dat, columns=head)

    #subtract IP6 'S' position
    df["S"] = df["S"] - df.loc[ df["NAME"]=="IP6" ]["S"].values[0]

    #invert sign of S for detector frame
    df["S"] = -df["S"]

    return df

#load_lattice

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
def set_grid(px, col="lime"):

    px.grid(True, color = col, linewidth = 0.5, linestyle = "--")

#set_grid

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







