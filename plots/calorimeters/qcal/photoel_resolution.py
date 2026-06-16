#!/usr/bin/env python3

#yields and resolution in terms of photoelectrons

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
import numpy as np

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, RDataFrame, TF1

from photoelectrons import photoelectrons

import sys
sys.path.append("../..")
import plot_utils as ut

#_____________________________________________________________________________
def main():


    iplot = 1

    func = {}
    func[0] = yields
    func[1] = resolution

    func[iplot]()

#main

#_____________________________________________________________________________
def yields(draw=True):

    #inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal3cx3/en_","/lmon.root"]
    inp = ["/home/jaroslav/sim/lmon2-data/qcal/qcal3cx4/en_","/lmon.root"]

    energy = [1, 5, 9, 14, 18]

    #photoelectron counts
    xmin = 0
    #xmax = 240
    #xmax = 700
    xmax = 2800
    #xbin = 2
    xbin = 20

    hx = ut.prepare_TH1D("hx", xbin, xmin, xmax)

    print(inp[0])

    func = []
    hist = []
    for i in energy:

        print()
        print("Energy:", i)

        df = RDataFrame("DetectorTree", inp[0]+str(i)+inp[1])
        df = df.Define("qcal_opdet_is_detected", "sens(qcal_opdet_phot_en)") # detected flag
        df = df.Define("qcal_nphotoel", "nphotoel_evt(qcal_opdet_is_detected)") # photoelectron counts

        h_en = rt.RDF.TH1DModel( hx.Clone("hx_"+str(i)) )
        h_en = df.Histo1D(h_en, "qcal_nphotoel").GetValue()

        #normalize to number of events
        h_en.Scale( 1./(df.Count().GetValue()) )
        print("Integral:", h_en.Integral())

        fg = TF1("fg_"+str(i), "gaus", xmin, xmax)
        fg.SetParameters(1, h_en.GetMean(), h_en.GetStdDev())
        fg.SetNpx(2000)

        h_en.Fit(fg)

        func.append(fg)
        hist.append(h_en)

    if not draw: return func

    can = ut.box_canvas()

    hx = hist[0]

    for i in hist: ut.line_h1(i, rt.kBlue, 2)

    hx.Draw()

    for i in hist: i.Draw("same")

    ut.put_yx_tit(hx, "Counts / event", "Fired microcells in event", 1.7, 1.2)

    ut.set_margin_lbtr(gPad, 0.12, 0.1, 0.03, 0.03)

    for i in func: i.Draw("same")

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#yields

#_____________________________________________________________________________
def resolution():

    energy = [1, 5, 9, 14, 18]

    func = yields(False)
    res = [func[i].GetParameter(2)/func[i].GetParameter(1) for i in range(len(func))]

    print(res)

    plt.style.use("dark_background")
    col = "lime"
    #col = "black"

    #fit the resolution
    pars, cov = curve_fit(resf2, energy, res)

    fig = plt.figure()
    fig.set_size_inches(5, 5)
    ax = fig.add_subplot(1, 1, 1)

    set_axes_color(ax, col)

    plt.grid(True, color = col, linewidth = 0.5, linestyle = "--")

    #resolution data
    plt.plot(energy, res, marker="o", linestyle="", color="blue")

    #plot the fit function
    x = np.linspace(energy[0], energy[-1], 300)
    y = resf2(x, pars[0], pars[1])

    plt.plot(x, y, "-", label="resf", color="red")

    plt.rc("text", usetex = True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
    ax.set_xlabel("Generated gamma energy $E_\gamma$ (GeV)")
    ax.set_ylabel(r"Resolution $\sigma_N/\langle N\rangle$")

    #fit parameters for legend
    fit_param = ""
    fit_param += r"\begin{align*}"
    fit_param += r"a &= {0:.4f} \pm {1:.4f}\\".format(pars[0], np.sqrt(cov[0,0]))
    fit_param += r"b &= {0:.4f} \pm {1:.4f}".format(pars[1], np.sqrt(cov[1,1]))
    fit_param += r"\end{align*}"

    leg = legend()
    leg.add_entry(leg_lin("red"), r"$\frac{\sigma_N}{\langle N\rangle} = \frac{a}{\sqrt{E}} \oplus\ b$")
    leg.add_entry(leg_txt(), fit_param)
    leg.draw(plt, col)

    plt.savefig("01fig.pdf", bbox_inches = "tight")

#resolution

#_____________________________________________________________________________
def resf2(E, a, b):

    #resolution function  sigma/E = sqrt( a^2/E + b^2 )
    return np.sqrt( (a**2)/E + b**2 )

#resf2

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
    gStyle.SetFrameLineWidth(2)

    #rt.EnableImplicitMT()

    #photoelectron functions
    photoelectrons()

    #call to main
    main()

