#!/usr/bin/env python3

#photoelectron distribution in time

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, gSystem, gInterpreter
from ROOT import RDataFrame, TPad, TLatex

from photoelectrons import photoelectrons

import sys
sys.path.append("../..")
import plot_utils as ut

#_____________________________________________________________________________
@rt.Numba.Declare(['const RVec<bool>'], "bool")
def my_func(i):
    print("hi from my_func")
    return 1

#_____________________________________________________________________________
def main():

    inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal3cx2/en_9/lmon.root"

    iev = 57; # 57  9  5  4
    xbin = 0.06
    xmin = 0.5
    xmax = 2

    #col = rt.kCyan
    col = rt.kBlue

    df = RDataFrame("DetectorTree", inp)
    df = df.Range(iev, iev+1)

    #print(df.GetColumnNames())

    #detected flag
    df = df.Define("qcal_opdet_is_detected", "sens(qcal_opdet_phot_en)")
    df = df.Define("time_xxx", "senstime(qcal_opdet_cell_id, qcal_opdet_time, qcal_opdet_is_detected)")
    #df = df.Define("time_xxx", "Numba::my_func(qcal_opdet_is_detected)")

    #dummy histogram to call sensor_time senstime call operator 
    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", 1, 0, 1) )
    hx = df.Histo1D(hx, "time_xxx").GetValue()

    #sensor ID and photoelectron times
    sens_time = rt.senstime.GetSens()

    #photons counts for sensor ID
    sensId_count = {}
    for i in sens_time:
        sensId_count[i[0]] = i[1].size()
    #sorted by photoelectron counts
    sensId_count_sort = dict(sorted(sensId_count.items(), key=lambda item: item[1], reverse=True))

    #time distribution for individual sensors
    hist = []
    for i in sensId_count_sort.keys():

        #time histogram for sensor i
        hist.append( ut.prepare_TH1D(str(i), xbin, xmin, xmax) )
        h = hist[len(hist)-1]
        h.SetYTitle("")

        #photoelectron loop
        for j in range(sens_time[i].size()):
            h.Fill( sens_time[i].at(j) )

        ut.line_h1(h, col)
        h.SetFillColor(col)

        #normalize by bin width
        #h.Scale(1./xbin)

    plot_mat(hist)
    return

    can = ut.box_canvas()

    ih = 0

    hist[ih].Draw("")

    print(hist[ih].GetEntries(),hist[ih].GetBinContent(0), hist[ih].GetBinContent(hist[ih].GetNbinsX()+1))

    gPad.SetGrid()

    #ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

    #awkward array holding photoelectron times for sensor IDs

    #for i,j in sensId_count_sort.items():
        #if j < 5: continue
        #print(i,j)
        #print(type(sens_time[i]))
        #pass

#main

#_____________________________________________________________________________
def plot_mat(hist):

    #several sensors as a matrix
    can = ut.box_canvas(800, 600)

    nx = 3
    ny = 2

    ytit = "Fired microcells"
    yofs = 1.5

    xtit = "Time (ns)"
    xofs = 1.2

    lmgg = 0.21
    rmgg = 0.03
    bmgg = 0.18
    tmgg = 0.03

    #lmgg = 0
    #rmgg = 0
    #bmgg = 0
    #tmgg = 0

    tsiz = 0.07

    #vertical maximum
    ymax = -1
    for i in hist:
        if i.GetMaximum() > ymax: ymax = i.GetMaximum()
    ymax = 1.1*ymax

    ijmap = {}
    ii = 1
    dx = 1./nx
    dy = 1./ny
    pads = []
    for i in range(ny):
        for j in range(nx):

            ijmap[ii] = (i, j)

            #print(ii, i, j)

            pnam = "pad_"+str(ii)
            pads.append(TPad(pnam, pnam, j*dx, (ny-i-1)*dy, dx*(j+1), (ny-i)*dy)) #xlo, ylo, xhi, yhi
            pads[len(pads)-1].Draw()

            ii += 1

    for i in range(1,nx*ny+1):

        pads[i-1].cd()

        hx = hist[i-1]

        hx.SetTitle("")

        ut.set_H1_text_size(hx, tsiz)
        hx.SetTitleSize(tsiz, "Z")
        hx.SetLabelSize(tsiz, "Z")

        hx.SetMinimum(0)
        hx.SetMaximum(ymax)

        lmg = lmgg
        bmg = bmgg
        tmg = tmgg
        rmg = rmgg

        ip = ijmap[i][1]
        jp = ijmap[i][0]

        if ip < nx-1:
            hx.GetZaxis().SetTickLength(0)

        if ip == 0 and jp == 0:
            hx.SetYTitle(ytit)
            hx.SetTitleOffset(yofs, "Y")
            hx.GetYaxis().ChangeLabel(1, -1, 0)

        if jp == ny-1 and ip == nx-1:
            hx.SetXTitle(xtit)
            hx.SetTitleOffset(xofs, "X")

        if jp == ny-1 and ip != nx-1:
            hx.GetXaxis().ChangeLabel(-1, -1, 0)

        if ip == 0: rmg = 0

        if ip > 0 and ip < nx-1:
            lmg = 0
            rmg = 0

        if ip == nx-1: lmg = 0

        if jp == 0: bmg = 0
        if jp == ny-1: tmg = 0
        if jp > 0 and jp < ny-1:
            bmg = 0
            tmg = 0

        ut.set_margin_lbtr(gPad, lmg, bmg, tmg, rmg)

        gPad.SetGrid()

        hx.Draw("")

        #label with sensor ID
        lab = TLatex()
        lab.SetTextSize(tsiz)
        xmin = hx.GetXaxis().GetXmin()
        xmax = hx.GetXaxis().GetXmax()
        lab.DrawLatex(xmin+0.6*(xmax-xmin), 0.85*ymax, "Sens# "+hx.GetName())

        ut.invert_col(gPad)

    can.SaveAs("01fig.pdf")

#plot_mat

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetLineWidth(2)
    gStyle.SetFrameLineWidth(2)
    gStyle.SetOptStat("")

    #rt.EnableImplicitMT()

    #photoelectron functions
    photoelectrons()

    #call to main
    main()














