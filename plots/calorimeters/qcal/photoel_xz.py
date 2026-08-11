#!/usr/bin/env python3

#cell photoelectron counts in x and z

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, gSystem, gInterpreter
from ROOT import RDataFrame, TMath

from photoelectrons import photoelectrons

import sys
sys.path.append("../..")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    #inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal3cx3/en_5/lmon.root"
    #inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal3ex3/en_1/lmon.root"
    #inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal3fx1/en_9/lmon.root"
    #inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal3gx3/en_18/lmon.root" # geom_fib_1mm_dx7mm.xml
    inp = "/home/jaroslav/sim/lmon2-data/qcal/qcal3hx3/en_9/lmon.root" # geom_fib_1mm_dx6mm.xml

    #geometry parameters
    #nx = 15
    #nz = 20
    #cell_xy = 12.2
    #cell_z = 263.571645
    #modz = 522.814704

    # geom_fib_1mm_dx10mm.xml:
    #nx = 15; nz = 20; cell_xy = 31; cell_z = 718.533188; modz = 1362.971788

    # geom_fib_1mm_dx8mm.xml:
    #nx = 21; nz = 22; cell_xy = 25; cell_z = 795.243866; modz = 1322.462120

    # geom_fib_1mm_dx7mm.xml
    #nx = 25; nz = 25; cell_xy = 21; cell_z = 799.043866; modz = 1299.551861

    # geom_fib_1mm_dx6mm.xml
    nx = 17; nz = 26; cell_xy = 18.1; cell_z = 498.087950; modz = 1004.931637

    cell_phi = 45*TMath.Pi()/180
    modx = cell_xy*nx

    #object of cell_pos_xz defined in photoelectrons.py
    par = str(nx)+", "+str(nz)+", "+str(cell_xy)+", "+str(cell_z)+", "+str(cell_phi)+", "+str(modx)+", "+str(modz)
    gInterpreter.Declare("cell_pos_xz pos_xz("+par+");");

    #bins in x and z
    xbin = cell_xy
    zstart = rt.pos_xz.xzpos.at(0).second
    zend = rt.pos_xz.xzpos.at(nx*nz-1).second
    zbin = (zstart-zend)/(nz-1)

    #range in x and z
    xstart = rt.pos_xz.xzpos.at(0).first
    xmin = xstart - 0.5*xbin
    xmax = xmin + nx*xbin
    zmin = zend - 0.5*zbin
    zmax = zmin + nz*zbin

    print(xbin, xmin, xmax)
    print(zbin, zmin, zmax)

    df = RDataFrame("DetectorTree", inp)
    #df = df.Range(1)

    #detected flag
    df = df.Define("qcal_opdet_is_detected", "sens(qcal_opdet_phot_en)")

    #photoelectrons x and z position
    df = df.Define("qcal_opdet_xpos_det", "pos_xz(qcal_opdet_cell_id, qcal_opdet_is_detected)")
    df = df.Define("qcal_opdet_zpos_det", "pos_xz(qcal_opdet_cell_id, qcal_opdet_is_detected, 1)")

    hx = rt.RDF.TH1DModel( ut.prepare_TH1D_n("hx", nx, xmin, xmax) )
    hz = rt.RDF.TH1DModel( ut.prepare_TH1D_n("hz", nz, zmin, zmax) )
    #hx = df.Histo1D(hx, "qcal_opdet_xpos_det").GetValue()
    #hz = df.Histo1D(hz, "qcal_opdet_zpos_det").GetValue()
    hxz = df.Histo2D(("hxz","", nx, xmin, xmax, nz, zmin, zmax), "qcal_opdet_xpos_det", "qcal_opdet_zpos_det").GetValue()
    #ut.line_h1(hx, rt.kBlue, 2)
    #ut.line_h1(hz, rt.kBlue, 2)

    #print(hx.GetBinWidth(1))
    #print(hz.GetBinWidth(1))

    print(hxz.GetXaxis().GetBinWidth(1))
    print(hxz.GetYaxis().GetBinWidth(1)) 

    #number of events processed
    nev = df.Count().GetValue()

    #normalize by number of events
    hxz.Scale(1./nev)
    print(hxz.Integral())

    hxz.SetMinimum(0.98/nev)

    can = ut.box_canvas() # 800, 600

    #hx.Draw("")
    #hz.Draw("")
    hxz.Draw("colz")

    ut.put_yx_tit(hxz, "#it{z} (mm)", "#it{x} (mm)", 1.6, 1.2)
    hxz.SetZTitle("Fired microcells / event")
    hxz.SetTitleOffset(1.4, "Z")

    ut.set_margin_lbtr(gPad, 0.11, 0.09, 0.02, 0.15)

    gPad.SetGrid()

    gPad.SetLogz()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

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


