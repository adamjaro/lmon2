#!/usr/bin/python3

# beam phase space, generated and at a given location

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, gSystem, TMath
from ROOT import RDataFrame

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    iplot = 3

    func = {}
    func[0] = generated_horizontal
    func[1] = gen_theta_x
    func[2] = Q3ER_xy
    func[3] = Q3ER_vertical
    func[4] = Q3ER_horizontal
    func[iplot]()

#main

#_____________________________________________________________________________
def generated_horizontal():

    inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/lmon.root"

    #mm
    xbin = 3
    xmin = -600
    xmax = 600

    #urad
    ybin = 3
    ymin = -800
    ymax = 800

    df = RDataFrame("DetectorTree", inp)

    #micro rad, mm
    df = df.Define("mcp0_theta_x", "1e6*TMath::ATan( TMath::Cos(mcp_phi[0])*TMath::Tan(mcp_theta[0]) )")
    df = df.Define("mcp0_vx", "1e3*mcp_vx[0]")

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, ybin, ymin, ymax) )

    hx = df.Histo2D(hx, "mcp0_vx", "mcp0_theta_x").GetValue()

    hvx = hx.ProjectionX()
    print("Projection along x:")
    print("mean:", hvx.GetMean(), "+/-", hvx.GetMeanError())
    print("sigma:", hvx.GetStdDev(), "+/-", hvx.GetStdDevError())

    htx = hx.ProjectionY()
    print("Projection along y:")
    print("mean:", htx.GetMean(), "+/-", htx.GetMeanError())
    print("sigma:", htx.GetStdDev(), "+/-", htx.GetStdDevError())

    hx.SetMinimum(0.8)

    hx.Draw("colz")

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#generated_horizontal

#_____________________________________________________________________________
def gen_theta_x():

    inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/lmon.root"

    xbin = 3
    xmin = -1200
    xmax = 1200

    df = RDataFrame("DetectorTree", inp)

    #micro rad
    #df = df.Define("mcp0_theta_x", "1e6*TMath::ATan( TMath::Sin(mcp_theta[0])*TMath::Cos(mcp_phi[0])/TMath::Cos(mcp_theta[0]) )")
    #df = df.Define("mcp0_theta_y", "1e6*TMath::ATan( TMath::Sin(mcp_theta[0])*TMath::Sin(mcp_phi[0])/TMath::Cos(mcp_theta[0]) )")
    df = df.Define("mcp0_theta_x", "1e6*TMath::ATan( TMath::Cos(mcp_phi[0])*TMath::Tan(mcp_theta[0]) )")
    df = df.Define("mcp0_theta_y", "1e6*TMath::ATan( TMath::Sin(mcp_phi[0])*TMath::Tan(mcp_theta[0]) )")

    can = ut.box_canvas()
    hx = rt.RDF.TH1DModel( ut.prepare_TH1D("hx", xbin, xmin, xmax) )

    #hx = df.Histo1D(hx, "mcp0_theta_x").GetValue()
    hx = df.Histo1D(hx, "mcp0_theta_y").GetValue()

    hx.Draw()

    print(hx.GetStdDev(), "+/-", hx.GetStdDevError())

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#gen_theta_x

#_____________________________________________________________________________
def Q3ER_xy(get_df=False):

    #inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/lmon.root"
    #inp = "/home/jaroslav/sim/lmon2-data/beam/v6p3_20240729/lmon.root"
    inp = "/home/jaroslav/sim/lmon2-data/beam/v6p3_20240802/lmon.root"


    #mm
    xbin = 0.1
    xmin = 5
    xmax = 15

    ymin = -5
    ymax = 5

    df = RDataFrame("DetectorTree", inp)

    #df = df.Range(12)

    gROOT.ProcessLine("Double_t zpos, xpos, theta, length;")
    #rt.zpos = -42760.222219 # Q3ER.zpos
    #rt.xpos = -465.263731 # Q3ER.xpos
    #rt.length = 600 # Q3ER.length
    rt.theta = 0.020000 # Q3ER.theta
    rt.zpos = -42460.282217
    rt.xpos = -459.264131
    rt.length = 0

    #Q3ER_SZ 1 -42460.3 -42460.282217
    #Q3ER_SX 1 -459.264 -459.264131

    gROOT.ProcessLine('auto to_local = [](const auto& zvec, const auto& xvec) {\
        vector<double> zxloc = {9999, 9999};\
        if(zvec.size() == 1) {\
          TVector2 v(zvec[0]-zpos-length/2, xvec[0]-xpos);\
          v = v.Rotate(-theta);\
          zxloc[0] = v.X();\
          zxloc[1] = v.Y();\
        }\
        /*cout << "hi " << zxloc[0] << " " << zxloc[1] << " " << xpos << endl;*/\
        return zxloc;}')

    df = df.Define("Q3ER_counter_zxloc", "to_local(Q3ER_counter_z, Q3ER_counter_x)")
    df = df.Define("Q3ER_counter_xloc", "Q3ER_counter_zxloc[1]").Filter("Q3ER_counter_xloc<9998")
    df = df.Define("Q3ER_counter_yloc", "Q3ER_counter_y[0]")

    if get_df is True:
        return df

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, xbin, ymin, ymax) )
    hx = df.Histo2D(hx, "Q3ER_counter_xloc", "Q3ER_counter_yloc").GetValue()

    #print(df.GetColumnNames())
    #print(df.GetColumnType("Q3ER_counter_x"))
    #print(df.GetColumnType("Q3ER_counter_HitZ"))

    hvx = hx.ProjectionX()
    print("Projection along x:")
    print("mean:", hvx.GetMean(), "+/-", hvx.GetMeanError())
    print("sigma:", hvx.GetStdDev(), "+/-", hvx.GetStdDevError())

    htx = hx.ProjectionY()
    print("Projection along y:")
    print("mean:", htx.GetMean(), "+/-", htx.GetMeanError())
    print("sigma:", htx.GetStdDev(), "+/-", htx.GetStdDevError())

    hx.SetMinimum(0.8)

    hx.Draw("colz")

    gPad.SetGrid()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#Q3ER_xy

#_____________________________________________________________________________
def Q3ER_vertical():

    #mm
    xbin = 0.02
    xmin = -2
    xmax = 2

    #urad
    ybin = 0.3
    ymin = -70
    ymax = 70

    df = Q3ER_xy(True)
    df = df.Define("Q3ER_theta_y", "1e6*TMath::ATan(Q3ER_counter_ydir[0]/Q3ER_counter_zdir[0])") # urad

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, ybin, ymin, ymax) )
    hx = df.Histo2D(hx, "Q3ER_counter_yloc", "Q3ER_theta_y").GetValue()

    ut.put_yx_tit(hx, "#it{#theta}_{y} (#murad)", "#it{y} (mm)", 1.7, 1.3)

    ut.set_margin_lbtr(gPad, 0.12, 0.1, 0.03, 0.11)

    ut.print_projections_xy(hx)

    hx.SetMinimum(0.8)
    hx.SetContour(300)

    hx.Draw("colz")

    gPad.SetGrid()

    gPad.SetLogz()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#Q3ER_vertical

#_____________________________________________________________________________
def Q3ER_horizontal():

    #mm
    xbin = 0.1
    #xmin = 3
    #xmax = 17
    xmin = -8
    xmax = 8

    #urad
    ybin = 1
    #ymin = -470
    #ymax = 200
    ymin = -350
    ymax = 350

    df = Q3ER_xy(True)
    df = df.Define("Q3ER_theta_x", "1e6*TMath::ATan(Q3ER_counter_xdir[0]/Q3ER_counter_zdir[0])-1e6*0.02") # urad

    can = ut.box_canvas()
    hx = rt.RDF.TH2DModel( ut.prepare_TH2D("hx", xbin, xmin, xmax, ybin, ymin, ymax) )
    hx = df.Histo2D(hx, "Q3ER_counter_xloc", "Q3ER_theta_x").GetValue()

    ut.put_yx_tit(hx, "#it{#theta}_{x} (#murad)", "#it{x} (mm)", 1.7, 1.3)

    ut.set_margin_lbtr(gPad, 0.12, 0.1, 0.03, 0.11)

    ut.print_projections_xy(hx)

    hx.SetMinimum(0.8)
    hx.SetContour(300)

    hx.Draw("colz")

    gPad.SetGrid()

    gPad.SetLogz()

    ut.invert_col(rt.gPad)
    can.SaveAs("01fig.pdf")

#Q3ER_horizontal

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetFrameLineWidth(2)

    main()











