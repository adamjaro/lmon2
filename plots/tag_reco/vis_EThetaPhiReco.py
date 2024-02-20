#!/usr/bin/python3

# visualization for EThetaPhiReco

import code
from ctypes import c_uint64, c_double

import ROOT as rt
from ROOT import gSystem, TFile, TIter, TMath
from ROOT import TGeoManager, TGeoMaterial, TGeoMedium, gGeoManager
from ROOT import TGeoRotation, TGeoTranslation, TGeoCombiTrans
from ROOT import TEveManager, TEveGeoTopNode, gEve, TEveArrow

#_____________________________________________________________________________
def main():

    #response input
    inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/tag_resp.root"
    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx2/tag_resp_pass2.root"

    #layer index
    ilay = 0

    #maximal number of points to show
    npoints = 2500

    #open the input
    infile = TFile.Open(inp)
    tree = infile.Get("s1_links_"+str(ilay))
    quant = infile.Get("s1_quantities_"+str(ilay))

    #maximal index number for visualization scale
    nmax_idx = c_uint64(1)

    qnext = TIter(quant)
    obj = qnext()
    while obj != None:
        nmax_idx.value *= obj.GetNbinsX()+1
        obj = qnext()

    print("Maximal index:", nmax_idx.value)


    idx = c_uint64(0)
    en = c_double(0)
    theta = c_double(0)
    phi = c_double(0)
    tree.SetBranchAddress("idx", idx)
    tree.SetBranchAddress("en", en)
    tree.SetBranchAddress("theta", theta)
    tree.SetBranchAddress("phi", phi)

    #top geometry
    geom = TGeoManager("geom", "geom")

    #vacuum material, done via TGeoManager to prevent segmentation violation
    geom.Material("Vacuum", 0, 0, 0, 0)
    med = geom.Medium("Vacuum", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    #top volume
    top = geom.MakeBox("TOP", med, 1000., 1000., 1000.);
    geom.SetTopVolume(top)

    #central column to represent the indices
    idx_col = geom.MakeTube("IDX", med, 0, 1, 50)
    idx_col.SetLineColor(rt.kRed)
    idx_col.SetTransparency(0)
    top.AddNode(idx_col, 1, TGeoCombiTrans("idx_trans", 0, 50, 0, TGeoRotation("idx_rot", 0, 90, 0)))


    print("Input entries:", tree.GetEntries())

    mod = 1
    if npoints < tree.GetEntries():
        mod = TMath.Ceil(tree.GetEntries()/npoints)
        print("mod:", mod)

    #loop over the tree
    ipoint = 2
    arrows = []
    for ient in range(tree.GetEntries()):

        if ient % mod != 0:
            continue

        tree.GetEntry(ient)

        #point on electron kinematics
        point = geom.MakeSphere("p_"+str(ipoint), med, 0, 1)
        point.SetLineColor(rt.kYellow)
        point.SetTransparency(0)
        #point_y = 100*(TMath.Pi()-theta.value)*1e3/16 # mrad, bottom-up
        point_y = 100-100*(TMath.Pi()-theta.value)*1e3/16 # mrad, top-down
        point_r = 100*(en.value-2)/(18-2) # radial distance by energy, GeV
        point_x = point_r*TMath.Sin(TMath.Pi()+phi.value) # xz plane by phi, rad
        point_z = point_r*TMath.Cos(TMath.Pi()+phi.value)
        top.AddNode(point, ipoint, TGeoTranslation(point_x, point_y, point_z))
        ipoint += 1

        #line from central index to electron kinematics
        line_start = 100*idx.value/nmax_idx.value
        arrows.append( TEveArrow(point_x, point_y-line_start, point_z, 0, line_start, 0) )

    print("Number of points shown:", len(arrows))

    #end of geometry definition
    geom.CloseGeometry()

    #add geometry to TEve
    TEveManager.Create()
    node = gGeoManager.GetTopNode()
    topnode = TEveGeoTopNode(gGeoManager, node)
    topnode.SetVisLevel(4)
    topnode.GetNode().GetVolume().SetVisibility(rt.kFALSE)
    gEve.AddGlobalElement(topnode)

    #put the connecting lines represented by arrows
    for i in arrows:
        i.SetMainColor(rt.kOrange)
        i.SetTubeR(0.003)
        i.SetConeR(0.003)
        i.SetConeL(0)
        gEve.AddGlobalElement(i)

    #end of objects definition
    gEve.FullRedraw3D(rt.kTRUE)
    viewer = gEve.GetDefaultGLViewer()

    viewer.SetClearColor(rt.kBlack)
    #viewer.SetClearColor(rt.kWhite)

    code.interact(local=locals())

#main

#_____________________________________________________________________________
if __name__ == "__main__":

    gSystem.Load("libGeom")

    main()

