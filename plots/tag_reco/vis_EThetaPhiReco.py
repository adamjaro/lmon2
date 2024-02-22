#!/usr/bin/python3

# visualization for EThetaPhiReco
# run as: 'python3 -i vis_EThetaPhiReco.py'

import code
from ctypes import c_uint64, c_double

import ROOT as rt
from ROOT import gSystem, TFile, TIter, TMath
from ROOT import TGeoManager, TGeoMaterial, TGeoMedium, gGeoManager
from ROOT import TGeoRotation, TGeoTranslation, TGeoCombiTrans
from ROOT import TEveManager, TEveGeoTopNode, gEve, TEveArrow, TEveText

#_____________________________________________________________________________
def main():

    #response input
    #inp = "/home/jaroslav/sim/lmon2/macro/low-Q2/tag_resp.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag7cx2/tag_resp_pass3.root"

    #layer index
    ilay = 0

    #maximal number of points to show
    npoints = 1000

    #energy range, GeV
    en_min = 2
    en_max = 18

    #theta range, mrad
    theta_max = 16

    #open the input
    infile = TFile.Open(inp)
    tree = infile.Get("s1_links_"+str(ilay))
    quant = infile.Get("s1_quantities_"+str(ilay))

    #maximal index number for visualization scale
    nmax_idx = c_uint64(1)

    #determine the maximal index from quantity histograms
    qnext = TIter(quant)
    obj = qnext()
    while obj != None:
        nmax_idx.value *= obj.GetNbinsX()+1
        obj = qnext()

    print("Maximal index:", nmax_idx.value)

    #connect the links
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

    #division for links selected to be drawn
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

        #get the link
        tree.GetEntry(ient)

        #point on electron kinematics
        point = geom.MakeSphere("p_"+str(ipoint), med, 0, 0.8)
        point.SetLineColor(rt.kGreen) # rt.kYellow
        point.SetTransparency(0)
        point_y = 100-100*(TMath.Pi()-theta.value)*1e3/theta_max # mrad, top-down
        point_r = 100*(en.value-en_min)/(en_max-en_min) # radial distance by energy, GeV
        point_x = point_r*TMath.Sin(TMath.Pi()+phi.value) # xz plane by phi, rad
        point_z = point_r*TMath.Cos(TMath.Pi()+phi.value)
        top.AddNode(point, ipoint, TGeoTranslation(point_x, point_y, point_z))
        ipoint += 1

        #line from central index to electron kinematics
        line_start = 100*idx.value/nmax_idx.value # central index, bottom-up
        arrows.append( TEveArrow(point_x, point_y-line_start, point_z, 0, line_start, 0) )

    print("Number of points shown:", len(arrows))

    #coordinate frame

    #rings representing radial coordinate (energy)
    nrings = 4
    for i in range(nrings):
        ring_r = (i+1)*100/nrings
        ring = geom.MakeTube("ring_"+str(ipoint), med, ring_r, ring_r+1, 0.5)
        ring.SetLineColor(rt.kBlue)
        ring.SetTransparency(0)
        top.AddNode(ring, ipoint, TGeoRotation("idx_rot", 0, 90, 0))
        ipoint += 1
        #top.AddNode(ring, ipoint, TGeoCombiTrans("ring_t", 0, 100, 0, TGeoRotation("ring_r", 0, 90, 0)))
        #ipoint += 1

    #rays representing azimuthal angles
    nrays = 8
    for i in range(nrays):
        ray = geom.MakeBox("ray_"+str(ipoint), med, 0.5, 0.5, 53)
        ray.SetLineColor(rt.kBlue)
        ray.SetTransparency(0)
        ray_phi = i*2*TMath.Pi()/nrays
        ray_x = 53*TMath.Sin(ray_phi)
        ray_z = 53*TMath.Cos(ray_phi)
        top.AddNode(ray, ipoint, TGeoCombiTrans("r_tr", ray_x, 0, ray_z, TGeoRotation("r_rot", 90, ray_phi*180/TMath.Pi(), 0)))
        ipoint += 1

    #vertical bar representing polar angles
    vbar = geom.MakeBox("vbar", med, 0.5, 50, 0.5)
    vbar.SetLineColor(rt.kBlue)
    vbar.SetTransparency(0)
    top.AddNode(vbar, ipoint, TGeoTranslation(0, 50, 100.5))
    ipoint += 1

    #tick marks for theta angles
    nticks = 5
    for i in range(nticks):
        vbar_tick = geom.MakeBox("vbar_tick", med, 0.5, 0.5, 3)
        vbar_tick.SetLineColor(rt.kBlue)
        vbar_tick.SetTransparency(0)
        top.AddNode(vbar_tick, ipoint, TGeoTranslation(0, 100-i*100/nticks, 100+3))
        ipoint += 1

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

    font_size = 20

    #tick labels for theta angles
    vbar_txt = TEveText("theta (mrad)")
    vbar_txt.SetFontSize( font_size )
    vbar_txt.SetMainColor( rt.kRed )
    vbar_txt.RefMainTrans().SetPos(0, 110, 100)
    gEve.AddGlobalElement(vbar_txt)
    for i in range(nticks):
        vbar_txt = TEveText( "{0:.1f}".format(i*theta_max/nticks) )
        vbar_txt.SetFontSize( font_size )
        vbar_txt.SetMainColor( rt.kRed )
        vbar_txt.RefMainTrans().SetPos(0, 100-i*100/nticks, 110)
        gEve.AddGlobalElement(vbar_txt)

    #tick labels for polar angle rays
    ray_txt = TEveText("phi (deg)")
    ray_txt.SetFontSize( font_size )
    ray_txt.SetMainColor( rt.kRed )
    ray_txt.RefMainTrans().SetPos(0, -4, 110)
    gEve.AddGlobalElement(ray_txt)
    for i in range(nrays):
        ray_txt = TEveText( "{0:.1f}".format(i*360/nrays) )
        ray_txt.SetFontSize( font_size )
        ray_txt.SetMainColor( rt.kRed )
        ray_phi = i*2*TMath.Pi()/nrays
        ray_x = 110*TMath.Sin(ray_phi)
        ray_z = 110*TMath.Cos(ray_phi)
        ray_txt.RefMainTrans().SetPos(ray_x, 3, ray_z)
        gEve.AddGlobalElement(ray_txt)

    #ring labels for energy
    en_txt = TEveText("Energy (GeV)")
    en_txt.SetFontSize( font_size )
    en_txt.SetMainColor( rt.kRed )
    #en_txt.RefMainTrans().SetPos(70, 10, 0)
    en_txt.RefMainTrans().SetPos(0, 10, 50)
    gEve.AddGlobalElement(en_txt)
    for i in range(nrings):
        en_txt = TEveText( "{0:.1f}".format(en_min+(i+1)*(en_max-en_min)/nrings) )
        en_txt.SetFontSize( font_size )
        en_txt.SetMainColor( rt.kRed )
        #en_txt.RefMainTrans().SetPos(((i+1)*100/nrings)-5, 3, 0)
        en_txt.RefMainTrans().SetPos(-5, 3, ((i+1)*100/nrings))
        gEve.AddGlobalElement(en_txt)


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

