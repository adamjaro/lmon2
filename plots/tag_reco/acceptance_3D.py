#!/usr/bin/python3

# acceptance in 3D as a function of energy, theta and phi
# run as: 'python3 -i acceptance_3D.py'

import code
from ctypes import c_uint64, c_double

import ROOT as rt
from ROOT import gSystem, TFile, TIter, TMath, RDataFrame
from ROOT import TGeoManager, TGeoMaterial, TGeoMedium, gGeoManager
from ROOT import TGeoRotation, TGeoTranslation, TGeoCombiTrans
from ROOT import TEveManager, TEveGeoTopNode, gEve, TEveArrow, TEveText

import sys
sys.path.append("../")
import plot_utils as ut

#_____________________________________________________________________________
def main():

    #inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax3/trk_v1.root"
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9ax4/trk_v0.root"

    #tagger 1 or 2
    #tag = 2

    #energy range, GeV
    ebin = 1
    en_min = 2
    en_max = 18

    #theta range, mrad
    pitheta_bin = 1
    #pitheta_max = 16
    pitheta_max = 14

    #phi (rad)
    phi_bin = 0.12

    sel = "s1_ntrk>0"
    #sel = "s2_ntrk>0"

    #open the input
    df = RDataFrame("event", inp)
    df = df.Define("pitheta", "(TMath::Pi()-true_el_theta)*1e3")
    df = df.Define("phi_0_2pi", "true_el_phi+TMath::Pi()")

    hx = rt.RDF.TH3DModel( ut.prepare_TH3D("hx", ebin, en_min, en_max, pitheta_bin, 0, pitheta_max, phi_bin, 0, 2*TMath.Pi()) )

    hAll = df.Histo3D(hx, "true_el_E", "pitheta", "phi_0_2pi")
    hSel = df.Filter(sel).Histo3D(hx, "true_el_E", "pitheta", "phi_0_2pi")

    hAcc = hSel.GetValue()
    hAcc.Divide(hAll.GetValue())

    #top geometry
    geom = TGeoManager("geom", "geom")

    #vacuum material, done via TGeoManager to prevent segmentation violation
    geom.Material("Vacuum", 0, 0, 0, 0)
    med = geom.Medium("Vacuum", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    #top volume
    top = geom.MakeBox("TOP", med, 1000., 1000., 1000.);
    geom.SetTopVolume(top)

    ipoint = 1

    for ix in range(1, hAcc.GetNbinsX()+1):
        for iy in range(1, hAcc.GetNbinsY()+1):
            for iz in range(1, hAcc.GetNbinsZ()+1):

                ibin = hAcc.GetBin(ix, iy, iz)
                acc = hAcc.GetBinContent(ibin)

                #acc = 1-acc

                #if( acc < 1e-12 ): continue
                #if( hAcc.GetBinContent(ibin) < 0.1 ): continue

                #print(ix, iy, iz, ibin, hAcc.GetBinContent(ibin))

                #print(ix, hAcc.GetXaxis().GetBinLowEdge(ix))

                #e1, e2, t1, t2, p1, p2, tp
                e1 = hAcc.GetXaxis().GetBinLowEdge(ix)
                e2 = e1 + hAcc.GetXaxis().GetBinWidth(ix)

                t1 = hAcc.GetYaxis().GetBinLowEdge(iy)
                t2 = t1 + hAcc.GetYaxis().GetBinWidth(iy)

                p1 = hAcc.GetZaxis().GetBinLowEdge(iz)
                p2 = p1 + hAcc.GetZaxis().GetBinWidth(iz)

                tp = int( 100*(1-acc) )

                #tp = 0

                ipoint = draw_bin(e1, e2, t1, t2, p1, p2, tp, ipoint, en_min, en_max, pitheta_max, geom, top)

    #TH3D bin  E = 10 - 14 (14 - 18), phi = 45 - 90 (0 - 45), theta: 9.6 - 12.8 (0 - 3.2)
    #ipoint = draw_bin(14, 18, 9.6, 12.8, 0, TMath.Pi()/4, 30, ipoint, en_min, en_max, pitheta_max, geom, top)

    #ipoint = draw_bin(14, 18, 9.6, 12.8, TMath.Pi()/4, TMath.Pi()/2, 0, ipoint, en_min, en_max, pitheta_max, geom, top)

    #ipoint = draw_bin(14, 18, 9.6, 12.8, TMath.Pi()/2, TMath.Pi(), 70, ipoint, en_min, en_max, pitheta_max, geom, top)




    #coordinate frame
    coordinates_draw(geom, top, med, ipoint)

    #end of geometry definition
    geom.CloseGeometry()

    #add geometry to TEve
    TEveManager.Create()
    node = gGeoManager.GetTopNode()
    topnode = TEveGeoTopNode(gGeoManager, node)
    topnode.SetVisLevel(4)
    topnode.GetNode().GetVolume().SetVisibility(rt.kFALSE)
    gEve.AddGlobalElement(topnode)

    #labels for coordinate frame
    coordinates_labels(gEve, en_min, en_max, pitheta_max)

    #end of objects definition
    gEve.FullRedraw3D(rt.kTRUE)
    viewer = gEve.GetDefaultGLViewer()

    viewer.SetClearColor(rt.kBlack)
    #viewer.SetClearColor(rt.kWhite)

    code.interact(local=locals())

#main

#_____________________________________________________________________________
def draw_bin(e1, e2, t1, t2, p1, p2, tp, ipoint, en_min, en_max, pitheta_max, geom, top):

    r1 = get_r(e1, en_min, en_max)
    r2 = get_r(e2, en_min, en_max)

    phi = get_phi(p1, p2)

    y1 = get_y(t1, pitheta_max)
    y2 = get_y(t2, pitheta_max)

    geom.Material("Vac_bin"+str(ipoint), 0, 0, 0, 0)
    med_bin = TGeoMedium("Vac_bin"+str(ipoint), ipoint, geom.GetMaterial("Vac_bin"+str(ipoint)))

    tbin = geom.MakeTubs("tbin"+str(ipoint), med_bin, r1, r2, (y1-y2)/2, phi[0], phi[1])
    tbin.SetLineColor(rt.kGreen)
    #tbin.SetLineColor(rt.kRed)
    tbin.SetTransparency(tp)
    top.AddNode(tbin, ipoint, TGeoCombiTrans("tbin_trans", 0, y1-(y1-y2)/2, 0, TGeoRotation("tbin_rot", 0, 90, 90)))
    ipoint += 1

    return ipoint

#draw_bin

#_____________________________________________________________________________
def get_r(en, en_min, en_max):

    #radial coordinate for energy in local scale 0 - 100

    return 100*(en-en_min)/(en_max-en_min)

#get_r

#_____________________________________________________________________________
def get_y(pitheta, pitheta_max):

    #vertical (y) coordinate for pi-theta in mrad, top-down

    return 100-100*pitheta/pitheta_max

#get_y

#_____________________________________________________________________________
def get_phi(phi1, phi2):

    #azimuthal in local coordinates, correct for two rotations to z-x plane

    #rad to deg
    phi1 = phi1*180/TMath.Pi()
    phi2 = phi2*180/TMath.Pi()

    #swap and invert to correct for rotations to z-x
    return -phi2, -phi1

#get_r

#_____________________________________________________________________________
def coordinates_draw(geom, top, med, ipoint):

    #rings representing radial coordinate (energy)
    nrings = 4
    for i in range(nrings):
        ring_r = (i+1)*100/nrings
        ring = geom.MakeTube("ring_"+str(ipoint), med, ring_r, ring_r+1, 0.5)
        ring.SetLineColor(rt.kBlue)
        ring.SetTransparency(0)
        top.AddNode(ring, ipoint, TGeoRotation("idx_rot", 0, 90, 0))
        ipoint += 1

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

#coordinates_draw

#_____________________________________________________________________________
def coordinates_labels(gEve, en_min, en_max, pitheta_max):

    #number of divisions, same as in coordinates_draw()
    nrings = 4
    nrays = 8
    nticks = 5

    font_size = 20

    #tick labels for theta angles
    vbar_txt = TEveText("theta (mrad)")
    vbar_txt.SetFontSize( font_size )
    vbar_txt.SetMainColor( rt.kRed )
    vbar_txt.RefMainTrans().SetPos(0, 110, 100)
    gEve.AddGlobalElement(vbar_txt)
    for i in range(nticks):
        vbar_txt = TEveText( "{0:.1f}".format(i*pitheta_max/nticks) )
        vbar_txt.SetFontSize( font_size )
        vbar_txt.SetMainColor( rt.kRed )
        vbar_txt.RefMainTrans().SetPos(0, 100-i*100/nticks, 110)
        gEve.AddGlobalElement(vbar_txt)

    #tick labels for azimuthal angle rays
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

#coordinates_labels

#_____________________________________________________________________________
if __name__ == "__main__":

    gSystem.Load("libGeom")

    main()















