#!/usr/bin/python3

from ctypes import c_double

import ROOT as rt
from ROOT import gSystem, gInterpreter, gROOT, std
from ROOT import TFile

#_____________________________________________________________________________
def main():

    #tracks input for LPS training
    inp = TFile("trk_training_400kevt.root", "read")

    #input tree with tracks and true kinematics
    event_tree = inp.Get("event")

    #attach the true kinematics (electron energy and theta and phi angles)
    true_en = c_double(0)
    true_theta = c_double(0)
    true_phi = c_double(0)
    event_tree.SetBranchAddress("true_el_E", true_en)
    event_tree.SetBranchAddress("true_el_theta", true_theta)
    event_tree.SetBranchAddress("true_el_phi", true_phi)

    #input tracks
    tracks = rt.TaggerTracks.Coll()
    tracks.ConnectInput("s1_tracks", event_tree)

    #program options from boost for LPS configuration
    gROOT.ProcessLine("boost::program_options::options_description opt(\"opt\");")

    #LPS instance, solution dimensionality is 3 (energy, theta, phi)
    solver = rt.LookupProblemSolver(3, "Tagger_1", rt.opt)

    #independent quantities as track geometry parameters
    solver.MakeQuantity("x") # x position, mm
    solver.MakeQuantity("y") # y position, mm
    solver.MakeQuantity("theta_x") # theta_x angle, rad
    solver.MakeQuantity("theta_y") # theta_y angle, rad

    #load 'conf.ini' with optional configuration for the layers
    gROOT.ProcessLine("boost::program_options::variables_map opt_map;")
    gROOT.ProcessLine("boost::program_options::store(parse_config_file(\"conf.ini\", opt), opt_map);")

    #optional non-uniform segmentation in quantities
    #solver.SetUseMedCellFinder();

    #initialize the LPS for training
    solver.Initialize(rt.opt_map)

    #input event loop
    for iev in range(event_tree.GetEntries()):

        #load the current event
        event_tree.GetEntry(iev)
        tracks.LoadInput()

        #select events with one track
        if tracks.GetN() != 1: continue

        #get the track (it is known there is one)
        trk = tracks.GetUnit(0)

        #set independent quantities for LPS from track geometry
        quant = std.vector("double")([trk.x, trk.y, trk.theta_x, trk.theta_y])

        #set the known solution from true kinematics, length is 3
        sol = std.vector("double")([true_en.value, true_theta.value, true_phi.value])

        #add training input to LPS
        solver.AddInput(quant, sol)

    #finalize the LPS for export
    solver.Finalize()

    #export to file
    out = TFile("lps_tracks.root", "recreate")

    solver.Export()

    out.Close()

    print("All done")

#_____________________________________________________________________________
if __name__ == "__main__":

    #load shared library for the LPS
    gSystem.Load("../LPS_library/build/libLPS.so")

    #set ROOT environment to work with LPS

    #boost
    gInterpreter.Declare('#include <boost/program_options.hpp>')

    #header for LPS
    gInterpreter.AddIncludePath("../../../reconstruction/include/")
    gInterpreter.Declare('#include "LookupProblemSolver.h"')

    #header for tracks representation
    gInterpreter.Declare('#include "TaggerTracks.h"')

    #run the 'main' function
    main()
