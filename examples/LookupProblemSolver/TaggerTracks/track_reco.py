#!/usr/bin/python3

from ctypes import c_double

import ROOT as rt
from ROOT import gSystem, gInterpreter, gROOT, std, gRandom
from ROOT import TFile, TTree, TObject

#_____________________________________________________________________________
def main():

    #tracks physics (photoproduction) input for kinematics reconstruction
    inp = TFile("trk_physics_20kevt.root", "read")

    #input tree
    event_tree = inp.Get("event")

    #attach true kinematics to compare with solutions by LPS
    true_en = c_double(0)
    true_theta = c_double(0)
    true_phi = c_double(0)
    event_tree.SetBranchAddress("true_el_E", true_en)
    event_tree.SetBranchAddress("true_el_theta", true_theta)
    event_tree.SetBranchAddress("true_el_phi", true_phi)

    #input tracks
    tracks = rt.TaggerTracks.Coll()
    tracks.ConnectInput("s1_tracks", event_tree)

    #LPS instance for obtaining solutions
    solver = rt.LookupProblemSolver(3, "Tagger_1")

    #import trained LPS
    inp_lps = TFile("lps_tracks.root", "read")
    solver.Import(inp_lps)

    #output file and tree on track kinematics by LPS solutions

    #create the file and tree
    out = TFile("rec_tracks.root", "recreate")
    otree = TTree("rec_trk", "rec_trk")

    #branches with reconstructed kinematics by LPS
    rec_en = c_double(0)
    rec_theta = c_double(0)
    rec_phi = c_double(0)
    otree.Branch("rec_en", rec_en, "rec_en/D")
    otree.Branch("rec_theta", rec_theta, "rec_theta/D")
    otree.Branch("rec_phi", rec_phi, "rec_phi/D")

    #branches with true kinematics for comparison
    otree.Branch("true_en", true_en, "true_en/D")
    otree.Branch("true_theta", true_theta, "true_theta/D")
    otree.Branch("true_phi", true_phi, "true_phi/D")

    #event loop
    for iev in range(event_tree.GetEntries()):

        #get the event and load the tracks
        event_tree.GetEntry(iev)
        tracks.LoadInput()

        #tracks loop
        for trk in tracks.GetReadData():

            #independent quantities for LPS from track geometry
            quant = std.vector("double")([trk.x, trk.y, trk.theta_x, trk.theta_y])

            #prepare the vector container for the solutions
            sol = std.vector("double")()
            sol.resize(3) # length is 3

            #obtain the solution and check its status (passed by ref to sol)
            stat = solver.Solve(quant, sol)
            if stat == False: continue

            #put the obtained solution to output tree
            rec_en.value = sol.at(0)
            rec_theta.value = sol.at(1)
            rec_phi.value = sol.at(2)

            #fill the output tree (also true kinematics will be there)
            otree.Fill()

            #print(sol.at(0), sol.at(1), sol.at(2))

    #write the output tree and close the output file
    otree.Write("", TObject.kOverwrite)

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
    gInterpreter.AddIncludePath("../../../reconstruction/include/");
    gInterpreter.Declare('#include "LookupProblemSolver.h"')

    #header for tracks representation
    gInterpreter.Declare('#include "TaggerTracks.h"')

    #run the 'main' function
    main()

