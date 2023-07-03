#!/usr/bin/python3

from ctypes import c_double

import ROOT as rt
from ROOT import gROOT, gSystem, TFile

#_____________________________________________________________________________
def main():

    gROOT.ProcessLine(".L ../../../base/include/MCParticles.h")

    #input file
    inp = "pwo.root"

    #open the input
    infile = TFile.Open(inp)
    tree = infile.Get("DetectorTree")

    #event attributes
    num_interactions = c_double(0)
    trueQ2 = c_double(0)
    tree.SetBranchAddress("num_interactions", num_interactions)
    tree.SetBranchAddress("true_Q2", trueQ2)

    #MC particles
    mcp = rt.MCParticles.Coll()
    mcp.ConnectInput("mcp", tree)

    #event loop
    for iev in range(tree.GetEntriesFast()):

        tree.GetEntry(iev)

        print("Next event:", iev)

        print("num_interactions:", num_interactions.value)
        print("trueQ2:", trueQ2.value)

        mcp.LoadInput()

        for imc in range(mcp.GetN()):

            mc = mcp.GetUnit(imc)

            #print("mcp:", mc.itrk, mc.pdg, mc.vx, mc.vy, mc.vz, mc.en)

#main

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()

    main()





