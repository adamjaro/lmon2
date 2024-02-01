#!/usr/bin/python3

# print pixel indices

import glob, sys

import ROOT as rt
from ROOT import gPad, gROOT, gStyle, TFile, TChain, std, double
from ROOT import TH2I, TH2D

#_____________________________________________________________________________
def main():

    #input and output
    inp = "/home/jaroslav/sim/lmon2-data/taggers/tag9a/????/lmon.root"
    out = "hits.root"

    tree = TChain("DetectorTree")

    for inp in glob.glob(inp):
        tree.Add(inp)

    print("Events:", tree.GetEntries())

    ipix_vec = std.vector(int)()
    irow_vec = std.vector(int)()
    en_vec = std.vector(double)()

    tree.SetBranchAddress("lowQ2_s1_4_ipix", ipix_vec)
    tree.SetBranchAddress("lowQ2_s1_4_irow", irow_vec)
    tree.SetBranchAddress("lowQ2_s1_4_en", en_vec)

    h_counts = TH2D("h_counts", "h_counts", 5455, 0, 5455, 3637, 0, 3637)

    #event loop
    for iev in range(tree.GetEntries()):

        #print("Next event:", iev)

        tree.GetEntry(iev)

        #hit loop
        for ihit in range(ipix_vec.size()):

            #threshold at 0.4 keV
            if en_vec.at(ihit) < 0.4: continue

            ipix = ipix_vec.at(ihit)
            irow = irow_vec.at(ihit)
            #print(ipix, irow, en.at(ihit))

            #increment hit counts at a givin ipix and irow
            h_counts.SetBinContent(ipix, irow, h_counts.GetBinContent(ipix, irow)+1)

        #hit loop

    #event loop

    out = TFile(out, "recreate")
    h_counts.Write("h_counts", 0)
    out.Close()

#main

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()
    gStyle.SetPadTickX(1)
    gStyle.SetFrameLineWidth(2)

    main()

