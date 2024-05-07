#!/usr/bin/python3

import ROOT as rt
from ROOT import gSystem, gInterpreter, gROOT, std, gRandom
from ROOT import TFile, TGraph2D, TCanvas

#_____________________________________________________________________________
def main():

    #LPS instance without program_options for obtaining solutions
    solver = rt.LookupProblemSolver(1, "Gauss2D") #dimensionality and name match training

    #import trained LPS
    inp = TFile("lps_gaussian_py.root", "read")
    #inp = TFile("lps_gaussian.root", "read")
    solver.Import(inp)

    #number of points to obtain solutions
    nxy = 1000

    #graph to draw the solutions
    g2 =TGraph2D(nxy)

    #point loop
    for i in range(nxy):

        #random x and y uniformly from -3 to 3
        x = -3 + 6*gRandom.Rndm()
        y = -3 + 6*gRandom.Rndm()

        #put x and y (independent quantities) to vector container
        quant = std.vector("double")([x, y])

        #prepare vector container for solution
        sol = std.vector("double")()
        sol.resize(1) # set length as 1

        #obtain the solution for x and y (passed to 'sol')
        stat = solver.Solve(quant, sol)

        #check status if solution was found
        if stat == False:
            continue

        #mark the obtained solution in output graph
        g2.SetPoint(i, x, y, sol[0])

    #draw the graph
    can = TCanvas("can", "can", 768, 768)

    g2.Draw("surf2")

    can.SaveAs("gauss2_py.pdf")

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

    #run the 'main' function
    main()
















