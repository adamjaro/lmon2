#!/usr/bin/python3

import ROOT as rt
from ROOT import gSystem, gInterpreter, gROOT, std
from ROOT import TFile, TF1, TF2, TCanvas

#_____________________________________________________________________________
def main():

    #two-dimensional Gaussian as a training function
    gauss2 = TF2("gauss2","exp(-(x^2)-(y^2))", -3, 3, -3, 3)

    #program options from boost for LPS configuration
    gROOT.ProcessLine("boost::program_options::options_description opt(\"opt\");")

    #LPS instance, solution dimensionality is 1 (single value)
    solver = rt.LookupProblemSolver(1, "Gauss2D", rt.opt)

    #independent quantities, x and y of the Gaussian
    solver.MakeQuantity("x")
    solver.MakeQuantity("y")

    #load 'conf.ini' with optional configuration for the layers
    gROOT.ProcessLine("boost::program_options::variables_map opt_map;")
    gROOT.ProcessLine("boost::program_options::store(parse_config_file(\"conf.ini\", opt), opt_map);")

    #optional non-uniform segmentation in quantities
    solver.SetUseMedCellFinder()

    #initialize the LPS for training
    solver.Initialize(rt.opt_map)

    #linear function to generate values for x and y,
    #non-uniform from 0.3 to 1 over x from -3 to 3
    linear = TF1("linear", "(1.3/2)+(0.7/6)*x", -3, 3)

    #number of training values
    nxy = 10000

    #training loop
    for i in range(nxy):

        #x and y as independent quantities
        x = linear.GetRandom()
        y = linear.GetRandom()

        #Gaussian value as known solution
        gv = gauss2.Eval(x, y)

        #set the quantities and solution to vector containers for LPS
        quant = std.vector("double")([x, y])
        sol = std.vector("double")([gv])

        #add training input to LPS
        solver.AddInput(quant, sol);

    #finalize the LPS for export
    solver.Finalize()

    #export to file
    out = TFile("lps_gaussian_py.root", "recreate")

    solver.Export()

    out.Close()

    #draw the Gaussian which was used for the training
    can = TCanvas("can", "can", 768, 768)

    gauss2.Draw("colz")

    can.SaveAs("gauss_py.pdf")

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

    #run the 'main' function
    main()



















