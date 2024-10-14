#!/usr/bin/python3

# runs TagTpix4Reco

import sys
import os

import ROOT as rt
from ROOT import gROOT, gSystem, gInterpreter, std

#_____________________________________________________________________________
def main():

    #analysis library
    gSystem.Load("liblmon2Reco.so")

    #includes to make the task instance
    gInterpreter.Declare('#include <boost/program_options.hpp>')
    gInterpreter.AddIncludePath( os.path.dirname(os.path.abspath(__file__))+"/../include/" )
    gInterpreter.Declare('#include "TagTpix4Reco.h"')

    #run the task    
    task = rt.TagTpix4Reco()
    task.Run(get_config())

#main

#_____________________________________________________________________________
def get_config():

    #command line options
    if len(sys.argv) < 2:
        print("No configuration specified.")
        quit()

    #store command line options in vector to pass to the task
    x=std.vector(str)()
    for i in sys.argv:
        x.push_back(i)

    return x

#get_config

#_____________________________________________________________________________
if __name__ == "__main__":

    gROOT.SetBatch()

    main()






