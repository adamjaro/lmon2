#!/usr/bin/python3

# runs tracking via TagTrackBasic

import sys
import os

import ROOT as rt
from ROOT import gSystem, gInterpreter

#_____________________________________________________________________________
def main():

    #configuration from command line
    config = get_config()

    #analysis library
    gSystem.Load("liblmon2Reco.so")

    #includes to make the task instance
    gInterpreter.Declare('#include <boost/program_options.hpp>')
    gInterpreter.AddIncludePath( os.path.dirname(os.path.abspath(__file__))+"/../include/" )
    gInterpreter.Declare('#include "TagTrackBasic.h"')

    #run the task    
    task = rt.TagTrackBasic()
    task.Run(config)

#main

#_____________________________________________________________________________
def get_config():

    #command line options
    args = sys.argv
    if len(args) < 2:
        print("No configuration specified.")
        quit()
    args.pop(0)

    return args.pop(0)

#get_config

#_____________________________________________________________________________
if __name__ == "__main__":

    main()

