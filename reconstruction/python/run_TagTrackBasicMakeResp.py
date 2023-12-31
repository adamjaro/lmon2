#!/usr/bin/python3

# reconstruction response creation by TagTrackBasicMakeResp

import sys
#import os

from ctypes import CDLL, c_char_p, c_void_p

#import ROOT as rt
#from ROOT import gSystem, gInterpreter

#_____________________________________________________________________________
def main():

    #configuration from command line
    config = get_config()

    #analysis library
    #gSystem.Load("liblmon2Reco.so")
    lib = CDLL("liblmon2Reco.so")

    #includes
    #gInterpreter.Declare('#include <boost/program_options.hpp>')
    #gInterpreter.AddIncludePath( os.path.dirname(os.path.abspath(__file__))+"/../include/" )
    #gInterpreter.Declare('#include "TagTrackBasicMakeResp.h"')

    #analysis task
    #task = rt.TagTrackBasicMakeResp()
    #task.Run(config)
    lib.make_TagTrackBasicMakeResp.restype = c_void_p
    task = c_void_p(lib.make_TagTrackBasicMakeResp())

    lib.run_TagTrackBasicMakeResp( task, c_char_p(bytes(config, "utf-8")) )

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

