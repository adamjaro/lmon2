
//_____________________________________________________________________________
//
// Calorimeter analysis for CalPWO in tagger station
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include <fstream>
#include "glob.h"

//Boost
#include <boost/program_options.hpp>

//ROOT
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

//lmon2 base
#include "GeoParser.h"
#include "LoadXML.h"

//local classes
#include "TagCalPWO.h"
#include "CalPWOClusterWavg.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
void TagCalPWO::Run(const char *conf) {

  cout << "TagCalPWO::Run, " << conf << endl;

  //configuration file
  program_options::options_description opt("opt");
  opt.add_options()
    ("main.input", program_options::value<string>(), "Analysis input")
    ("main.geo", program_options::value<string>(), "Geometry configuration")
    ("main.outfile", program_options::value<string>(), "Output from the analysis")
  ;

  //load the configuration file
  ifstream config(conf);
  program_options::variables_map opt_map;
  program_options::store(program_options::parse_config_file(config, opt), opt_map);

  //geometry
  GeoParser geo;
  LoadXML xml(geo);
  xml.ReadInput( GetStr(opt_map, "main.geo") );

  //inputs
  string input = GetStr(opt_map, "main.input");
  glob_t glob_inputs;
  glob(input.c_str(), GLOB_TILDE, NULL, &glob_inputs);

  //input tree
  TChain tree("DetectorTree");
  for(size_t i=0; i<glob_inputs.gl_pathc; i++) {
    cout << "Adding input: " << glob_inputs.gl_pathv[i] << endl;
    tree.Add( glob_inputs.gl_pathv[i] );
  }

  //calorimeter clusters
  CalPWOClusterWavg cls_s1("lowQ2_s1_pwo");
  cls_s1.ConnectInput(&tree);
  cls_s1.SetGeometry("vac_S1", "lowQ2_s1_pwo", &geo);

  //event loop
  Long64_t nev = tree.GetEntries();
  //Long64_t nev = 124;
  Long64_t iprint = nev/12;
  for(Long64_t iev=0; iev<nev; iev++) {

    tree.GetEntry(iev);

    //cout << endl;

    cls_s1.ProcessEvent();

  }//event loop

  cout << "Events: " << nev << endl;

}//Run

//_____________________________________________________________________________
string TagCalPWO::GetStr(program_options::variables_map& opt_map, std::string par) {

  //string for a given program option

  string res = opt_map[par].as<string>();
  res.erase(remove(res.begin(), res.end(), '\"'), res.end());

  return res;

}//GetStr






















