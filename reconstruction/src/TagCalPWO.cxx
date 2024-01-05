
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
    ("main.cal_name", program_options::value<string>()->default_value("lowQ2_s1_pwo"), "Calorimeter module name")
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

  //output
  string outfile = GetStr(opt_map, "main.outfile");
  cout << "Output: " << outfile << endl;
  TFile out(outfile.c_str(), "recreate");

  //interaction (event) output tree
  TTree otree("event", "event");

  //calorimeter clusters
  string cal_name = GetStr(opt_map, "main.cal_name"); // name for calorimeter volume
  CalPWOClusterWavg cls_s1(cal_name); // name as in the geometry
  cls_s1.ConnectInput(&tree);
  cls_s1.SetGeometry("vac_S1", cal_name, &geo);
  cls_s1.CreateOutput(&otree);

  //reconstruction counters
  Long64_t ncls_s1=0;

  //event loop
  Long64_t nev = tree.GetEntries();
  //Long64_t nev = 124;
  Long64_t iprint = nev/12;
  for(Long64_t iev=0; iev<nev; iev++) {

    tree.GetEntry(iev);

    //cout << endl;

    cls_s1.ProcessEvent();

    //update the counters
    ncls_s1 += cls_s1.GetClusters().size();

    //fill event tree
    otree.Fill();

  }//event loop

  //close the output
  otree.Write(0, TObject::kOverwrite);
  out.Close();

  cout << "Events: " << nev << endl;
  cout << "Calorimeter clusters: " << ncls_s1 << endl;

}//Run

//_____________________________________________________________________________
string TagCalPWO::GetStr(program_options::variables_map& opt_map, std::string par) {

  //string for a given program option

  string res = opt_map[par].as<string>();
  res.erase(remove(res.begin(), res.end(), '\"'), res.end());

  return res;

}//GetStr






















