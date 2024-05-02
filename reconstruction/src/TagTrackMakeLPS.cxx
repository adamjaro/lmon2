
//_____________________________________________________________________________
//
// Task to create reconstruction response by LookupProblemSolver (LPS)
// using tracks by TagTrackFindBasic
//
// Runs by  run_TagTrackMakeLPS.py
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
#include "TFile.h"

//lmon2 base
#include "GeoParser.h"
#include "LoadXML.h"

//local classes
#include "TagTrackMakeLPS.h"
#include "TagTrackFindBasic.h"
#include "LookupProblemSolver.h"

using namespace std;
using namespace boost;
using namespace boost::program_options;

//_____________________________________________________________________________
void TagTrackMakeLPS::Run(const char *conf) {

  cout << "TagTrackMakeLPS::Run, " << conf << endl;

  //configuration file
  options_description opt("opt");
  opt.add_options()
    ("main.input", program_options::value<string>(), "Analysis input")
    ("main.geo", program_options::value<string>(), "Geometry configuration")
    ("main.outfile", program_options::value<string>(), "Output from the analysis")
    ("main.max_chi2ndf", program_options::value<double>(), "Maximal tracks Chi2/NDF")
    ("main.min_cls_dist", program_options::value<double>(), "Minimal cluster distance")
    ("main.tracks_output", program_options::value<bool>()->default_value(true), "Write output for tracks")
    ("main.use_FoamCellFinder", program_options::value<bool>()->default_value(false), "Binning by TFoam")
    ("main.use_MedCellFinder", program_options::value<bool>()->default_value(false), "Binning by median split")
    ("main.use_SumCellFinder", program_options::value<bool>()->default_value(false), "Binning by even sum per bin")
  ;

  //reconstruction for tagger stations
  LookupProblemSolver *s1_rec = new LookupProblemSolver(3, "s1", opt);

  //quantities measured by the taggers
  s1_rec->MakeQuantity("x"); // x position, mm
  s1_rec->MakeQuantity("y"); // y position, mm
  s1_rec->MakeQuantity("theta_x"); // theta_x angle, rad
  s1_rec->MakeQuantity("theta_y"); // theta_y angle, rad

  //load the configuration file
  variables_map opt_map;
  store(parse_config_file(conf, opt), opt_map);

  //initialize the reconstruction
  s1_rec->Initialize(opt_map);

  s1_rec->SetUseFoamCellFinder( opt_map["main.use_FoamCellFinder"].as<bool>() );
  s1_rec->SetUseMedCellFinder( opt_map["main.use_MedCellFinder"].as<bool>() );
  s1_rec->SetUseSumCellFinder( opt_map["main.use_SumCellFinder"].as<bool>() );

  //inputs
  string input = GetStr(opt_map, "main.input");
  cout << "Input: " << input << endl;
  glob_t glob_inputs;
  glob(input.c_str(), GLOB_TILDE, NULL, &glob_inputs);

  //input tree
  TChain tree("DetectorTree");
  for(size_t i=0; i<glob_inputs.gl_pathc; i++) {
    cout << "Adding input: " << glob_inputs.gl_pathv[i] << endl;
    tree.Add( glob_inputs.gl_pathv[i] );
  }

  //input true event kinematics
  tree.SetBranchAddress("true_el_E", &fTrueEn);
  tree.SetBranchAddress("true_el_theta", &fTrueTheta);
  tree.SetBranchAddress("true_el_phi", &fTruePhi);

  //geometry
  GeoParser geo;
  LoadXML xml(geo);
  xml.ReadInput( GetStr(opt_map, "main.geo") );
  //geo.PrintAll();

  //output file
  string outfile = GetStr(opt_map, "main.outfile");
  cout << "Output: " << outfile << endl;
  TFile out(outfile.c_str(), "recreate");

  //interaction (event) output tree
  TTree otree("event", "event");

  //track finder
  TagTrackFindBasic s1("s1");
  TagTrackFindBasic s2("s2");

  s1.SetGeometry(&geo);
  s2.SetGeometry(&geo);

  s1.ConnectHitsInput(&tree);
  s2.ConnectHitsInput(&tree);

  bool tracks_output = opt_map["main.tracks_output"].as<bool>();
  cout << "Tracks output: " << tracks_output << endl;
  s1.CreateOutput(&otree, false, tracks_output); // no output for planes, tracks as required
  s2.CreateOutput(&otree, false, tracks_output);

  //selection criteria for tracking
  if( opt_map.find("main.max_chi2ndf") != opt_map.end() ) {
    Double_t x = opt_map["main.max_chi2ndf"].as<double>();
    s1.SetMaxChi2Ndf(x);
    s2.SetMaxChi2Ndf(x);
  }
  if( opt_map.find("main.min_cls_dist") != opt_map.end() ) {
    Double_t x = opt_map["main.min_cls_dist"].as<double>();
    s1.SetClsLimMdist(x);
    s2.SetClsLimMdist(x);
  }

  //counters
  Long64_t ncls_s1=0, ncls_s2=0, ntrk_s1=0, ntrk_s2=0;

  //event loop
  Long64_t nev = tree.GetEntries();
  //Long64_t nev = 1200;
  Long64_t iprint = nev/12;
  for(Long64_t iev=0; iev<nev; iev++) {

    tree.GetEntry(iev);

    if( iev > 0 and iev%iprint == 0 ) {
      cout << Form("%.1f", 100.*iev/nev) << "%" << endl;
    }

    //cout << "Next event" << endl;

    s1.ProcessEvent();
    s2.ProcessEvent();

    AddInput(&s1, s1_rec);
    //AddInput(&s2, s2_rec);

    s1.FinishEvent();
    s2.FinishEvent();

    //counters
    ncls_s1 += s1.GetNumberOfClusters();
    ncls_s2 += s2.GetNumberOfClusters();
    ntrk_s1 += s1.GetTracks().size();
    ntrk_s2 += s2.GetTracks().size();

    //fill event tree
    otree.Fill();

  }//event loop

  cout << "Finalizing..." << endl;
  s1_rec->Finalize();
  //s2_rec->Finalize();

  otree.Write(0, TObject::kOverwrite);

  cout << "Exporting..." << endl;
  s1_rec->Export();
  //s2_rec->Export();

  out.Close();

  cout << "Events: " << nev << endl;
  cout << "Clusters, s1: " << ncls_s1 << ", s2: " << ncls_s2 << endl;
  cout << "Tracks, s1: " << ntrk_s1 << ", s2: " << ntrk_s2 << endl;

}//Run

//_____________________________________________________________________________
void TagTrackMakeLPS::AddInput(TagTrackFindBasic *tag, LookupProblemSolver *rec) {

  //event with just one track
  if( tag->GetTracks().size() != 1 ) return;

  //track loop
  for(TagTrackFindBasic::Track& i: tag->GetTracks()) {

    //add input for reconstruction
    vector<double> quant = {i.x, i.y, i.theta_x, i.theta_y};
    vector<double> sol = {fTrueEn, fTrueTheta, fTruePhi};
    rec->AddInput(quant, sol);

  }//track loop

}//AddInput

//_____________________________________________________________________________
string TagTrackMakeLPS::GetStr(program_options::variables_map& opt_map, std::string par) {

  string res = opt_map[par].as<string>();
  res.erase(remove(res.begin(), res.end(), '\"'), res.end());

  return res;

}//GetStr




















