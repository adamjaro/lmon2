
//_____________________________________________________________________________
//
// Reconstruction for tagger detectors with Timepix4
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include <array>
#include <memory>

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
#include "TrkPlaneBasicHits.h"
#include "IOProgramOptions.h"
#include "TrkHitsTransform.h"
#include "TrkClusterFinder.h"
#include "TagTpix4Reco.h"
#include "TagTrackFinder.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
void TagTpix4Reco::Run(const std::vector<std::string>& argvv) {

  cout << "TagTpix4Reco::Run" << endl;

  //configuration
  program_options::options_description opt("opt");
  opt.add_options()
    ("input", program_options::value<string>(), "Analysis input")
    ("output", program_options::value<string>()->default_value("rec.root"), "Output from the analysis")
    ("config", program_options::value<string>(), "Configuration file")
    ("geo", program_options::value<string>(), "Geometry configuration")
    ("write_hits", program_options::value<bool>()->default_value(false), "Write hits to output")
    ("write_clusters", program_options::value<bool>()->default_value(false), "Write clusters to output")
    ("max_chi2ndf", program_options::value<double>()->default_value(0.01), "Maximal tracks Chi2/NDF")
    ("min_cls_dist", program_options::value<double>()->default_value(0), "Minimal cluster distance")
  ;
  program_options::variables_map opt_map;

  //read the command line options
  const char *argv[argvv.size()]; //char* from vector for the parser
  for(size_t i=0; i<argvv.size(); i++) {
    argv[i] = argvv[i].c_str();
  }
  program_options::store(program_options::parse_command_line(argvv.size(), argv, opt), opt_map);

  //optional configuration file
  if( opt_map.count("config") ) {
    cout << "Config: " << opt_map["config"].as<string>() << endl;
    program_options::store(program_options::parse_config_file(opt_map["config"].as<string>().c_str(), opt), opt_map);
  }

  //input/output from program_options
  IOProgramOptions io(opt_map);

  //input tree
  unique_ptr<TChain> tree = io.MakeChain("input");

  //true event kinematics
  tree->SetBranchAddress("true_Q2", &fTrueQ2);
  tree->SetBranchAddress("true_x", &fTrueX);

  //output file
  unique_ptr<TFile> out = io.MakeFile("output");

  //output tree
  TTree otree("event", "event");
  otree.Branch("true_Q2", &fTrueQ2, "true_Q2/D");
  otree.Branch("true_x", &fTrueX, "true_x/D");

  //geometry
  GeoParser geo;
  LoadXML xml(geo);
  xml.ReadInput( io.GetStr("geo") );

  //hits
  std::array<TrkHitsTransform, 8> hits {
    TrkHitsTransform("lowQ2_s1_1", tree.get()), // tagger 1
    TrkHitsTransform("lowQ2_s1_2", tree.get()),
    TrkHitsTransform("lowQ2_s1_3", tree.get()),
    TrkHitsTransform("lowQ2_s1_4", tree.get()),
    TrkHitsTransform("lowQ2_s2_1", tree.get()), // tagger 2
    TrkHitsTransform("lowQ2_s2_2", tree.get()),
    TrkHitsTransform("lowQ2_s2_3", tree.get()),
    TrkHitsTransform("lowQ2_s2_4", tree.get())
  };
  //output for hits if requested
  if( opt_map["write_hits"].as<bool>() ) {
    for(TrkHitsTransform& i: hits) { i.CreateOutput(&otree); }
  }
  //geometry for hits
  for(size_t i=0; i<4; i++) { hits[i].LoadV6p3_rev3(geo, "vac_S1"); } // tagger 1 section
  for(size_t i=4; i<8; i++) { hits[i].LoadV6p3_rev3(geo, "vac_S2"); } // tagger 2 section

  //cluster finder for all planes
  vector<std::shared_ptr<TrkClusterFinder>> cls;
  cout << "min_cls_dist: " << opt_map["min_cls_dist"].as<double>() << endl;
  for(TrkHitsTransform& i: hits) {
    cls.push_back( make_shared<TrkClusterFinder>(i.GetName()+"_cls", i.GetHits()) ); // cluster names with _cls postfix
    cls.back()->SetLimMdist( opt_map["min_cls_dist"].as<double>() );
    if( opt_map["write_clusters"].as<bool>() ) { cls.back()->CreateOutput(&otree); }
  }

  //track finder
  cout << "max_chi2ndf: " << opt_map["max_chi2ndf"].as<double>() << endl;
  TagTrackFinder s1("s1_tracks", cls);
  s1.LoadV6p3_rev3(geo, "lowQ2_s1_1", "lowQ2_s1_2", "zpos");
  s1.SetMaxChi2Ndf( opt_map["max_chi2ndf"].as<double>() );
  s1.CreateOutput(&otree);
  TagTrackFinder s2("s2_tracks", cls, 4); // starting offset for planes in tagger 2
  s2.LoadV6p3_rev3(geo, "lowQ2_s2_1", "lowQ2_s2_2", "zpos");
  s2.SetMaxChi2Ndf( opt_map["max_chi2ndf"].as<double>() );
  s2.CreateOutput(&otree);

  Long64_t nev = tree->GetEntries();
  //nev = 12;

  cout << "Number of events: " << nev << endl;

  Long64_t iprint = nev/12;

  //event loop
  for(Long64_t iev=0; iev<nev; iev++) {

    if( iev > 0 and iev%iprint == 0 ) {
      cout << Form("%.1f", 100.*iev/nev) << "%" << endl;
    }

    tree->GetEntry(iev);

    //call to ProcessEvent for all hits
    for_each(hits.begin(), hits.end(), mem_fn(&TrkHitsTransform::ProcessEvent));

    //cls.ProcessEvent();
    for_each(cls.begin(), cls.end(), mem_fn(&TrkClusterFinder::ProcessEvent));

    s1.ProcessEvent();
    s2.ProcessEvent();

    //FinishEvent for clusters after track finder for cluster information about tracks
    for_each(cls.begin(), cls.end(), mem_fn(&TrkClusterFinder::FinishEvent));

    otree.Fill();

  }//event loop

  otree.Write(0, TObject::kOverwrite);

  out->Close();

  //print hit counts
  for_each(hits.begin(), hits.end(), mem_fn(&TrkHitsTransform::PrintCounts));


}//Run

























