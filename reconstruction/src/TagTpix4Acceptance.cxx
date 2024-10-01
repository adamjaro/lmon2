
//_____________________________________________________________________________
//
// Acceptance analysis for tagger detectors with Timepix4
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include <array>

//Boost
#include <boost/program_options.hpp>

//ROOT
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

//local classes
#include "TrkPlaneBasicHits.h"
#include "IOProgramOptions.h"
#include "TagTpix4Acceptance.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
void TagTpix4Acceptance::Run(const std::vector<std::string>& argvv) {

  cout << "Tpix4HitAnalysis::Run" << endl;

  //configuration
  program_options::options_description opt("opt");
  opt.add_options()
    ("input", program_options::value<string>(), "Analysis input")
    ("output", program_options::value<string>()->default_value("acc.root"), "Output from the analysis")
    ("config", program_options::value<string>(), "Configuration file")
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

  Bool_t s1_allhit, s2_allhit;
  otree.Branch("s1_allhit", &s1_allhit, "s1_allhit/O");
  otree.Branch("s2_allhit", &s2_allhit, "s2_allhit/O");

  //plane hits, tagger 1
  std::array<plane_hits, 4> planes_s1{
    plane_hits("lowQ2_s1_1", tree.get()),
    plane_hits("lowQ2_s1_2", tree.get()),
    plane_hits("lowQ2_s1_3", tree.get()),
    plane_hits("lowQ2_s1_4", tree.get())
  };

  //plane hits, tagger 2
  std::array<plane_hits, 4> planes_s2{
    plane_hits("lowQ2_s2_1", tree.get()),
    plane_hits("lowQ2_s2_2", tree.get()),
    plane_hits("lowQ2_s2_3", tree.get()),
    plane_hits("lowQ2_s2_4", tree.get())
  };

  Long64_t nev = tree->GetEntries();
  //nev = 12;

  cout << "Number of events: " << nev << endl;

  Long64_t nsel1=0, nsel2=0;

  //event loop
  for(Long64_t iev=0; iev<nev; iev++) {

    //cout << "iev: " << iev << endl;

    tree->GetEntry(iev);

    //hits in tagger 1
    unsigned long ns1, ns2, ns3, ns4;
    ns1 = planes_s1[0].GetNhits();
    ns2 = planes_s1[1].GetNhits();
    ns3 = planes_s1[2].GetNhits();
    ns4 = planes_s1[3].GetNhits();

    //hit in all tagger 1 planes
    s1_allhit = kFALSE;
    if( ns1!=0 and ns2!=0 and ns3!=0 and ns4!=0 ) {
      s1_allhit = kTRUE;
      nsel1++;
    }

    //hits in tagger2
    ns1 = planes_s2[0].GetNhits();
    ns2 = planes_s2[1].GetNhits();
    ns3 = planes_s2[2].GetNhits();
    ns4 = planes_s2[3].GetNhits();

    //hit in all tagger 2 planes
    s2_allhit = kFALSE;
    if( ns1!=0 and ns2!=0 and ns3!=0 and ns4!=0 ) {
      s2_allhit = kTRUE;
      nsel2++;
    }

    otree.Fill();

  }//event loop

  otree.Write(0, TObject::kOverwrite);

  out->Close();

  cout << "Selected: " << nsel1 << " " << nsel2 << endl;

}//Run

//_____________________________________________________________________________
TagTpix4Acceptance::plane_hits::plane_hits(std::string nam, TTree *in_tree) {

  //input hits
  hits = make_unique<TrkPlaneBasicHits::Coll>();
  hits->ConnectInput(nam, in_tree);

}//plane_hits

//_____________________________________________________________________________
unsigned long TagTpix4Acceptance::plane_hits::GetNhits() {

  hits->LoadInput();

  return hits->GetN();

}//GetNhits























