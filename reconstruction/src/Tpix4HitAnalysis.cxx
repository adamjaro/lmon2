
//_____________________________________________________________________________
//
// Hit analysis for Timepix4
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include <fstream>
#include "glob.h"
#include <array>
#include <time.h>
#include <atomic>
#include <future>
#include <thread>

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
#include "Tpix4HitAnalysis.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
void Tpix4HitAnalysis::Run(const char *conf) {

  cout << "Tpix4HitAnalysis::Run, " << conf << endl;

  //configuration file
  program_options::options_description opt("opt");
  opt.add_options()
    ("main.input", program_options::value<string>(), "Analysis input")
    ("main.geo", program_options::value<string>(), "Geometry configuration")
    ("main.outfile", program_options::value<string>(), "Output from the analysis")
    ("main.nev", program_options::value<long>()->default_value(0), "Maximal number of events")
    ("main.nofs", program_options::value<long>()->default_value(0), "Offset in processed events")
    ("main.make_tree", program_options::value<bool>()->default_value(true), "Make also output tree")
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

  //output file
  string outfile = GetStr(opt_map, "main.outfile");
  cout << "Output: " << outfile << endl;
  TFile out(outfile.c_str(), "recreate");

  //flag for output tree
  bool make_tree = opt_map["main.make_tree"].as<bool>();

  //plane hits
  std::array<plane_hits, 8> ahits{
    plane_hits("lowQ2_s1_1", &geo, &tree, make_tree),
    plane_hits("lowQ2_s1_2", &geo, &tree, make_tree),
    plane_hits("lowQ2_s1_3", &geo, &tree, make_tree),
    plane_hits("lowQ2_s1_4", &geo, &tree, make_tree),
    plane_hits("lowQ2_s2_1", &geo, &tree, make_tree),
    plane_hits("lowQ2_s2_2", &geo, &tree, make_tree),
    plane_hits("lowQ2_s2_3", &geo, &tree, make_tree),
    plane_hits("lowQ2_s2_4", &geo, &tree, make_tree)
  };

  //number of events to analyze
  Long64_t nev = opt_map["main.nev"].as<long>();
  if( nev == 0 ) {
    //number of events from tree when not requested specific (left default at zero)
    nev = tree.GetEntries();
  }

  //offset in events
  Long64_t nofs = opt_map["main.nofs"].as<long>();

  cout << "Number of events: " << nev << endl;
  if( nofs != 0 ) { cout << "Offset: " << nofs << endl; }

  //print progress bar at regular intervals
  atomic<Long64_t> prog_iev(0);
  atomic<bool> prog_keep(true);
  auto prog_bar = [this, nev, &prog_iev, &prog_keep]() {
    while(prog_keep.load() == true) {
      this_thread::sleep_for(1s);
      this->ShowProgress(prog_iev.load(), nev);
    }
    this->ShowProgress(nev, nev);
    cout << endl;
  };

  //async task for the progress bar
  future prog_fut = async(prog_bar);

  //event loop
  for(Long64_t iev=nofs; iev<nev+nofs; iev++) {

    //current event for progress bar
    prog_iev.store(iev-nofs);

    tree.GetEntry(iev);

    //cout << "iev: " << iev << endl;

    for_each(ahits.begin(), ahits.end(), std::mem_fn(&plane_hits::run_event));

  }//event loop

  //stop the task for progress bar
  prog_keep.store(false);

  //close the output
  for(const auto& i: ahits) {
    i.h_counts.Write(i.h_counts.GetName(), 0);
    if( i.h_tree != nullptr ) i.h_tree->Write(0, TObject::kOverwrite);
  }
  out.Close();

}//Run

//_____________________________________________________________________________
Tpix4HitAnalysis::plane_hits::plane_hits(std::string nam, GeoParser *geo, TTree *in_tree, bool make_tree) {

  //input hits
  hits = make_unique<TrkPlaneBasicHits::Coll>();
  hits->ConnectInput(nam, in_tree);

  //counts in individual pixels
  Int_t nx = geo->GetD(nam, "nx");
  Int_t ny = geo->GetD(nam, "ny");
  //cout << "nx, ny: " << nx << " " << ny << endl;
  h_counts.SetBins(nx, 0, nx, ny, 0, ny);
  h_counts.SetNameTitle((nam+"_h_counts").c_str(), (nam+"_h_counts").c_str());

  //tree on individual hits
  if( make_tree ) {
    h_tree = make_unique<TTree>();
    h_tree->SetNameTitle((nam+"_h_tree").c_str(), (nam+"_h_tree").c_str());

    auto make_branch = [&](string nam, string type, auto *x) {
      h_tree->Branch(nam.c_str(), x, (nam+"/"+type).c_str());
    };
    make_branch("ipix", "I", &ipix);
    make_branch("irow", "I", &irow);
    make_branch("en", "D", &en);
  }

}//plane_hits::plane_hits

//_____________________________________________________________________________
void Tpix4HitAnalysis::plane_hits::run_event() {

  hits->LoadInput();

  //hits loop
  for(const TrkPlaneBasicHits::Hit& i: hits->GetReadData()) {
    if( i.en < 1e-12 ) continue;

    //pixel location and energy
    ipix = i.ipix;
    irow = i.irow;
    en = i.en;

    //fill the output tree when present
    if( h_tree != nullptr ) {
      h_tree->Fill();
    }

    //cout << ipix << " " << irow << " " << en << endl;

    //threshold at 0.4 keV
    if( en > 0.4 ) {
      //increment hit counts at a givin ipix and irow
      h_counts.SetBinContent(ipix, irow, h_counts.GetBinContent(ipix, irow)+1);
    }

  }//hits loop

}//plane_hits::run_event

//_____________________________________________________________________________
string Tpix4HitAnalysis::GetStr(program_options::variables_map& opt_map, std::string par) {

  //string for a given program option

  string res = opt_map[par].as<string>();
  res.erase(remove(res.begin(), res.end(), '\"'), res.end());

  return res;

}//GetStr

//_____________________________________________________________________________
void Tpix4HitAnalysis::ShowProgress(Double_t xi, Double_t xall) {

  //print proggress bar
  Int_t pbar = 43;

  Double_t prog = xi/xall;
  Int_t pos = prog*pbar;

  string pdone = "\u2588"; //full block
  string phead = "\u2591"; //light shade
  string prem = "\u2591"; //light shade
  cout << " [";
  //proggress loop
  for(Int_t i=0; i<pbar; ++i) {
    if( i < pos ) {
      cout << pdone;
    } else {
      if( i == pos ) {
        cout << phead;
      } else {
        cout << prem;
      }
    }
  }//proggress loop
  //cout << "] " << Form("%.1f", prog*100.)<<" %\r";
  cout << Form("] %.1f", prog*100.)<<" %\r";
  cout.flush();

}//ShowProgress












