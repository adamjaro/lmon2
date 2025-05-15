
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
#include "MCParticles.h"

//local classes
#include "TrkPlaneBasicHits.h"
#include "IOProgramOptions.h"
#include "TrkHitsTransform.h"
#include "TrkClusterFinder.h"
#include "TagTrackFinder.h"
#include "LookupProblemSolver.h"
#include "TagTpix4Reco.h"

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
    ("beam_energy", program_options::value<double>()->default_value(0), "Beam energy (GeV), used in reconstruction")
    ("write_hits", program_options::value<bool>()->default_value(false), "Write hits to output")
    ("write_clusters", program_options::value<bool>()->default_value(false), "Write clusters to output")
    ("nev", program_options::value<long>()->default_value(0), "Maximal number of events")
    ("max_chi2ndf", program_options::value<double>()->default_value(0.01), "Maximal tracks Chi2/NDF")
    ("min_cls_dist", program_options::value<double>()->default_value(0), "Minimal cluster distance")
    ("write_event_tree", program_options::value<bool>()->default_value(true), "Write main event tree")
    ("lps_mode", program_options::value<int>()->default_value(0), "Mode to use LPS, 1:make, 2:rec")
    ("use_FoamCellFinder", program_options::value<bool>()->default_value(false), "Binning by TFoam")
    ("use_MedCellFinder", program_options::value<bool>()->default_value(false), "Binning by median split")
    ("use_SumCellFinder", program_options::value<bool>()->default_value(false), "Binning by even sum per bin")
    ("input_resp", program_options::value<string>(), "Input response for reconstruction")
    ("min_resp_ninp", program_options::value<int>()->default_value(-1), "Minimal number of inputs for reconstruction")
    ("max_resp_rel_en", program_options::value<double>()->default_value(-1), "Maximal relative error in energy")
    ("max_resp_rel_phi", program_options::value<double>()->default_value(-1), "Maximal relative error in phi")
    ("max_resp_err_theta", program_options::value<double>()->default_value(-1), "Maximal error in theta (rad)")
  ;
  program_options::variables_map opt_map;

  //LPS
  unique_ptr<LookupProblemSolver> s1_rec = make_unique<LookupProblemSolver>(3, "s1", &opt);
  unique_ptr<LookupProblemSolver> s2_rec = make_unique<LookupProblemSolver>(3, "s2", &opt);

  //read the command line options
  const char *argv[argvv.size()]; //char* from vector for the parser
  for(size_t i=0; i<argvv.size(); i++) {
    argv[i] = argvv[i].c_str();
  }
  program_options::store(program_options::parse_command_line(argvv.size(), argv, opt), opt_map);

  //optional configuration file
  if( opt_map.count("config") ) {
    cout << "Config: " << opt_map["config"].as<string>() << endl;
    program_options::parsed_options popt = program_options::parse_config_file(opt_map["config"].as<string>().c_str(), opt);
    cout << "Parsed options from config: " << popt.options.size() << endl;
    program_options::store(popt, opt_map);
  }

  //LPS configuration
  int lps_mode = opt_map["lps_mode"].as<int>();
  cout << "lps_mode: " << lps_mode << endl;
  if( lps_mode == 1 ) {
    SetMakeLPS(s1_rec, opt_map);
    SetMakeLPS(s2_rec, opt_map);
  }
  if( lps_mode == 2 ) {
    ImportLPS(s1_rec, opt_map);
    ImportLPS(s2_rec, opt_map);
  }

  //input/output from program_options
  IOProgramOptions io(opt_map);

  //input tree
  unique_ptr<TChain> tree = io.MakeChain("input");

  //true event kinematics
  tree->SetBranchAddress("true_el_E", &fTrueEn);
  tree->SetBranchAddress("true_el_theta", &fTrueTheta);
  tree->SetBranchAddress("true_el_phi", &fTruePhi);
  tree->SetBranchAddress("true_Q2", &fTrueQ2);
  tree->SetBranchAddress("true_x", &fTrueX);
  tree->SetBranchAddress("num_interactions", &fNumInteractions);

  //MC particles
  fMC = make_shared<MCParticles::Coll>();
  fMC->ConnectInput("mcp", tree.get());

  //output file
  unique_ptr<TFile> out = io.MakeFile("output");

  //output event tree
  TTree otree("event", "event");
  otree.Branch("true_el_E", &fTrueEn, "true_el_E/D");
  otree.Branch("true_el_theta", &fTrueTheta, "true_el_theta/D");
  otree.Branch("true_el_phi", &fTruePhi, "true_el_phi/D");
  otree.Branch("true_Q2", &fTrueQ2, "true_Q2/D");
  otree.Branch("true_x", &fTrueX, "true_x/D");
  otree.Branch("num_interactions", &fNumInteractions, "num_interactions/D");

  //flag to write output event tree (assumed to use with creating LPS response)
  bool write_event_tree = opt_map["write_event_tree"].as<bool>();

  //geometry
  GeoParser geo;
  LoadXML xml(geo);
  xml.ReadInput( io.GetStr("geo") );

  //beam energy
  fBeamEn = opt_map["beam_energy"].as<double>();

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

  //number of events to analyze, from tree by default
  Long64_t nev = opt_map["nev"].as<long>();
  if( nev == 0 ) nev = tree->GetEntries();

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

    if( lps_mode == 1 ) {

      //inputs to create LPS response
      AddInput(s1.GetTracks(), s1_rec);
      AddInput(s2.GetTracks(), s2_rec);
    } else {

      //load the MC particles
      fMC->LoadInput();

      //Q2 reconstruction and MC-track assignment
      ElectronRec(s1.GetTracks(), s1_rec);
      ElectronRec(s2.GetTracks(), s2_rec);
    }

    //FinishEvent for clusters after track finder for cluster information about tracks
    for_each(cls.begin(), cls.end(), mem_fn(&TrkClusterFinder::FinishEvent));

    //finish for tracks
    s1.FinishEvent();
    s2.FinishEvent();

    if(write_event_tree) otree.Fill();

  }//event loop

  if(write_event_tree) otree.Write(0, TObject::kOverwrite);

  if( lps_mode == 1 ) {

    cout << "Finalizing..." << endl;
    s1_rec->Finalize();
    s2_rec->Finalize();

    cout << "Exporting..." << endl;
    s1_rec->Export();
    s2_rec->Export();
  }

  out->Close();

  //print counts
  cout << "All done, hit counts:" << endl;
  for_each(hits.begin(), hits.end(), mem_fn(&TrkHitsTransform::PrintCounts));
  cout << "Cluster counts:" << endl;
  for_each(cls.begin(), cls.end(), mem_fn(&TrkClusterFinder::PrintCounts));
  cout << "Track counts:" << endl;
  s1.PrintCounts();
  s2.PrintCounts();

}//Run

//_____________________________________________________________________________
void TagTpix4Reco::ElectronRec(vector<TagTracks::Track>& trk, unique_ptr<LookupProblemSolver>& rec) {

  //track loop
  for(TagTracks::Track& i: trk) {

    //MC loop
    for(const MCParticles::Part& mcp: fMC->GetReadData()) {

      if( mcp.itrk != i.prim_id ) continue;

      //track is paired with the MC particle
      i.has_mcp = kTRUE;

      //track MC particle kinematics
      i.mcp_en = mcp.en;
      i.mcp_theta = mcp.theta;
      i.mcp_phi = mcp.phi;
      i.mcp_Q2 = 2*fBeamEn*mcp.en*(1-TMath::Cos(TMath::Pi()-mcp.theta));

      break;

    }//MC loop

    //reconstruction for the track
    if( !rec ) continue;

    // input tagger quantity
    vector<double> quant = {i.x, i.y, i.theta_x, i.theta_y};

    // reconstructed electron by reference
    vector<double> sol = {0, 0, 0};
    vector<double> sol_err = {0, 0, 0};
    Int_t ipar[2];

    //run the reconstruction
    bool stat = rec->Solve(quant, sol, &sol_err, ipar);
    if( !stat ) continue;

    //cout << "Track:" << endl;
    //cout << sol[0] << " " << sol[1] << " " << sol[2] << endl;
    //cout << i.true_en << " " << i.true_theta << " " << i.true_phi << endl;

    //set the track parameters
    i.is_rec = kTRUE;
    i.rec_en = sol[0];
    i.rec_theta = sol[1];
    i.rec_phi = sol[2];
    i.rec_en_err = sol_err[0];
    i.rec_theta_err = sol_err[1];
    i.rec_phi_err = sol_err[2];
    i.ninp = ipar[0]; // at index 0 from Solve
    i.ilay = ipar[1]; // at index 1

    //reconstructed electron Q^2 by beam energy, electron energy and polar angle
    i.rec_Q2 = 2*fBeamEn*i.rec_en*(1-TMath::Cos(TMath::Pi()-i.rec_theta));

  }//track loop

}//ElectronRec

//_____________________________________________________________________________
void TagTpix4Reco::SetMakeLPS(unique_ptr<LookupProblemSolver>& rec, program_options::variables_map& opt_map) {

  //quantities measured by the taggers
  rec->MakeQuantity("x"); // x position, mm
  rec->MakeQuantity("y"); // y position, mm
  rec->MakeQuantity("theta_x"); // theta_x angle, rad
  rec->MakeQuantity("theta_y"); // theta_y angle, rad

  //initialize the reconstruction
  rec->Initialize(opt_map);

  rec->SetUseFoamCellFinder( opt_map["use_FoamCellFinder"].as<bool>() );
  rec->SetUseMedCellFinder( opt_map["use_MedCellFinder"].as<bool>() );
  rec->SetUseSumCellFinder( opt_map["use_SumCellFinder"].as<bool>() );

}//SetMakeLPS

//_____________________________________________________________________________
void TagTpix4Reco::AddInput(vector<TagTracks::Track>& trk, unique_ptr<LookupProblemSolver>& rec) {

  //event with just one track
  if( trk.size() != 1 ) return;

  //track loop
  for(TagTracks::Track& i: trk) {

    //add input for reconstruction
    vector<double> quant = {i.x, i.y, i.theta_x, i.theta_y};
    vector<double> sol = {fTrueEn, fTrueTheta, fTruePhi};
    rec->AddInput(quant, sol);

  }//track loop

}//AddInput

//_____________________________________________________________________________
void TagTpix4Reco::ImportLPS(unique_ptr<LookupProblemSolver>& rec, program_options::variables_map& opt_map) {

  //initialize the response from trained input
  IOProgramOptions io(opt_map);
  TFile in_resp(io.GetStr("input_resp").c_str(), "read");
  rec->Import(in_resp);

  //error limits in reconstruction (optional)
  if( opt_map["max_resp_err_theta"].as<double>() > 0 ) {
    rec->SetMaxErr(1, opt_map["max_resp_err_theta"].as<double>());
  }
  if( opt_map["max_resp_rel_en"].as<double>() > 0 ) {
    rec->SetMaxRel(0, opt_map["max_resp_rel_en"].as<double>());
  }
  if( opt_map["max_resp_rel_phi"].as<double>() > 0 ) {
    rec->SetMaxRel(2, opt_map["max_resp_rel_phi"].as<double>());
  }

  if( opt_map["min_resp_ninp"].as<int>() > 0 ) {
    rec->SetMinNinp( opt_map["min_resp_ninp"].as<int>() );
  }

}//ImportLPS















