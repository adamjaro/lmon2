
//_____________________________________________________________________________
//
// Reconstruction task to run TagTrackFindBasic tracker finder
// and provide kinematics reconstruction by LookupProblemSolver
//
// Runs by  run_TagTrackRecLPS.py
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
#include "TMath.h"

//lmon2 base
#include "GeoParser.h"
#include "LoadIN.h"
#include "LoadXML.h"
#include "MCParticles.h"

//local classes
#include "TagTrackRecLPS.h"
#include "TagTrackFindBasic.h"
#include "TrkPlaneBasicHits.h"
#include "LookupProblemSolver.h"
#include "CalPWOClusterWavg.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
void TagTrackRecLPS::Run(const char *conf) {

  cout << "TagTrackRecLPS::Run, " << conf << endl;

  //configuration file
  program_options::options_description opt("opt");
  opt.add_options()
    ("main.input", program_options::value<string>(), "Analysis input")
    ("main.geo", program_options::value<string>(), "Geometry configuration")
    ("main.outfile", program_options::value<string>(), "Output from the analysis")
    ("main.max_chi2ndf", program_options::value<double>(), "Maximal tracks Chi2/NDF")
    ("main.min_cls_dist", program_options::value<double>(), "Minimal cluster distance")
    ("main.input_resp", program_options::value<string>(), "Input response for reconstruction")
    ("main.min_resp_ninp", program_options::value<int>(), "Minimal number of inputs for reconstruction")
    ("main.max_resp_rel_en", program_options::value<double>()->default_value(-1), "Maximal relative error in energy")
    ("main.max_resp_rel_phi", program_options::value<double>()->default_value(-1), "Maximal relative error in phi")
    ("main.max_resp_err_theta", program_options::value<double>()->default_value(-1), "Maximal error in theta (rad)")
    ("main.planes_output", program_options::value<bool>(), "Write output for planes")
    ("main.cal_s1_name", program_options::value<string>()->default_value("lowQ2_s1_pwo"), "Calorimeter name in S1")
  ;

  //load the configuration file
  ifstream config(conf);
  program_options::variables_map opt_map;
  program_options::store(program_options::parse_config_file(config, opt), opt_map);

  //string from program_options
  auto GetStr = [&opt_map](string par) {
    string res = opt_map[par].as<string>();
    res.erase(remove(res.begin(), res.end(), '\"'), res.end());
    return res;
  };

  //geometry
  GeoParser geo;
  LoadXML xml(geo);
  xml.ReadInput( GetStr("main.geo") );
  //geo.PrintAll();

  //beam energy
  fBeamEn = 18; //GeV

  //inputs
  string input = GetStr("main.input");
  glob_t glob_inputs;
  glob(input.c_str(), GLOB_TILDE, NULL, &glob_inputs);

  //input tree
  TChain tree("DetectorTree");
  for(size_t i=0; i<glob_inputs.gl_pathc; i++) {
    cout << "Adding input: " << glob_inputs.gl_pathv[i] << endl;
    tree.Add( glob_inputs.gl_pathv[i] );
  }

  //true event kinematics
  tree.SetBranchAddress("true_el_E", &fTrueEn);
  tree.SetBranchAddress("true_el_theta", &fTrueTheta);
  tree.SetBranchAddress("true_el_phi", &fTruePhi);
  tree.SetBranchAddress("true_Q2", &fTrueQ2);
  tree.SetBranchAddress("true_x", &fTrueX);

  //input MC particles
  fMC = new MCParticles::Coll();
  fMC->ConnectInput("mcp", &tree);

  //electron reconstruction
  unique_ptr<LookupProblemSolver> s1_rec, s2_rec;

  if( opt_map.find("main.input_resp") != opt_map.end() ) {
    string input_resp = GetStr("main.input_resp");
    cout << "Response input: " << input_resp << endl;

    //create the response for both taggers
    s1_rec = make_unique<LookupProblemSolver>(3, "s1");
    s2_rec = make_unique<LookupProblemSolver>(3, "s2");

    //initialize the response from trained input
    TFile in_resp(input_resp.c_str(), "read");
    s1_rec->Import(in_resp);
    s2_rec->Import(in_resp);

    //error limits in reconstruction (optional)
    if( opt_map["main.max_resp_err_theta"].as<double>() > 0 ) {
      s1_rec->SetMaxErr(1, opt_map["main.max_resp_err_theta"].as<double>());
      s2_rec->SetMaxErr(1, opt_map["main.max_resp_err_theta"].as<double>());
    }
    if( opt_map["main.max_resp_rel_en"].as<double>() > 0 ) {
      s1_rec->SetMaxRel(0, opt_map["main.max_resp_rel_en"].as<double>());
      s2_rec->SetMaxRel(0, opt_map["main.max_resp_rel_en"].as<double>());
    }
    if( opt_map["main.max_resp_rel_phi"].as<double>() > 0 ) {
      s1_rec->SetMaxRel(2, opt_map["main.max_resp_rel_phi"].as<double>());
      s2_rec->SetMaxRel(2, opt_map["main.max_resp_rel_phi"].as<double>());
    }
  }

  if( s1_rec and opt_map.find("main.min_resp_ninp") != opt_map.end() ) {

    s1_rec->SetMinNinp( opt_map["main.min_resp_ninp"].as<int>() );
    s2_rec->SetMinNinp( opt_map["main.min_resp_ninp"].as<int>() );
  }

  //output
  string outfile = GetStr("main.outfile");
  cout << "Output: " << outfile << endl;
  TFile out(outfile.c_str(), "recreate");

  //interaction (event) output tree
  TTree otree("event", "event");

  //true event kinematics in output tree
  otree.Branch("true_el_E", &fTrueEn, "true_el_E/D");
  otree.Branch("true_el_theta", &fTrueTheta, "true_el_theta/D");
  otree.Branch("true_el_phi", &fTruePhi, "true_el_phi/D");
  otree.Branch("true_Q2", &fTrueQ2, "true_Q2/D");
  otree.Branch("true_x", &fTrueX, "true_x/D");

  //event counters
  Int_t evt_s1_ntrk=0, evt_s2_ntrk=0;
  otree.Branch("s1_ntrk", &evt_s1_ntrk, "s1_ntrk/I");
  otree.Branch("s2_ntrk", &evt_s2_ntrk, "s2_ntrk/I");

  //tracking
  TagTrackFindBasic s1("s1");
  TagTrackFindBasic s2("s2");

  s1.SetGeometry(&geo);
  s2.SetGeometry(&geo);

  s1.ConnectHitsInput(&tree);
  s2.ConnectHitsInput(&tree);

  //output for tracking clusters from individual planes
  bool planes_output = true;
  if( opt_map.find("main.planes_output") != opt_map.end() ) {
    planes_output = opt_map["main.planes_output"].as<bool>();
  }

  s1.CreateOutput(&otree, planes_output);
  s2.CreateOutput(&otree, planes_output);

  //track selection criteria
  if( opt_map.find("main.max_chi2ndf") != opt_map.end() ) {
    Double_t max_chi2ndf = opt_map["main.max_chi2ndf"].as<double>();
    cout << "Using Chi2/NDF = " << max_chi2ndf << endl;

    s1.SetMaxChi2Ndf( max_chi2ndf );
    s2.SetMaxChi2Ndf( max_chi2ndf );
  }
  if( opt_map.find("main.min_cls_dist") != opt_map.end() ) {
    Double_t min_cls_dist = opt_map["main.min_cls_dist"].as<double>();
    cout << "Using min_cls_dist = " << min_cls_dist << endl;

    s1.SetClsLimMdist(min_cls_dist);
    s2.SetClsLimMdist(min_cls_dist);
  }

  //calorimeter clusters
  string cal_s1_name = GetStr("main.cal_s1_name"); // name for calorimeter volume
  CalPWOClusterWavg cls_s1(cal_s1_name); // name as in the geometry
  cls_s1.ConnectInput(&tree);
  cls_s1.SetGeometry("vac_S1", cal_s1_name, &geo);
  cls_s1.CreateOutput(&otree);

  //total reconstruction counters
  Long64_t ncls_s1=0, ncls_s2=0, ntrk_s1=0, ntrk_s2=0; // tracking counters
  Long64_t ncal_cls_s1=0; // calorimeter counters

  //event loop
  Long64_t nev = tree.GetEntries();
  //Long64_t nev = 320;
  Long64_t iprint = nev/12;
  for(Long64_t iev=0; iev<nev; iev++) {

    tree.GetEntry(iev);

    if( iev > 0 and iev%iprint == 0 ) {
      cout << Form("%.1f", 100.*iev/nev) << "%" << endl;
    }

    //cout << "Next event" << endl;

    //load the mc
    fMC->LoadInput();

    //run the reconstruction
    s1.ProcessEvent();
    s2.ProcessEvent();

    ElectronRec(s1, s1_rec);
    ElectronRec(s2, s2_rec);

    //calorimeter clusters
    cls_s1.ProcessEvent();

    //track and calorimeter cluster matching
    //TrackCalMatch(&s1, &cls_s1);

    s1.FinishEvent();
    s2.FinishEvent();

    //set the event counters
    evt_s1_ntrk = s1.GetTracks().size();
    evt_s2_ntrk = s2.GetTracks().size();
    ncal_cls_s1 += cls_s1.GetClusters().size();

    //fill event tree
    otree.Fill();

    //counters
    ncls_s1 += s1.GetNumberOfClusters();
    ncls_s2 += s2.GetNumberOfClusters();
    ntrk_s1 += s1.GetTracks().size();
    ntrk_s2 += s2.GetTracks().size();

  }//event loop

  otree.Write(0, TObject::kOverwrite);

  out.Close();

  cout << "Events: " << nev << endl;
  cout << "Clusters, s1: " << ncls_s1 << ", s2: " << ncls_s2 << endl;
  cout << "Tracks, s1: " << ntrk_s1 << ", s2: " << ntrk_s2 << endl;
  cout << "Calorimeter clusters, s1: " << ncal_cls_s1 << endl;

}//Run

//_____________________________________________________________________________
void TagTrackRecLPS::ElectronRec(TagTrackFindBasic& tag, std::unique_ptr<LookupProblemSolver>& rec) {

  if( !rec ) return;

  //cout << "TagTrackRecLPS::ElectronRec" << endl;

  //track loop
  for(TagTrackFindBasic::Track& i: tag.GetTracks()) {

    //true event kinematics
    i.true_en = fTrueEn;
    i.true_theta = fTrueTheta;
    i.true_phi = fTruePhi;
    i.true_Q2 = fTrueQ2;

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

    // input tagger quantity
    vector<double> quant = {i.x, i.y, i.theta_x, i.theta_y};
    vector<double> sol = {0, 0, 0}; // reconstructed electron by reference
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

  }//tracks loop

}//ElectronRec


/*

//_____________________________________________________________________________
void TagTrackBasic::TrackCalMatch(TagTrackFindBasic *tag, CalPWOClusterWavg *cal) {

  //track and calorimeter cluster matching, version for CalPWOClusterWavg
  //which is assuming one track (and one calorimeter cluster) per event

  //event with cluster
  if( cal->GetClusters().size() != 1 ) return;

  //get the calorimeter cluster
  const CaloCluster::Cls& cls = cal->GetClusters().at(0);

  //cout << "cluster: " << cls.x << " " << cls.y << endl;

  //track loop
  for(TagTrackFindBasic::Track& i: tag->GetTracks()) {

    //set the matched calorimeter cluster
    i.has_cal = kTRUE;
    i.cal_x = cls.x;
    i.cal_y = cls.y;
    i.cal_en = cls.en;

    //extrapolate the track to calorimeter position
    Double_t xe=0, ye=0;
    i.ExtrapolateZ(-650, xe, ye); // distance from middle tracker to calorimeter front
    i.cal_ext_x = xe+3.9; // offset to calorimeter coordinates
    i.cal_ext_y = ye;

    //cout << "track: " << i.cal_x << " " << i.cal_y << " " << i.cal_ext_x << " " << i.cal_ext_y << endl;

  }//track loop

}//TrackCalMatch

*/











