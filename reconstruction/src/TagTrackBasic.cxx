
//_____________________________________________________________________________
//
// Reconstruction task to run TagTrackFindBasic tracker finder
//
// Runs by  run_TagTrackBasic.py  and  run_TagTrackBasicVis.py
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
#include "MCParticles.h"

//local classes
#include "TagTrackBasic.h"
#include "TagTrackFindBasic.h"
#include "TrkPlaneBasicHits.h"
#include "EThetaPhiReco.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
void TagTrackBasic::Run(const char *conf) {

  cout << "TagTrackBasic::Run, " << conf << endl;

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
    ("main.planes_output", program_options::value<bool>(), "Write output for planes")
  ;

  //load the configuration file
  ifstream config(conf);
  program_options::variables_map opt_map;
  program_options::store(program_options::parse_config_file(config, opt), opt_map);

  //geometry
  GeoParser geo;
  LoadIN in(geo);
  in.ReadInput( GetStr(opt_map, "main.geo") );
  //geo.PrintAll();

  //beam energy
  fBeamEn = 18; //GeV

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

  //true event kinematics
  tree.SetBranchAddress("true_el_E", &fTrueEn);
  tree.SetBranchAddress("true_el_theta", &fTrueTheta);
  tree.SetBranchAddress("true_el_phi", &fTruePhi);
  tree.SetBranchAddress("true_Q2", &fTrueQ2);

  //input MC particles
  fMC = new MCParticles::Coll();
  fMC->ConnectInput("mcp", &tree);

  //electron reconstruction
  EThetaPhiReco *s1_rec = 0x0;
  EThetaPhiReco *s2_rec = 0x0;
  if( opt_map.find("main.input_resp") != opt_map.end() ) {
    string input_resp = GetStr(opt_map, "main.input_resp");
    cout << "Response input: " << input_resp << endl;

    //create the response for both taggers
    s1_rec = new EThetaPhiReco("s1");
    s2_rec = new EThetaPhiReco("s2");

    //initialize the response from trained input
    TFile in_resp(input_resp.c_str(), "read");
    s1_rec->Import(&in_resp);
    s2_rec->Import(&in_resp);
  }
  if( s1_rec and opt_map.find("main.min_resp_ninp") != opt_map.end() ) {

    s1_rec->SetMinNinp( opt_map["main.min_resp_ninp"].as<int>() );
    s2_rec->SetMinNinp( opt_map["main.min_resp_ninp"].as<int>() );

  }

  //output
  string outfile = GetStr(opt_map, "main.outfile");
  cout << "Output: " << outfile << endl;
  TFile out(outfile.c_str(), "recreate");

  //interaction (event) output tree
  TTree otree("event", "event");

  TagTrackFindBasic s1("s1");
  TagTrackFindBasic s2("s2");

  s1.SetGeometry(&geo);
  s2.SetGeometry(&geo);

  s1.ConnectHitsInput(&tree);
  s2.ConnectHitsInput(&tree);

  //output for clusters from individual planes
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

  Long64_t ncls_s1=0, ncls_s2=0, ntrk_s1=0, ntrk_s2=0;

  //event loop
  Long64_t nev = tree.GetEntries();
  //Long64_t nev = 4;
  Long64_t iprint = nev/12;
  for(Long64_t iev=0; iev<nev; iev++) {

    tree.GetEntry(iev);

    if( iev > 0 and iev%iprint == 0 ) {
      cout << Form("%.1f", 100.*iev/nev) << "%" << endl;
    }

    //cout << "Next event" << endl;

    fMC->LoadInput();

    s1.ProcessEvent();
    s2.ProcessEvent();

    ElectronRec(&s1, s1_rec);
    ElectronRec(&s2, s2_rec);

    s1.FinishEvent();
    s2.FinishEvent();

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

}//Run

//_____________________________________________________________________________
void TagTrackBasic::ElectronRec(TagTrackFindBasic *tag, EThetaPhiReco *rec) {

  //track loop
  for(TagTrackFindBasic::Track& i: tag->GetTracks()) {

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
    Double_t quant[4]{i.x, i.y, i.theta_x, i.theta_y}; // input tagger quantity
    Double_t rec_en=0, rec_theta=0, rec_phi=0; // reconstructed electron by reference
    Int_t ninp=0;
    Bool_t stat = rec->Reconstruct(quant, rec_en, rec_theta, rec_phi, &ninp); // perform the reconstruction
    if( !stat ) continue;

    //set the track parameters
    i.is_rec = kTRUE;
    i.rec_en = rec_en;
    i.rec_theta = rec_theta;
    i.rec_phi = rec_phi;
    i.ninp = ninp;

    //reconstructed electron Q^2 by beam energy, electron energy and polar angle
    i.rec_Q2 = 2*fBeamEn*i.rec_en*(1-TMath::Cos(TMath::Pi()-i.rec_theta));

  }//tracks loop

}//ElectronRec

//_____________________________________________________________________________
string TagTrackBasic::GetStr(program_options::variables_map& opt_map, std::string par) {

  string res = opt_map[par].as<string>();
  res.erase(remove(res.begin(), res.end(), '\"'), res.end());

  return res;

}//GetStr

//_____________________________________________________________________________
extern "C" {

  TagTrackBasic* make_TagTrackBasic() { return new TagTrackBasic(); }

  void run_TagTrackBasic(void *task, const char *conf) {
    reinterpret_cast<TagTrackBasic*>(task)->Run(conf);
  }

}//extern

