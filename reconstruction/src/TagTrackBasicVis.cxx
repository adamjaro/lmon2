
//C++
#include <iostream>
#include "glob.h"

//Boost
#include <boost/program_options.hpp>

//ROOT
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

//Geant
#include "G4Step.hh"
#include "G4String.hh"

//local classes
#include "TagTrackFindBasic.h"
#include "GeoParser.h"
#include "LoadIN.h"
#include "LoadXML.h"
//#include "RefCounter.h"
#include "TagClustersBasic.h"
#include "TagTrackBasicVis.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
TagTrackBasicVis::TagTrackBasicVis(const char *conf): iev(-1), min_ntrk(0),
    min_ncls(0), min_ncnt(0), min_etrk(0), min_sig_trk(0) {

  //configuration file
  program_options::options_description opt("opt");
  opt.add_options()
    ("main.input", program_options::value<string>(), "Analysis input")
    ("main.geo", program_options::value<string>(), "Geometry configuration")
    ("main.outfile", program_options::value<string>(), "Output from the analysis")
    ("main.max_chi2ndf", program_options::value<double>(), "Maximal tracks Chi2/NDF")
    ("main.min_cls_dist", program_options::value<double>(), "Minimal cluster distance")
    ("main.input_resp", program_options::value<string>(), "Input response for reconstruction")
    ("main.planes_output", program_options::value<bool>(), "Write output for planes")
  ;

  //load the configuration file
  ifstream config(conf);
  program_options::variables_map opt_map;
  program_options::store(program_options::parse_config_file(config, opt), opt_map);

  //inputs
  string input = GetStr(opt_map, "main.input");
  glob_t glob_inputs;
  glob(input.c_str(), GLOB_TILDE, NULL, &glob_inputs);

  //input tree
  tree = new TChain("DetectorTree");
  for(size_t i=0; i<glob_inputs.gl_pathc; i++) {
    //cout << "Adding input: " << glob_inputs.gl_pathv[i] << endl;
    tree->Add( glob_inputs.gl_pathv[i] );
  }

  //input MC particles
  //tree->SetBranchAddress("mcp_itrk", &fMCItrk); // track ID for the particle
  //tree->SetBranchAddress("mcp_en", &fMCEn); // particle energy, GeV
  //tree->SetBranchAddress("mcp_theta", &fMCTheta); // particle polar angle, rad
  //tree->SetBranchAddress("mcp_phi", &fMCPhi); // particle azimuthal angle, rad

  //geometry
  string geo_nam = GetStr(opt_map, "main.geo");
  //cout << "Geometry: " << geo_nam << endl;
  GeoParser geo;
  //LoadIN in(geo);
  //in.ReadInput(geo_nam);
  LoadXML xml(geo);
  xml.ReadInput( GetStr(opt_map, "main.geo") );

  //outputs
  //string outfile = GetStr(opt_map, "main.outfile");
  //cout << "Output: " << outfile << endl;
  out = new TFile("vis.root", "recreate");

  //interaction (event) output tree
  otree = new TTree("event", "event");

  //tagger stations
  s1 = new TagTrackFindBasic("s1"); //, tree, &geo, otree
  s2 = new TagTrackFindBasic("s2");

  s1->CreateOutput(otree);
  s2->CreateOutput(otree);

  s1->SetGeometry(&geo);
  s2->SetGeometry(&geo);

  s1->ConnectHitsInput(tree);
  s2->ConnectHitsInput(tree);

  //s1->SetMCParticles(fMCItrk, fMCEn, fMCTheta, fMCPhi);
  //s2->SetMCParticles(fMCItrk, fMCEn, fMCTheta, fMCPhi);

  //track selection for tagger stations
  if( opt_map.find("main.max_chi2ndf") != opt_map.end() ) {

    Double_t max_chi2ndf = opt_map["main.max_chi2ndf"].as<double>();
    //cout << "Using Chi2/NDF = " << max_chi2ndf << endl;

    s1->SetMaxChi2Ndf( max_chi2ndf );
    s2->SetMaxChi2Ndf( max_chi2ndf );
  }
  if( opt_map.find("main.min_cls_dist") != opt_map.end() ) {
    Double_t min_cls_dist = opt_map["main.min_cls_dist"].as<double>();
    cout << "Using min_cls_dist = " << min_cls_dist << endl;

    s1->SetClsLimMdist(min_cls_dist);
    s2->SetClsLimMdist(min_cls_dist);
  }

  //reference counters
  //cnt_s1 = new RefCounter("cnt_s1", tree, &geo, otree);
  //cnt_s2 = new RefCounter("cnt_s2", tree, &geo, otree);

  tag = s1;
  //cnt = cnt_s1;

}//TagTrackBasicVis

//_____________________________________________________________________________
int TagTrackBasicVis::ProcessEvent(bool *stat) {

  if( iev < 0 ) iev = 0;
  if( iev >= tree->GetEntries() ) iev = tree->GetEntries()-1;

  tree->GetEntry(iev);

  //process the event for both taggers
  s1->ProcessEvent();
  s2->ProcessEvent();

  //cnt_s1->ProcessEvent();
  //cnt_s2->ProcessEvent();

  //AssociateMC(*s1, *cnt_s1);
  //AssociateMC(*s2, *cnt_s2);

  s1->FinishEvent();
  s2->FinishEvent();

  //event selection
  if(stat) {

    *stat = true;

    if( GetNumberOfTracks() < min_ntrk ) *stat = false;
    if( tag->GetNumberOfClusters() < min_ncls ) *stat = false;
    if( GetNumberOfRefTracks() < min_ncnt ) *stat = false;
    if( GetNumberOfTracks()-GetNumberOfRefTracks() < min_etrk ) *stat = false;
    if( GetSigTracks() < min_sig_trk ) *stat = false;

  }

  //fill event tree
  //otree->Fill();

  return iev;

}//ProcessEvent

//_____________________________________________________________________________
int TagTrackBasicVis::NextEvent(int di) {

  //load and analyze next event

  //selection loop
  while(true) {

    iev += di;

    bool stat = true;
    ProcessEvent( &stat );

    //range in events
    if( di > 0 and iev >= tree->GetEntries()-1 ) stat = true;
    if( di < 0 and iev <= 0 ) stat = true;

    if(stat) break; // event criteria are satisfied

  }//selection loop

  return iev;

}//NextEvent

//_____________________________________________________________________________
int TagTrackBasicVis::PreviousEvent() {

  return NextEvent(-1);

}//PreviousEvent

//_____________________________________________________________________________
void TagTrackBasicVis::SetMaxChi2ndf(double chi2) {

  //set chi^2 limit for both taggers

  s1->SetMaxChi2Ndf(chi2);
  s2->SetMaxChi2Ndf(chi2);

}//SetMaxChi2ndf

//_____________________________________________________________________________
void TagTrackBasicVis::SetClsLimMdist(double d) {

  //set limit on cluster distance for both taggers

  s1->SetClsLimMdist(d);
  s2->SetClsLimMdist(d);

}//SetClsLimMdist

//_____________________________________________________________________________
int TagTrackBasicVis::GetNumberOfClusters() {

  return tag->GetNumberOfClusters();

}//GetNumberOfClusters

//_____________________________________________________________________________
int TagTrackBasicVis::GetNumberOfClusters(int iplane) {

  //number of clusters

  return tag->GetPlane(iplane)->GetClusters().GetStore().size();

}//GetNumberOfClusters

//_____________________________________________________________________________
void TagTrackBasicVis::GetCluster(int iplane, int icls, double& x, double& y, double& z, double& md) {

  //cluster on a given plane

  TagClustersBasic::Cluster& cls = tag->GetPlane(iplane)->GetClusters().GetStore()[icls];

  //position
  x = cls.x;
  y = cls.y;
  z = tag->GetPlaneZ(iplane);

  //minimal distance to another cluster on the plane
  md = cls.min_dist;

  //tag->GetCluster(iplane, icls, x, y, z, md);

}//GetCluster

//_____________________________________________________________________________
int TagTrackBasicVis::GetNumberOfTracks() {

  //number of reconstructed tracks

  return tag->GetTracks().size();

}//GetNumberOfTracks

//_____________________________________________________________________________
int TagTrackBasicVis::GetSigTracks() {

  //number of signal tracks having itrk == 1

  int nsig = 0;

  for(const auto& i: tag->GetTracks()) {
    if( i.itrk == 1 ) nsig++;
  }

  return nsig;

}//GetSigTracks

//_____________________________________________________________________________
void TagTrackBasicVis::GetTrack(int i, double& x0, double& y0, double& slope_x, double& slope_y, double& chi2, int& itrk) {

  const vector<TagTrackFindBasic::Track>& tracks = tag->GetTracks();

  x0 = tracks[i].x;
  y0 = tracks[i].y;

  slope_x = tracks[i].slope_x;
  slope_y = tracks[i].slope_y;

  chi2 = tracks[i].chi2_xy;

  itrk = tracks[i].itrk;

}//GetTrack

//_____________________________________________________________________________
int TagTrackBasicVis::GetNumberOfRefTracks() {

  //number of reference tracks

  //return cnt->GetTracks().size();

  return 0;

}//GetNumberOfRefTracks

//_____________________________________________________________________________
void TagTrackBasicVis::SetDet(int i) {

  //select the active tagger detector
  if(i == 0) {
    tag = s1;
    //cnt = cnt_s1;
  }
  if(i == 1) {
    tag = s2;
    //cnt = cnt_s2;
  }

}//SetDet

//_____________________________________________________________________________
string TagTrackBasicVis::GetStr(program_options::variables_map& opt_map, std::string par) {

  string res = opt_map[par].as<string>();
  res.erase(remove(res.begin(), res.end(), '\"'), res.end());

  return res;

}//GetStr

//_____________________________________________________________________________
extern "C" {

  //make the instance
  TagTrackBasicVis* make_TagTrackBasicVis(const char *c) { return new TagTrackBasicVis(c); }

  //detector station
  const char* task_TagTrackBasicVis_det_nam(TagTrackBasicVis& t) { return t.GetDetName().c_str(); }
  void task_TagTrackBasicVis_set_det(TagTrackBasicVis& t, int i) { t.SetDet(i); }

  //event navigation
  int task_TagTrackBasicVis_next_event(TagTrackBasicVis& t) { return t.NextEvent(); }
  int task_TagTrackBasicVis_prev_event(TagTrackBasicVis& t) { return t.PreviousEvent(); }
  void task_TagTrackBasicVis_set_event(TagTrackBasicVis& t, int i) { t.SetEvent(i); }
  int task_TagTrackBasicVis_process_event(TagTrackBasicVis& t) { return t.ProcessEvent(); }

  //clusters
  int task_TagTrackBasicVis_ncls(TagTrackBasicVis& t, int i) { return t.GetNumberOfClusters(i); }
  void task_TagTrackBasicVis_cluster(TagTrackBasicVis& t, int iplane, int icls, double& x, double& y, double& z, double& md) {
    return t.GetCluster(iplane, icls, x, y, z, md);
  }

  //tracks
  int task_TagTrackBasicVis_ntrk(TagTrackBasicVis& t) { return t.GetNumberOfTracks(); }
  void task_TagTrackBasicVis_track(TagTrackBasicVis& t, int i, double& x0, double& y0, double& slope_x, double& slope_y,
    double& chi2, int& itrk) {
    return t.GetTrack(i, x0, y0, slope_x, slope_y, chi2, itrk);
  }
  void task_TagTrackBasicVis_set_max_chi2(TagTrackBasicVis& t, double c) { return t.SetMaxChi2ndf(c); }
  double task_TagTrackBasicVis_get_max_chi2(TagTrackBasicVis& t) { return t.GetMaxChi2ndf(); }
  void task_TagTrackBasicVis_set_lim_mdist(TagTrackBasicVis& t, double d) { t.SetClsLimMdist(d); }
  double task_TagTrackBasicVis_get_lim_mdist(TagTrackBasicVis& t) { return t.GetClsLimMdist(); }
  int task_TagTrackBasicVis_ntrk_ref(TagTrackBasicVis& t) { return t.GetNumberOfRefTracks(); }

  //event selection
  void task_TagTrackBasicVis_set_min_ntrk(TagTrackBasicVis& t, int n) { t.SetMinNtrk(n); }
  void task_TagTrackBasicVis_set_min_ncls(TagTrackBasicVis& t, int n) { t.SetMinNcls(n); }
  void task_TagTrackBasicVis_set_min_ncnt(TagTrackBasicVis& t, int n) { t.SetMinNcnt(n); }
  void task_TagTrackBasicVis_set_min_etrk(TagTrackBasicVis& t, int n) { t.SetMinEtrk(n); }
  void task_TagTrackBasicVis_set_min_sig_trk(TagTrackBasicVis& t, int n) { t.SetMinSigTrk(n); }

}





















