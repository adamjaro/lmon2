
//_____________________________________________________________________________
//
// Track finder for clusters by TagClustersBasic
//
//_____________________________________________________________________________

//C++
#include <vector>

//ROOT
#include "TTree.h"
#include "TMath.h"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//local classes
#include "TagTrackFindBasic.h"
//#include "TagClustersBasic.h"
#include "GeoParser.h"

using namespace std;

//_____________________________________________________________________________
TagTrackFindBasic::TagTrackFindBasic(string nam) {

  fNam = nam;

  //planes for the station
  fPlanes.push_back( new TagClustersBasic(nam+"_1") );
  fPlanes.push_back( new TagClustersBasic(nam+"_2") );
  fPlanes.push_back( new TagClustersBasic(nam+"_3") );
  fPlanes.push_back( new TagClustersBasic(nam+"_4") );

}//TagTrackFindBasic

//_____________________________________________________________________________
void TagTrackFindBasic::SetGeometry(GeoParser *geo) {

  //individual planes
  for(TagClustersBasic *i: fPlanes) i->SetGeometry(geo);

  //plane spacing, mm
  fL = geo->GetD("lowQ2_"+fNam+"_2", "zpos") - geo->GetD("lowQ2_"+fNam+"_1", "zpos");
  //cout << "L: " << fL << endl;

  //local z positions for planes, mm
  fZ[0] = (-3./2)*fL;
  fZ[1] = (-1./2)*fL;
  fZ[2] = (1./2)*fL;
  fZ[3] = (3./2)*fL;

}//SetGeometry

//_____________________________________________________________________________
void TagTrackFindBasic::CreateOutput(TTree *tree, bool create_planes) {

  //individual planes
  if(create_planes) {
    for(TagClustersBasic *i: fPlanes) i->CreateOutput(tree);
  }

  //tracks
  fTracks.CreateOutput(fNam+"_tracks", tree);

}//CreateOutput

//_____________________________________________________________________________
void TagTrackFindBasic::ProcessEvent() {

  fTracks.ClearEvent();

  //load hits and make clusters for individual planes
  for_each(fPlanes.begin(), fPlanes.end(), mem_fn( &TagClustersBasic::ProcessEvent ));

  //clusters in planes
  vector<TagClustersBasic::Cluster>& cls1 = fPlanes[0]->GetClusters().GetStore();
  vector<TagClustersBasic::Cluster>& cls2 = fPlanes[1]->GetClusters().GetStore();
  vector<TagClustersBasic::Cluster>& cls3 = fPlanes[2]->GetClusters().GetStore();
  vector<TagClustersBasic::Cluster>& cls4 = fPlanes[3]->GetClusters().GetStore();

  //plane 1
  for(TagClustersBasic::Cluster& c1: cls1) {
    c1.iplane = 1;
    if( !c1.stat ) continue;

    //plane 2
    for(TagClustersBasic::Cluster& c2: cls2) {
      c2.iplane = 2;
      if( !c2.stat ) continue;

      //plane 3
      for(TagClustersBasic::Cluster& c3: cls3) {
        c3.iplane = 3;
        if( !c3.stat ) continue;

        //plane 4
        for(TagClustersBasic::Cluster& c4: cls4) {
          c4.iplane = 4;
          if( !c4.stat ) continue;

          //make the track from the clusters

          //track points in x and y
          Double_t x[] = {c1.x, c2.x, c3.x, c4.x};
          Double_t y[] = {c1.y, c2.y, c3.y, c4.y};

          //track parameters in x and y for the track candidate
          Track init_trk;
          MakeTrack(x, init_trk.x, init_trk.slope_x, init_trk.theta_x, init_trk.chi2_x);
          MakeTrack(y, init_trk.y, init_trk.slope_y, init_trk.theta_y, init_trk.chi2_y);

          //maximal tracks reduced chi2
          //if( init_trk.chi2_x > 2.*fChi2ndfMax ) continue; // 2 degrees of freedom
          //if( init_trk.chi2_y > 2.*fChi2ndfMax ) continue; // 2 degrees of freedom

          //maximal track reduced chi2 in the xy plane
          if( TrackChi2(x, y, init_trk) > 4.*fChi2ndfMax ) continue; // 4 degrees of freedom in the xy plane

          //the track is selected, add it for the event
          Track& trk = fTracks.Add(init_trk);

          //evaluate clusters making the track
          TagClustersBasic::Cluster *cls[] = {&c1, &c2, &c3, &c4};
          trk.cls.assign(cls, cls+4); // clusters in the track
          ClusterAnalysis(cls, trk);

          //track for primary particle
          trk.is_prim = c1.is_prim and c2.is_prim and c3.is_prim and c4.is_prim;

          //MC particle corresponding to the track, set for all clusters from the same MC particle
          if( (c1.itrk == c2.itrk) and (c2.itrk == c3.itrk) and (c3.itrk == c4.itrk) ) {

            trk.itrk = c1.itrk;
            trk.pdg = c1.pdg;
          }

          //primary particle ID, set to the track when the primary ID is the same for all clusters
          if( (c1.prim_id == c2.prim_id) and (c2.prim_id == c3.prim_id) and (c3.prim_id == c4.prim_id) ) {

            trk.prim_id = c1.prim_id;
          }

        }//plane 4
      }//plane 3
    }//plane 2
  }//plane 1

}//ProcessEvent

//_____________________________________________________________________________
void TagTrackFindBasic::FinishEvent() {

  fTracks.FinishEvent();

  //finish for planes
  for_each(fPlanes.begin(), fPlanes.end(), mem_fn( &TagClustersBasic::FinishEvent ));

}//FinishEvent

//_____________________________________________________________________________
void TagTrackFindBasic::SetClsLimMdist(Double_t d) {

  for(TagClustersBasic *c: fPlanes) {
    c->SetLimMdist(d);
  }

}//SetClsLimMdist

//_____________________________________________________________________________
void TagTrackFindBasic::MakeTrack(Double_t *x, Double_t& pos, Double_t& slope, Double_t& theta, Double_t& chi2) {

  //track position (mm) and slope
  pos = (x[0]+x[1]+x[2]+x[3])/4.;

  slope = (-3*x[0]-x[1]+x[2]+3*x[3])/(10*fL);

  //angle, rad
  theta = TMath::ATan(slope);

  //chi^2
  chi2 = 0;
  for(int i=0; i<4; i++) {
    chi2 += (pos + slope*fZ[i] - x[i])*(pos + slope*fZ[i] - x[i]);
  }

  //cout << "Track: " << pos << " " << slope << " " << theta << " " << chi2 << endl;

}//MakeTrack

//_____________________________________________________________________________
Double_t TagTrackFindBasic::TrackChi2(Double_t *x, Double_t *y, Track& trk) {

  //calculate track chi^2 for its points
  Double_t chi2_xy = 0;

  //points loop
  for(int i=0; i<4; i++) {

    //track position in x and y
    Double_t track_x_i = trk.x + trk.slope_x*fZ[i];
    Double_t track_y_i = trk.y + trk.slope_y*fZ[i];

    //square distance between the track and measured point along x and y
    Double_t dx2 = (track_x_i-x[i])*(track_x_i-x[i]);
    Double_t dy2 = (track_y_i-y[i])*(track_y_i-y[i]);

    //add the square distance in the xy plane to the chi^2
    chi2_xy += dx2 + dy2;

    //cout << i << ": " << track_x_i-x[i] << " " << track_y_i-y[i] << " ";
    //cout << i << ": " << dx2 << " " << dy2 << " ";

  }//points loop

  //set the chi^2 for the track
  trk.chi2_xy = chi2_xy;

  return chi2_xy;

}//TrackChi2

//_____________________________________________________________________________
template<size_t N> void TagTrackFindBasic::ClusterAnalysis(TagClustersBasic::Cluster* (&cls)[N], Track& trk) {

  //cout << "Clusters:";

  unordered_set<Int_t> itrk; // MC track indices from the clusters

  //cluster loop
  for(size_t icls=0; icls<N; icls++) {

    //increment track counts for the cluster
    cls[icls]->ntrk += 1;

    itrk.insert( cls[icls]->itrk ); // add the MC track index

    //cout << " " << cls[icls]->iplane << " " << cls[icls]->id;
    //cout << " " << cls[icls]->itrk;

  }//cluster loop

  trk.num_diff_itrk = itrk.size(); // number of unique MC track indices

  //cout << " " << itrk.size();
  //cout << endl;

}//ClusterAnalysis

//_____________________________________________________________________________
int TagTrackFindBasic::GetNumberOfClusters() {

  int n = 0;

  for(TagClustersBasic *i: fPlanes) {
    n += i->GetClusters().GetStore().size();
  }

  return n;

}//GetNumberOfClusters

//_____________________________________________________________________________
TagTrackFindBasic::Coll::Coll() {

  //data members of Track to be written in the output from Coll of Track objects
  DATA_ADD_UNIT_ATTR( x )
  DATA_ADD_UNIT_ATTR( y )
  DATA_ADD_UNIT_ATTR( slope_x )
  DATA_ADD_UNIT_ATTR( slope_y )
  DATA_ADD_UNIT_ATTR( theta_x )
  DATA_ADD_UNIT_ATTR( theta_y )
  DATA_ADD_UNIT_ATTR( chi2_x )
  DATA_ADD_UNIT_ATTR( chi2_y )
  DATA_ADD_UNIT_ATTR( chi2_xy )
  DATA_ADD_UNIT_ATTR( is_prim )
  DATA_ADD_UNIT_ATTR( itrk )
  DATA_ADD_UNIT_ATTR( pdg )
  DATA_ADD_UNIT_ATTR( prim_id )

  DATA_ADD_UNIT_ATTR( is_rec )
  DATA_ADD_UNIT_ATTR( rec_en )
  DATA_ADD_UNIT_ATTR( rec_theta )
  DATA_ADD_UNIT_ATTR( rec_phi )
  DATA_ADD_UNIT_ATTR( rec_Q2 )
  DATA_ADD_UNIT_ATTR( ninp )

  DATA_ADD_UNIT_ATTR( has_mcp )
  DATA_ADD_UNIT_ATTR( mcp_en )
  DATA_ADD_UNIT_ATTR( mcp_theta )
  DATA_ADD_UNIT_ATTR( mcp_phi )
  DATA_ADD_UNIT_ATTR( mcp_Q2 )

  DATA_ADD_UNIT_ATTR( true_en )
  DATA_ADD_UNIT_ATTR( true_theta )
  DATA_ADD_UNIT_ATTR( true_phi )
  DATA_ADD_UNIT_ATTR( true_Q2 )

}//Coll



























