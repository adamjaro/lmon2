
//_____________________________________________________________________________
//
// Track finder for timepix tracker
//
//_____________________________________________________________________________

//C++
#include <iostream>

//local classes
#include "TrkClusterFinder.h"
#include "GeoParser.h"
#include "TagTrackFinder.h"

using namespace std;

//_____________________________________________________________________________
TagTrackFinder::TagTrackFinder(std::string nam, vector<shared_ptr<TrkClusterFinder>>& cls, size_t iofs):
    fNam(nam), fPlanes(4) {

  cout << "TagTrackFinder::TagTrackFinder: " << fNam << endl;

  //planes used in the finder starting at offset position iofs
  for(size_t i=iofs; i<iofs+4; i++) {
    fPlanes[i-iofs] = cls[i];
    cout << fPlanes[i-iofs]->GetName() << endl;
  }

}//TagTrackFinder

//_____________________________________________________________________________
void TagTrackFinder::LoadV6p3_rev3(GeoParser& geo, std::string p1, std::string p2, std::string par) {

  //load geometry corresponding to lattice V6.3, rev3 layout

  //distance between planes along z
  fL = geo.GetD(p2, par) - geo.GetD(p1, par);

  //local z positions for planes, mm
  fZ[0] = (-3./2)*fL;
  fZ[1] = (-1./2)*fL;
  fZ[2] = (1./2)*fL;
  fZ[3] = (3./2)*fL;

  cout << "TagTrackFinder::LoadV6p3_rev3, " << fNam << ", fL: " << fL << endl;

}//LoadV6p3_rev3

//_____________________________________________________________________________
void TagTrackFinder::ProcessEvent() {

  fTracks.ClearEvent();

  //clusters in planes
  vector<TrkPlaneClusters::Cluster>& cls1 = fPlanes[0]->GetClusters().GetStore();
  vector<TrkPlaneClusters::Cluster>& cls2 = fPlanes[1]->GetClusters().GetStore();
  vector<TrkPlaneClusters::Cluster>& cls3 = fPlanes[2]->GetClusters().GetStore();
  vector<TrkPlaneClusters::Cluster>& cls4 = fPlanes[3]->GetClusters().GetStore();

  //plane 1
  for(TrkPlaneClusters::Cluster& c1: cls1) {
    c1.iplane = 1;
    if( !c1.stat ) continue;

    //plane 2
    for(TrkPlaneClusters::Cluster& c2: cls2) {
      c2.iplane = 2;
      if( !c2.stat ) continue;

      //plane 3
      for(TrkPlaneClusters::Cluster& c3: cls3) {
        c3.iplane = 3;
        if( !c3.stat ) continue;

        //plane 4
        for(TrkPlaneClusters::Cluster& c4: cls4) {
          c4.iplane = 4;
          if( !c4.stat ) continue;

          //make the track from the clusters

          //track points in x and y
          Double_t x[] = {c1.x, c2.x, c3.x, c4.x};
          Double_t y[] = {c1.y, c2.y, c3.y, c4.y};

          //track parameters in x and y for the track candidate
          TagTracks::Track init_trk;
          MakeTrack(x, init_trk.x, init_trk.slope_x, init_trk.theta_x, init_trk.chi2_x); // track init in x
          MakeTrack(y, init_trk.y, init_trk.slope_y, init_trk.theta_y, init_trk.chi2_y); // track init in y

          //maximal track reduced chi2 in xy plane
          if( TrackChi2(x, y, init_trk) > 4.*fChi2ndfMax ) continue; // 4 degrees of freedom in the xy plane

          //cout << "Track " << init_trk.chi2_xy << endl;

          //the track is selected, add it for the event
          TagTracks::Track& trk = fTracks.Add(init_trk);

          //evaluate clusters making the track
          TrkPlaneClusters::Cluster *cls[] = {&c1, &c2, &c3, &c4};
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

  fTracks.FinishEvent();

}//ProcessEvent

//_____________________________________________________________________________
void TagTrackFinder::MakeTrack(Double_t *x, Double_t& pos, Double_t& slope, Double_t& theta, Double_t& chi2) {

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
Double_t TagTrackFinder::TrackChi2(Double_t *x, Double_t *y, TagTracks::Track& trk) {

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
template<size_t N> void TagTrackFinder::ClusterAnalysis(TrkPlaneClusters::Cluster* (&cls)[N], TagTracks::Track& trk) {

  unordered_set<Int_t> itrk; // MC track indices from the clusters

  //cluster loop
  for(size_t icls=0; icls<N; icls++) {

    //increment track counts for the cluster
    cls[icls]->ntrk += 1;

    itrk.insert( cls[icls]->itrk ); // add the MC track index

  }//cluster loop

  trk.num_diff_itrk = itrk.size(); // number of unique MC track indices

}//ClusterAnalysis

































