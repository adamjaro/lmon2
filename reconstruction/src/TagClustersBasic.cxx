
//C++
#include <iostream>

//ROOT
#include "TTree.h"
#include "TMath.h"

//Geant
#include "G4String.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//local classes
#include "TagClustersBasic.h"
#include "GeoParser.h"

using namespace std;

//_____________________________________________________________________________
TagClustersBasic::TagClustersBasic(std::string nam) {

  fNam = nam;

  fEmin = 0.4; // keV, threshold energy

}//TagClustersBasic

//_____________________________________________________________________________
void TagClustersBasic::SetGeometry(GeoParser *geo) {

  G4String geo_nam = "vac_S1";

  //plane position
  geo->GetOptD(geo_nam, "xpos", fXpos, GeoParser::Unit(mm));
  geo->GetOptD(geo_nam, "ypos", fYpos, GeoParser::Unit(mm));
  geo->GetOptD(geo_nam, "zpos", fZpos, GeoParser::Unit(mm));

  //plane rotation
  G4double theta = 0;
  geo->GetOptD(geo_nam, "theta", theta, GeoParser::Unit(rad));

  G4bool rotate_x = false;
  geo->GetOptB(geo_nam, "rotate_x", rotate_x);

  if( rotate_x ) {
    fThetaX = theta;
  } else {
    fThetaY = theta;
  }

  G4String det_nam = "lowQ2_"+fNam;

  //shift in x position
  fXpos += geo->GetD(det_nam, "xpos");

  cout << "xyz: " << fXpos << " " << fYpos << " " << fZpos << " " << fThetaX << " " << fThetaY << endl;


}//SetGeometry

//_____________________________________________________________________________
void TagClustersBasic::ProcessEvent() {

  //cout << "TagClustersBasic::ProcessEvent" << endl;

  fHits.LoadInput();
  fClusters.ClearEvent();

  //hits loop to init
  for(unsigned long ihit=0; ihit<fHits.GetN(); ihit++) {

    //get the hit
    TrkPlaneBasicHits::Hit& hit = fHits.GetUnit(ihit);

    //transform to plane local coordinates
    hit.Translate(-fXpos, -fYpos, -fZpos);
    hit.RotateXY(-fThetaX, -fThetaY);

    //create status entry for the hit
    unordered_map<unsigned long, bool>::iterator istat = fHitStat.insert( make_pair(ihit, false) ).first;

    //energy threshold
    if( hit.en < fEmin ) continue;

    //mark hit as accepted
    (*istat).second = true;

  }//hits loop to init

  //hits loop to make clusters
  while(GetHitsCount() > 0) {

    //hit with most energy deposition as a seed for the cluster
    unsigned long ih = FindHitEmax();
    fHitStat[ih] = false; // mark the hit as used

    //initialize the cluster to the seed
    Cluster& cls = fClusters.Add(Cluster());
    cls.hits.push_back(ih);

    //MC track and pdg from the seed hit
    const TrkPlaneBasicHits::Hit& seed_hit = fHits.GetUnit( ih );
    cls.itrk = seed_hit.itrk;
    cls.pdg = seed_hit.pdg;
    cls.prim_id = seed_hit.prim_id;

    //adjacent hits for the hit with the most energy
    int nfound = FindAdjHits(ih, cls.hits);

    //cout << nfound << " " << seed_hit.ipix << " " << seed_hit.irow << " ";
    //cout << seed_hit.x << " " << seed_hit.y << " " << seed_hit.z << endl;

    //search for adjacent hits from the first set
    list<unsigned long> adj_hits_sec;
    do {
      nfound = 0;
      adj_hits_sec.clear();

      for(list<unsigned long>::iterator ith = cls.hits.begin(); ith != cls.hits.end(); ith++) {

        nfound += FindAdjHits(*ith, adj_hits_sec);
      }

      //append the next set of adjacent hits to the cluster hits
      for(list<unsigned long>::iterator ith = adj_hits_sec.begin(); ith != adj_hits_sec.end(); ith++) {
        cls.hits.push_back(*ith);
      }

    } while( nfound > 0 );

    //test for uniform gradient in cluster hits
    GradientTest(cls);

  }//hits loop to make clusters

  //clusters loop
  Int_t cls_id = 0;
  for(vector<Cluster>::iterator icls = fClusters.GetStore().begin(); icls != fClusters.GetStore().end(); icls++) {
    Cluster& cls = *icls;
    cls.id = cls_id++; // set the cluster ID

    //number of hits in the cluster
    cls.nhits = cls.hits.size();

    //cluster position by energy-weighted average
    for(list<unsigned long>::iterator ith = cls.hits.begin(); ith != cls.hits.end(); ith++) {
      const TrkPlaneBasicHits::Hit& hit = fHits.GetUnit( *ith );

      cls.en += hit.en;
      cls.x += hit.x*hit.en;
      cls.y += hit.y*hit.en;

      cls.is_prim = cls.is_prim && hit.is_prim;

      cls.sigma_x += hit.en*hit.x*hit.x;
      cls.sigma_y += hit.en*hit.y*hit.y;
    }

    cls.x = cls.x/cls.en;
    cls.y = cls.y/cls.en;

    cls.sigma_x = cls.GetSigma(cls.sigma_x, cls.x);
    cls.sigma_y = cls.GetSigma(cls.sigma_y, cls.y);

  }//clusters loop

  //cout << "clusters:" << endl;

  //cluster mutual distances and status
  for(vector<Cluster>::iterator icls = fClusters.GetStore().begin(); icls != fClusters.GetStore().end(); icls++) {
    Cluster& cls = *icls;

    //inner loop
    for(vector<Cluster>::iterator i2 = fClusters.GetStore().begin(); i2 != fClusters.GetStore().end(); i2++) {
      if( icls == i2 ) continue; // same cluster

      Cluster& cls2 = *i2;

      //distance between the clusters
      Double_t dist = TMath::Sqrt( (cls.x-cls2.x)*(cls.x-cls2.x) + (cls.y-cls2.y)*(cls.y-cls2.y) );

      //update the minimal distance for the cluster
      if( dist < cls.min_dist or cls.min_dist < 0 ) cls.min_dist = dist;

    }//inner loop

    //cluster status
    if( cls.min_dist > 0 and cls.min_dist < fClsMinLimMdist ) cls.stat = kFALSE; // limit on minimal distance to another cluster

    //cout << cls.x << " " << cls.y << " " << cls.en << " " << cls.min_dist << " " << cls.stat << endl;

  }//cluster mutual distances and status

}//ProcessEvent

//_____________________________________________________________________________
void TagClustersBasic::FinishEvent() {

  fClusters.FinishEvent();

}//FinishEvent

//_____________________________________________________________________________
int TagClustersBasic::GetHitsCount() {

  //count hits with active status

  int nhit = 0;
  for(unsigned long ihit=0; ihit<fHitStat.size(); ihit++) {

    if( fHitStat[ihit] ) nhit++;
  }

  return nhit;

}//GetHitsCount

//_____________________________________________________________________________
unsigned long TagClustersBasic::FindHitEmax() {

  //hit with most energy deposition

  //initial values
  unsigned long ih = 0;
  Double_t en = -1.;

  //compare to all other hits
  for(unsigned long ihit=0; ihit<fHits.GetN(); ihit++) {
    if( !fHitStat[ihit] ) continue;

    const TrkPlaneBasicHits::Hit& hit = fHits.GetUnit(ihit);

    if( hit.en < en ) continue;

    //update for hit at larger energy
    en = hit.en;
    ih = ihit;

  }

  return ih;

}//FindHitEmax

//_____________________________________________________________________________
int TagClustersBasic::FindAdjHits(unsigned long ih, list<unsigned long>& adj_hits) {

  //find hits adjacent to a hit at a given ih
  //along pixel and row direction

  const TrkPlaneBasicHits::Hit& hit_start = fHits.GetUnit(ih);
  Int_t ipix_start = hit_start.ipix;
  Int_t irow_start = hit_start.irow;

  //range in pixel and row indices for adjacent hits
  Int_t pmin = ipix_start-1;
  Int_t pmax = ipix_start+1;
  Int_t rmin = irow_start-1;
  Int_t rmax = irow_start+1;

  //counter for found hits
  int nfound = 0;

  for(unsigned long ihit=0; ihit<fHits.GetN(); ihit++) {
    if( !fHitStat[ihit] ) continue;

    //hit indices
    const TrkPlaneBasicHits::Hit& hit = fHits.GetUnit(ihit);
    Int_t ipix = hit.ipix;
    Int_t irow = hit.irow;

    //apply the range for adjacent hit
    if( ipix < pmin or ipix > pmax ) continue;
    if( irow < rmin or irow > rmax ) continue;

    //adjacent hit found
    adj_hits.push_back(ihit);

    //mark the hit as used
    fHitStat[ihit] = false;

    nfound++;

  }

  return nfound;

}//FindAdjHits

//_____________________________________________________________________________
void TagClustersBasic::GradientTest(Cluster& cls) {

  //test hits in cluster for decreasing gradient from the seed
  //and remove those hits which would cause an increase in the gradient

  //cout << fNam << " grad test: " << cls.is_prim << " " << cls.nhits;
  //cout << " " << cls.x << " " << cls.y << endl;

  //seed hit for the cluster
  const TrkPlaneBasicHits::Hit& seed_hit = fHits.GetUnit( cls.hits.front() );

  list<unsigned long> hits_to_remove;

  //hit loop
  for(list<unsigned long>::iterator ith = cls.hits.begin(); ith != cls.hits.end(); ith++) {
    const TrkPlaneBasicHits::Hit& hit = fHits.GetUnit( *ith );

    //distance to the seed
    Int_t dx = seed_hit.ipix-hit.ipix;
    Int_t dy = seed_hit.irow-hit.irow;

    //adjacent hits over two or more pixels
    if( labs(dx) < 2 and labs(dy) < 2 ) continue;

    //direction to the seed in x and y
    Int_t dxy[2] = {0, 0};
    if( dx != 0 ) dxy[0] = dx>0?1:-1;
    if( dy != 0 ) dxy[1] = dy>0?1:-1;

    //cout << "  hit: " << hit.ipix << " " << hit.irow;
    //cout << " " << Form("%.0f", hit.en) << " " << hit.is_prim;
    //cout << " Dx: " << dx << " Dy: " << dy << " dxy: " << dxy[0] << " " << dxy[1] << endl;

    //hit status to satisfy the gradient, true = keep, false = remove from the cluster
    bool grad_stat = false;

    //all hits in the direction to the seed by a unit distance
    for(auto ihs = cls.hits.begin(); ihs != cls.hits.end(); ihs++) {
      const TrkPlaneBasicHits::Hit& hit_to_seed = fHits.GetUnit( *ihs );
      if( labs(hit_to_seed.ipix-hit.ipix) > 1 or labs(hit_to_seed.irow-hit.irow) > 1 ) continue; //unit distance
      if( hit_to_seed.ipix == hit.ipix and hit_to_seed.irow == hit.irow ) continue; //same hit

      //hit is in the direction to seed by dxy
      if( hit.ipix+dxy[0] != hit_to_seed.ipix and hit.irow+dxy[1] != hit_to_seed.irow ) continue;

      //cout << "    hit_to_seed: " << Form("%.0f", hit_to_seed.en) << endl;

      //hit of larger energy in the direction to the seed
      if( hit_to_seed.en > hit.en ) grad_stat = true;
    }

    //cout << "    grad_stat: " << grad_stat << endl;

    //keep the hits satisfying the gradient
    if( grad_stat ) continue;

    //cout << "    to remove: " << Form("%.0f", hit.en) << endl;

    //remove the hit not in gradient and all hits in the direction away from the seed
    for(auto ihs = cls.hits.begin(); ihs != cls.hits.end(); ihs++) {
      const TrkPlaneBasicHits::Hit& hit_from_seed = fHits.GetUnit( *ihs );
      if( hit_from_seed.ipix == hit.ipix and hit_from_seed.irow == hit.irow ) continue; //same hit

      //unit distance about the seed
      if( labs(hit_from_seed.ipix-seed_hit.ipix) < 2
        and labs(hit_from_seed.irow-seed_hit.irow) < 2 ) continue;

      //direction from the hit not in gradient
      Int_t gdx = SignI( hit_from_seed.ipix - hit.ipix );
      Int_t gdy = SignI( hit_from_seed.irow - hit.irow );

      //cout << "    hit_from_seed: " << Form("%.0f", hit_from_seed.en) << " ";
      //cout << gdx << " " << gdy << " " << dxy[0] << " " << dxy[1] << endl;

      //diagonal direction
      if( (labs(dxy[0]) + labs(dxy[1])) > 1 ) {
        if( dxy[0] == gdx ) continue;
        if( dxy[1] == gdy ) continue;
      }

      //direction along x or y
      if( dxy[0] == 0 ) {
        if( gdy == 0 ) continue;
        if( dxy[1] == gdy ) continue;
      }
      if( dxy[1] == 0 ) {
        if( gdx == 0 ) continue;
        if( dxy[0] == gdx ) continue;
      }

      //cout << "    remove hit_from_seed: " << Form("%.0f", hit_from_seed.en) << " ";
      //cout << gdx << " " << gdy << " " << dxy[0] << " " << dxy[1] << endl;

      //remove the hit, direction away from the seed
      hits_to_remove.push_back( *ihs );
    }

    //remove the original hit not in gradient
    hits_to_remove.push_back( *ith );

  }//hit loop

  //removing the hits from the cluster
  for(list<unsigned long>::iterator ith = hits_to_remove.begin(); ith != hits_to_remove.end(); ith++) {

    cls.hits.remove( *ith );
    fHitStat[ *ith ] = true;
  }

}//GradientTest

//_____________________________________________________________________________
Double_t TagClustersBasic::Cluster::GetSigma(Double_t swx2, Double_t pos) {

  //Bevington, Robinson, eq. 4.22, p. 58

  if( hits.size() < 2 ) return 0;

  Double_t sig2 = ((swx2/en) - pos*pos)/(hits.size() - 1);

  if(sig2 > 0) return TMath::Sqrt(sig2);

  return 0;

}//Cluster::GetSigma

//_____________________________________________________________________________
TagClustersBasic::Coll::Coll() {

  AddUnitAttr("_x", fUnitIO.x);
  AddUnitAttr("_y", fUnitIO.y);
  AddUnitAttr("_en", fUnitIO.en);
  AddUnitAttr("_nhits", fUnitIO.nhits);
  AddUnitAttr("_is_prim", fUnitIO.is_prim);
  AddUnitAttr("_sigma_x", fUnitIO.sigma_x);
  AddUnitAttr("_sigma_y", fUnitIO.sigma_y);
  AddUnitAttr("_ntrk", fUnitIO.ntrk);
  AddUnitAttr("_min_dist", fUnitIO.min_dist);

}//Coll





































