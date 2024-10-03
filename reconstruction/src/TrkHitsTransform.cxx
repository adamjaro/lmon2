
//_____________________________________________________________________________
//
// Geometry transformation for TrkPlaneBasicHits
//
//_____________________________________________________________________________

//C++
#include <iostream>

//ROOT
#include "TVector2.h"

//lmon2 base
#include "GeoParser.h"

//local classes
#include "TrkHitsTransform.h"

using namespace std;

//_____________________________________________________________________________
TrkHitsTransform::TrkHitsTransform(string nam, TTree *tree): fNam(nam) {

  fHits.ConnectInput(fNam, tree);

  fEmin = 0.4; // keV, threshold energy

}//TrkHitsTransform

//_____________________________________________________________________________
void TrkHitsTransform::CreateOutput(TTree *tree) {

  fWriteOut=true;

  fHits.CreateOutput(fNam, tree);

  cout << "TrkHitsTransform::CreateOutput, " << fNam << endl;

}//CreateOutput

//_____________________________________________________________________________
void TrkHitsTransform::LoadV6p3_rev3(GeoParser& geo, string section) {

  //load geometry corresponding to lattice V6.3, rev3 layout

  //final tracking plane position, initialize with drift section from B2eR to Q3eR
  TVector2 plane_pos( geo.GetD("vac_B2Q3", "zpos"), geo.GetD("vac_B2Q3", "xpos") );

  //position of tagger section inside the drift section
  TVector2 sec_pos( geo.GetD(section, "zpos"), geo.GetD(section, "xpos") );

  //rotation to the angle of drift section
  sec_pos = sec_pos.Rotate( geo.GetD("vac_B2Q3", "theta") );

  //move with plane position to the tagger section
  plane_pos += sec_pos;

  //local plane position inside the tagger section
  TVector2 lp_pos( geo.GetD(fNam, "zpos"), 0 );

  //rotation to the angle of tagger section inside the drift section
  Double_t lp_theta = geo.GetD("vac_B2Q3", "theta") + geo.GetD(section, "theta");
  lp_pos = lp_pos.Rotate( lp_theta );

  //move with plane position to the local plane position
  plane_pos += lp_pos;

  //set the final plane position and rotation
  fX = plane_pos.Y(); // TVector2 is in z-x plane
  fZ = plane_pos.X(); // TVector2 is in z-x plane
  fTheta = lp_theta;

  cout << "TrkHitsTransform::LoadV6p3_rev3 for " << fNam << ": " << fX << ", "<< fY << ", " << fZ << ", " << fTheta << endl;

}//LoadV6p3_rev3

//_____________________________________________________________________________
void TrkHitsTransform::ProcessEvent() {

  fHits.LoadInput(); // load the hits
  fHits.ClearEvent(); // clear io vectors

  //counter for all hits
  fNall += fHits.GetN();

  //hits loop
  for(TrkPlaneBasicHits::Hit& i: fHits.GetReadData()) {

    //initial flag
    i.rtstat = kFALSE;

    //energy threshold
    if( i.en < fEmin ) continue;

    //valid status for further processing
    i.rtstat = kTRUE;

    //counter for selected hits
    fNsel++;

    //apply the geometry transformation
    i.TranslateXZ(-fX, -fZ);
    i.RotateZX(-fTheta);

    //cout << i.z-0.25 << endl;

    //write copy of the hit to the output if requested
    if(fWriteOut) {
      fHits.ConstructedAt(make_pair(i.ipix, i.irow)) = i;
    }

  }//hits loop

  //output if requested
  if(fWriteOut) fHits.FinishEvent();

}//ProcessEvent

//_____________________________________________________________________________
void TrkHitsTransform::PrintCounts() {

  cout << "TrkHitsTransform::PrintCounts for " << fNam << ", all, sel: ";
  cout << fNall << ", " << fNsel << ", fraction: " << (float)fNsel/(float)fNall << endl;

}//PrintCounts



















