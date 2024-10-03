
//_____________________________________________________________________________
//
// Hits for TrkPlaneBasic
//
//_____________________________________________________________________________

//ROOT
#include "TVector3.h"

//local classes
#include "TrkPlaneBasicHits.h"

using namespace std;

//_____________________________________________________________________________
TrkPlaneBasicHits::Hit::Hit(Int_t ip, Int_t ir, Double_t xp, Double_t yp, Double_t zp, Int_t it, Int_t pd, Bool_t prim, Int_t pid):
  ipix(ip), irow(ir), x(xp), y(yp), z(zp), en(0), itrk(it), pdg(pd), is_prim(prim), prim_id(pid) {

}//Hit

//_____________________________________________________________________________
void TrkPlaneBasicHits::Hit::Translate(Double_t xp, Double_t yp, Double_t zp) {

  //hit translation in all coordinates
  x += xp;
  y += yp;
  z += zp;

}//Hit::Translate

//_____________________________________________________________________________
void TrkPlaneBasicHits::Hit::TranslateXZ(Double_t xp, Double_t zp) {

  //hit translation in x and z
  x += xp;
  z += zp;

}//Hit::TranslateXZ

//_____________________________________________________________________________
void TrkPlaneBasicHits::Hit::RotateZX(Double_t theta) {

  //hit rotation in z-x plane
  TVector2 zx = TVector2(z, x).Rotate(theta);

  //the TVector2 is given in z-x plane
  x = zx.Y();
  z = zx.X();

}//RotateZX

//_____________________________________________________________________________
void TrkPlaneBasicHits::Hit::RotateXY(Double_t tx, Double_t ty) {

  TVector3 pos(x, y, z);
  pos.RotateY(ty);
  pos.RotateX(tx);

  x = pos.X();
  y = pos.Y();
  z = pos.Z();

}//Hit::RotateXY

//_____________________________________________________________________________
TrkPlaneBasicHits::Coll::Coll() {

  DATA_ADD_UNIT_ATTR( ipix )
  DATA_ADD_UNIT_ATTR( irow )

  DATA_ADD_UNIT_ATTR( x )
  DATA_ADD_UNIT_ATTR( y )
  DATA_ADD_UNIT_ATTR( z )
  DATA_ADD_UNIT_ATTR( en )

  DATA_ADD_UNIT_ATTR( itrk )
  DATA_ADD_UNIT_ATTR( pdg )
  DATA_ADD_UNIT_ATTR( is_prim )
  DATA_ADD_UNIT_ATTR( prim_id )

}//Coll

