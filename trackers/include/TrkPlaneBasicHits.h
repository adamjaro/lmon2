
#ifndef TrkPlaneBasicHits_h
#define TrkPlaneBasicHits_h

#include "DetectorData.h"
#include <map>
#include "Rtypes.h"

namespace TrkPlaneBasicHits {

//hit representation
struct Hit {

  Hit() {}
  Hit(Int_t ip, Int_t ir, Double_t xp, Double_t yp, Double_t zp, Int_t it, Int_t pd, Bool_t prim, Int_t pid);

  void Translate(Double_t xp, Double_t yp, Double_t zp);
  void TranslateXZ(Double_t xp, Double_t zp);
  void RotateXY(Double_t tx, Double_t ty);
  void RotateZX(Double_t theta);

  Int_t ipix=0; // pixel index in the row
  Int_t irow=0; // row index in the layer
  Double_t x=0; // hit position in x, mm
  Double_t y=0; // hit position in y, mm
  Double_t z=0; // hit position in z, mm
  Double_t en=0; // hit energy, keV
  Int_t itrk=0; // track index
  Int_t pdg=0; // track PDG code
  Bool_t is_prim=0; // hit by primary particle
  Int_t prim_id=0; // ID of primary particle associated with the hit
  Bool_t rtstat=0; // run-time status used in cluster finder

};//Hit

//hits collection
class Coll : public DetectorData<Hit, std::map<std::pair<Int_t, Int_t>, Hit>> {

  public:

    Coll();

};//Coll

}

#endif

