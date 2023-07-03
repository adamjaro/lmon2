
#ifndef ParticleCounterHits_h
#define ParticleCounterHits_h

// hits for ParticleCounter

#include "DetectorData.h"

namespace ParticleCounterHits {

//hit representation
struct Hit {

  Int_t pdg; // particle pdg
  Double_t en; // hit energy, GeV
  Double_t x; // hit position in x, mm
  Double_t y; // hit position in y, mm
  Double_t z; // hit position in z, mm
  Int_t parentID; // parent ID for track in hit
  Int_t itrk; // track index
  Bool_t is_prim; // hit by primary particle

};//Hit

//hits collection
class Coll : public DetectorData<Hit> {

  public:

    Coll() {

      //hits memory representation
      AddUnitAttr("_HitPdg", fUnitIO.pdg);
      AddUnitAttr("_HitEn", fUnitIO.en);
      AddUnitAttr("_HitX", fUnitIO.x);
      AddUnitAttr("_HitY", fUnitIO.y);
      AddUnitAttr("_HitZ", fUnitIO.z);
      AddUnitAttr("_HitParentID", fUnitIO.parentID);
      AddUnitAttr("_HitItrk", fUnitIO.itrk);
      AddUnitAttr("_HitPrim", fUnitIO.is_prim);
    }
};//Coll

}//ParticleCounterHits

#endif
















