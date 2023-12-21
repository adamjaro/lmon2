
#ifndef CalPWOHits_h
#define CalPWOHits_h

#include "DetectorData.h"
#include <unordered_map>
#include "Rtypes.h"

namespace CalPWOHits {

//hit representation
struct Hit {

  //geometry transformations
  void Translate(Double_t xp, Double_t yp, Double_t zp);
  void RotateXY(Double_t tx, Double_t ty);

  //hit members to appear in fUnitIO
  Int_t cell_id=0; // cell ID in the module
  Double_t x=0; // cell position in x, mm
  Double_t y=0; // cell position in y, mm
  Double_t z=0; // cell position in z, mm
  Double_t en=0; // hit energy, GeV
  Int_t prim_id=-1; // ID of primary particle associated with the hit

  std::unordered_map<Int_t, Double_t> prim_energy; // energy by individual primaries, transient member

  Hit() {}
  Hit(Int_t i, Double_t xp, Double_t yp, Double_t zp): cell_id(i), x(xp), y(yp), z(zp) {}

};//Hit

//hits collection
class Coll : public DetectorData<Hit, std::unordered_map<Int_t, Hit>> {

  public:

  Coll();

  void FinishEvent();

};//Coll

}//CalPWOHits

#endif

