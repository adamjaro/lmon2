
#ifndef CaloBPCHits_h
#define CaloBPCHits_h

#include "DetectorData.h"
#include <map>
#include "Rtypes.h"

namespace CaloBPCHits {

//hit representation
struct Hit {

  Hit() {}
  Hit(Int_t is, Int_t il, Double_t xp, Double_t yp, Double_t zp):
    iscin(is), ilay(il), x(xp), y(yp), z(zp) {}

  Int_t iscin=-1; // scintillator index
  Int_t ilay=-1; // layer index
  Double_t en=0; // hit energy, GeV
  Double_t x=0; // hit position in x, mm
  Double_t y=0; // hit position in y, mm
  Double_t z=0; // hit position in z, mm

};//Hit

//hits collection
class Coll : public DetectorData<Hit, std::map<std::pair<Int_t, Int_t>, Hit>> {

  public:

    Coll() {

      DATA_ADD_UNIT_ATTR( iscin )
      DATA_ADD_UNIT_ATTR( ilay )

      DATA_ADD_UNIT_ATTR( en )
      DATA_ADD_UNIT_ATTR( x )
      DATA_ADD_UNIT_ATTR( y )
      DATA_ADD_UNIT_ATTR( z )

    }//Coll

};//Coll

}//CaloBPCHits

#endif

