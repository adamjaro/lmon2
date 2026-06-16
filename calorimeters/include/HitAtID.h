
#ifndef HitAtID_h
#define HitAtID_h

#include <unordered_map>

#include "DetectorData.h"
#include "Rtypes.h"

namespace HitAtID {

//hit representation
struct Hit {

  //hit members to appear in fUnitIO
  Int_t cell_id=0; // cell ID in the module
  Double_t en=0; // hit energy, GeV

  Hit() {}
  Hit(Int_t i): cell_id(i) {}

};//Hit

//hits collection
class Coll : public DetectorData<Hit, std::unordered_map<Int_t, Hit>> {

  public:

  Coll() {

    //hits memory representation, the names will be a suffix to the detector name
    DATA_ADD_UNIT_ATTR( cell_id )
    DATA_ADD_UNIT_ATTR( en )

  }//Coll

};//Coll

}//HitAtID

#endif

