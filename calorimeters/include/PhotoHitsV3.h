
#ifndef PhotoHitsV3_h
#define PhotoHitsV3_h

#include "Rtypes.h"

#include "DetectorData.h"

namespace PhotoHitsV3 {

//hit representation
struct Hit {

  //hit members to appear in fUnitIO
  Double_t time; // time of the hit, ns
  Double_t phot_en; // energy of optical photons, eV
  Int_t cell_id; // cell index

};//Hit

//hits collection
class Coll : public DetectorData<Hit> {

  public:

  Coll() {

    //hits memory representation, the names will be a suffix to the detector name
    DATA_ADD_UNIT_ATTR( time )
    DATA_ADD_UNIT_ATTR( phot_en )
    DATA_ADD_UNIT_ATTR( cell_id )

  }//Coll

};//Coll

}//PhotoHitsV3

#endif

