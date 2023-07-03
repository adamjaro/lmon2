
#ifndef PhotoHitsV2_h
#define PhotoHitsV2_h

#include "DetectorData.h"

namespace PhotoHitsV2 {

//hit representation
struct Hit {

  //hit members to appear in fUnitIO
  Double_t pos_x; // hit position in x, mm
  Double_t pos_y; // hit position in y, mm
  Double_t pos_z; // hit position in z, mm
  Double_t time; // time of the hit, ns
  Double_t pmt_x; // global PMT position in x, mm
  Double_t pmt_y; // global PMT position in y, mm
  Double_t pmt_z; // global PMT position in z, mm
  Int_t cell_id; // cell index with the pmt in the module
  Int_t prim_id; // ID of primary particle associated with the hit

};//Hit

//hits collection
class Coll : public DetectorData<Hit> {

  public:

    Coll();

};//Coll

}//PhotoHitsV2

#endif

