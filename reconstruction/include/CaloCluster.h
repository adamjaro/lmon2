
#ifndef CaloCluster_h
#define CaloCluster_h

//general calorimeter cluster

#include "DetectorData.h"

namespace CaloCluster {

//cluster representation
struct Cls {

  Double_t x=0; // position in x, mm
  Double_t y=0; // position in y, mm
  Double_t en=0; // cluster energy, GeV
  Int_t nhit=0; // number of contributing hits

};//Cls

//cluster collection
class Coll: public DetectorData<Cls> {

  public:

    Coll() {

      //clusters memory representation
      DATA_ADD_UNIT_ATTR( x )
      DATA_ADD_UNIT_ATTR( y )
      DATA_ADD_UNIT_ATTR( en )
      DATA_ADD_UNIT_ATTR( nhit )

    }

    //access to created clusters during reconstruction
    std::vector<Cls>& GetStore() { return fStorage; }

};//Coll

}//CaloCluster

#endif










