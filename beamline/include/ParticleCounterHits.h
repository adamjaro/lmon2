
#ifndef ParticleCounterHits_h
#define ParticleCounterHits_h

// hits for ParticleCounter

#include "DetectorData.h"

namespace ParticleCounterHits {

//hit representation
struct Hit {

  Int_t pdg=0; // particle pdg
  Double_t en=0; // hit energy, GeV
  Double_t x=0; // hit position in x, mm
  Double_t y=0; // hit position in y, mm
  Double_t z=0; // hit position in z, mm
  Double_t xdir=0; // track direction in x
  Double_t ydir=0; // track direction in y
  Double_t zdir=0; // track direction in z
  Int_t parentID=0; // parent ID for track in hit
  Int_t itrk=0; // track index
  Bool_t is_prim=0; // hit by primary particle

};//Hit

//hits collection
class Coll : public DetectorData<Hit> {

  public:

    Coll() {

      //hits memory representation
      DATA_ADD_UNIT_ATTR( pdg )
      DATA_ADD_UNIT_ATTR( en )
      DATA_ADD_UNIT_ATTR( x )
      DATA_ADD_UNIT_ATTR( y )
      DATA_ADD_UNIT_ATTR( z )
      DATA_ADD_UNIT_ATTR( xdir )
      DATA_ADD_UNIT_ATTR( ydir )
      DATA_ADD_UNIT_ATTR( zdir )
      DATA_ADD_UNIT_ATTR( parentID )
      DATA_ADD_UNIT_ATTR( itrk )
      DATA_ADD_UNIT_ATTR( is_prim )

    }
};//Coll

}//ParticleCounterHits

#endif
















