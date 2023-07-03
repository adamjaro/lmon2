
#ifndef MCParticles_h
#define MCParticles_h

#include "DetectorData.h"

namespace MCParticles {

//particle representation
struct Part {

  Int_t pdg; // particle pdg
  Int_t itrk; // track ID for the particle
  Double_t en; // particle energy, GeV
  Double_t theta; // particle polar angle, rad
  Double_t phi; // particle azimuthal angle, rad
  Double_t vx; // particle vertex position in x, mm
  Double_t vy; // particle vertex position in y, mm
  Double_t vz; // particle vertex position in z, mm

};//Part

//particles collection
class Coll : public DetectorData<Part> {

  public:

    Coll() {

      //particles memory representation
      DATA_ADD_UNIT_ATTR( pdg )
      DATA_ADD_UNIT_ATTR( itrk )
      DATA_ADD_UNIT_ATTR( en )
      DATA_ADD_UNIT_ATTR( theta )
      DATA_ADD_UNIT_ATTR( phi )
      DATA_ADD_UNIT_ATTR( vx )
      DATA_ADD_UNIT_ATTR( vy )
      DATA_ADD_UNIT_ATTR( vz )
    }

};//Coll

}//MCParticles

#endif














