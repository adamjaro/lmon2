
#ifndef TaggerTracks_h
#define TaggerTracks_h

#include "DetectorData.h"

namespace TaggerTracks {

  //track representation
  struct Track {

    Double_t x=0; // track position in x, mm
    Double_t y=0; // track position in y, mm
    Double_t theta_x=0; // track angle along x, rad
    Double_t theta_y=0; // track angle along y, rad
    Double_t chi2_xy=0; // track chi^2 in xy plane

  };//Track

  //tracks collection
  class Coll : public DetectorData<Track> {

    public:

      Coll() {

  //data members of Track to be written in the output from Coll of Track objects
  DATA_ADD_UNIT_ATTR( x )
  DATA_ADD_UNIT_ATTR( y )
  DATA_ADD_UNIT_ATTR( theta_x )
  DATA_ADD_UNIT_ATTR( theta_y )
  DATA_ADD_UNIT_ATTR( chi2_xy )

      }

      //std::vector<Track>& GetStore() { return fStorage; }

    };//Coll

}//TaggerTracks

#endif

