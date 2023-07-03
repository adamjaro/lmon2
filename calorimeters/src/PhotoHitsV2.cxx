
//_____________________________________________________________________________
//
// Hits for optical detector, PMT or SiPM, version V2
// derived from the DetectorData template
//_____________________________________________________________________________

//C++
#include <iostream>

//ROOT
#include "TTree.h"

//Geant
//#include "G4ios.hh"

//local classes
#include "PhotoHitsV2.h"

using namespace std;

//_____________________________________________________________________________
PhotoHitsV2::Coll::Coll() {

  //G4cout << "PhotoHitsV2::Coll::Coll" << G4endl;

  //hits memory representation, the names will be a suffix to the detector name
  DATA_ADD_UNIT_ATTR( pos_x )
  DATA_ADD_UNIT_ATTR( pos_y )
  DATA_ADD_UNIT_ATTR( pos_z )
  DATA_ADD_UNIT_ATTR( time )
  DATA_ADD_UNIT_ATTR( pmt_x )
  DATA_ADD_UNIT_ATTR( pmt_y )
  DATA_ADD_UNIT_ATTR( pmt_z )
  DATA_ADD_UNIT_ATTR( cell_id )
  DATA_ADD_UNIT_ATTR( prim_id )

}//Coll













