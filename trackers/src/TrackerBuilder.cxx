
//_____________________________________________________________________________
//
// construction of calorimeters
//
//_____________________________________________________________________________

//C++
#include <vector>

//ROOT
#include "Rtypes.h"

//Geant
#include "G4LogicalVolume.hh"

//local classes
#include "GeoParser.h"
#include "TrackerBuilder.h"

//local detectors and components
#include "TrkPlaneBasic.h"
#include "Timepix4v1.h"

//macros
#define ADD_DETECTOR(det) (fDets.insert( make_pair(#det, &TrackerBuilder::MakeDet<det>) ))

//_____________________________________________________________________________
TrackerBuilder::TrackerBuilder(G4LogicalVolume *top, GeoParser *geo): BuilderBase(),
  fTop(top), fGeo(geo) {

  G4cout << "TrackerBuilder::TrackerBuilder" << G4endl;

  ADD_DETECTOR( TrkPlaneBasic );
  ADD_DETECTOR( Timepix4v1 );

}//TrackerBuilder

//_____________________________________________________________________________
Detector* TrackerBuilder::FindAndLoad(G4String type, G4String name) {

  //G4cout << "TrackerBuilder::FindAndLoad, " << type << " " << name << G4endl;

  //factory detector construction
  std::map<G4String, MakeDetPtr>::iterator idet = fDets.find(type);
  if( idet == fDets.end() ) return 0x0;

  //detector found, make its instance
  Detector *det = (this->*(*idet).second)(name, fGeo, fTop);

  return det;

}//FindAndLoad


















