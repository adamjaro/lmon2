
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
#include "CaloBuilder.h"

//local detectors and components
#include "CaloBPC.h"
#include "HcalA262.h"
#include "UcalA290.h"
#include "WScFiZXv3.h"
#include "CalPWO.h"
#include "QCal1.h"
#include "FiberXYCal.h"

//macros
#define ADD_DETECTOR(det) (fDets.insert( make_pair(#det, &CaloBuilder::MakeDet<det>) ))

//_____________________________________________________________________________
CaloBuilder::CaloBuilder(G4LogicalVolume *top, GeoParser *geo): BuilderBase(),
  fTop(top), fGeo(geo) {

  G4cout << "CaloBuilder::CaloBuilder" << G4endl;

  //individual detectors as defined here
  ADD_DETECTOR( CaloBPC );
  ADD_DETECTOR( HcalA262 );
  ADD_DETECTOR( UcalA290 );
  ADD_DETECTOR( WScFiZXv3 );
  ADD_DETECTOR( CalPWO );
  ADD_DETECTOR( QCal1 );
  ADD_DETECTOR( FiberXYCal );

}//CaloBuilder

//_____________________________________________________________________________
Detector* CaloBuilder::FindAndLoad(G4String type, G4String name) {

  //G4cout << "CaloBuilder::FindAndLoad, " << type << " " << name << G4endl;

  //factory detector construction
  std::map<G4String, MakeDetPtr>::iterator idet = fDets.find(type);
  if( idet == fDets.end() ) return 0x0;

  //detector found, make its instance
  Detector *det = (this->*(*idet).second)(name, fGeo, fTop);

  return det;

}//FindAndLoad


















