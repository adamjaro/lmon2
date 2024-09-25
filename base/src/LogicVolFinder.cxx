
//_____________________________________________________________________________
//
// Helper to locate a given volume in hierarchy of volumes
//
//_____________________________________________________________________________

//Geant
#include "G4LogicalVolume.hh"

//local classes
#include "GeoParser.h"
#include "LogicVolFinder.h"

//_____________________________________________________________________________
G4LogicalVolume* LogicVolFinder::GetMotherVolume(G4String place_into, G4LogicalVolume *top, GeoParser *geo, G4String det_name) {

  //optional mother volume in geometry

  G4LogicalVolume *mvol = top;
  G4String mother_nam;
  if( geo->GetOptS(det_name, place_into, mother_nam) ) {
    mvol = GetMotherVolume(mother_nam, top);
  }

  return mvol;

}//GetMotherVolume

//_____________________________________________________________________________
G4LogicalVolume* LogicVolFinder::GetMotherVolume2(G4String p1, G4String p2, G4LogicalVolume *top, GeoParser *geo, G4String det_name) {

  //optional recursive search for volume named 'p2' in 'p1' with p1 itself in top

  G4LogicalVolume *mvol = top;

  if( geo->HasParameter(det_name, p2) ) {
    mvol = GetMotherVolume2( geo->GetS(det_name, p1), geo->GetS(det_name, p2), top );
  }

  return mvol;

}//GetMotherVolume2

//_____________________________________________________________________________
G4LogicalVolume* LogicVolFinder::GetMotherVolume(G4String mother_nam, G4LogicalVolume *top) {

  //get volume named 'mother_nam' in volume provided by 'top'

  for(size_t i=0; i<top->GetNoDaughters(); i++) {

    G4LogicalVolume *dv = top->GetDaughter(i)->GetLogicalVolume();

    if( dv->GetName() == mother_nam ) {
      return dv;
    }
  }

  return nullptr;

}//GetMotherVolume

//_____________________________________________________________________________
G4LogicalVolume* LogicVolFinder::GetMotherVolume2(G4String m1, G4String m2, G4LogicalVolume *top) {

  //recursive search for volume m2 in m1 with m1 itself in top

  G4LogicalVolume *vol1 = GetMotherVolume(m1, top);

  if(!vol1) return nullptr;

  return GetMotherVolume(m2, vol1);

}//GetMotherVolume2

