
//_____________________________________________________________________________
//
// Cherenkov fiber calorimeter with fibers at an angle to its axis
//
//_____________________________________________________________________________

//C++

//ROOT
#include "TTree.h"

//Geant
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4VisAttributes.hh"

//local classes
#include "GeoParser.h"
#include "LogicVolFinder.h"
#include "QCal2.h"

//_____________________________________________________________________________
QCal2::QCal2(const G4String& nam, GeoParser *geo, G4LogicalVolume *top): Detector(),
   G4VSensitiveDetector(nam), fNam(nam) {

  G4cout << "QCal2: " << fNam << G4endl;

  //module size
  G4double modx=10*mm, mody=10*mm, modz=10*mm;

  //calorimeter module
  //G4Box *mod_shape = new G4Box(fNam+"_mod", modx/2., mody/2., modz/2.);
  G4Box *mod_shape = new G4Box(fNam, modx/2., mody/2., modz/2.);
  G4Material *mod_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *mod_vol = new G4LogicalVolume(mod_shape, mod_mat, mod_shape->GetName());




  //module center position
  G4double xpos=0, ypos=0, zpos=0;

  //module in its mother volume
  G4LogicalVolume *mother_vol = LogicVolFinder::GetMotherVolume("place_into", top, geo, fNam);
  //case of two levels of volume to place, volume 'place_into2' inside 'place_into1'
  if( geo->HasParameter(fNam, "place_into2") ) {
    mother_vol = LogicVolFinder::GetMotherVolume2("place_into1", "place_into2", top, geo, fNam);
  }

  new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos), mod_vol, mod_vol->GetName(), mother_vol, false, 0);

}//QCal2























