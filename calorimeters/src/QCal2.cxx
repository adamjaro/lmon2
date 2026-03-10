
//_____________________________________________________________________________
//
// Cherenkov fiber calorimeter with fibers at an angle to its axis
//
//_____________________________________________________________________________

//C++

//ROOT
#include "TTree.h"
#include "TMath.h"

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
#include "ColorDecoder.h"
#include "LogicVolFinder.h"
#include "QCal2.h"

//_____________________________________________________________________________
QCal2::QCal2(const G4String& nam, GeoParser *geo, G4LogicalVolume *top): Detector(),
   G4VSensitiveDetector(nam), fNam(nam) {

  G4cout << "QCal2: " << fNam << G4endl;

  //module size
  G4double modx = geo->GetD(fNam, "modx")*mm; // full size in x, mm
  G4double mody = geo->GetD(fNam, "mody")*mm; // full size in y, mm
  G4double modz = geo->GetD(fNam, "modz")*mm; // full size in z, mm

  //calorimeter module
  //G4Box *mod_shape = new G4Box(fNam+"_mod", modx/2., mody/2., modz/2.);
  G4Box *mod_shape = new G4Box(fNam, modx/2., mody/2., modz/2.);
  G4Material *mod_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *mod_vol = new G4LogicalVolume(mod_shape, mod_mat, mod_shape->GetName());

  //module visibility
  ColorDecoder mod_vis("0:0:1:2");
  mod_vol->SetVisAttributes(mod_vis.MakeVis(geo, fNam, "mod_vis"));

  //cell volume
  G4LogicalVolume *cell_vol = MakeCell(geo);

  //cells in module
  G4RotationMatrix *cell_rot = new G4RotationMatrix(); // is HepRotation
  cell_rot->rotateX((90+geo->GetD(fNam, "cell_phi"))*deg);

  G4int cell_cnt = 0; // cell count in module

  //initial position for the first cell along z
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm;
  G4double cell_z = geo->GetD(fNam, "cell_z")*mm;
  G4double cell_phi = geo->GetD(fNam, "cell_phi")*deg; // the 'cell_phi' value is in rad after geo-> call is *deg
  G4double cell_posz_init = 0.5*(cell_xy*TMath::Cos(cell_phi) + cell_z*TMath::Sin(cell_phi));

  //number of cells in module along z and x
  G4int nz = geo->GetI(fNam, "nz"); // num cells in z

  G4double cell_posz = 0.5*modz - cell_posz_init;

  //z-loop
  for(G4int iz=0; iz<nz; iz++) {

    new G4PVPlacement(cell_rot, G4ThreeVector(0, 0, cell_posz), cell_vol, cell_vol->GetName(), mod_vol, false, cell_cnt++);
    cell_posz -= cell_xy/TMath::Cos(cell_phi);

  }//z-loop

/*

  new G4PVPlacement(cell_rot, G4ThreeVector(0, 0, cell_posz), cell_vol, cell_vol->GetName(), mod_vol, false, cell_cnt++);
*/


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

//_____________________________________________________________________________
G4LogicalVolume* QCal2::MakeCell(GeoParser *geo) {

  //cell (absorber) size
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm;
  G4double cell_z = geo->GetD(fNam, "cell_z")*mm;






  //local materials for the cell
  if( !G4Material::GetMaterial("QCal2_Copper", false) ) {
    new G4Material("QCal2_Copper", 29., 63.55*g/mole, 8.960*g/cm3); //z, a, density
  }

  //name for cell material
  G4String cell_mat_name = "G4_W";
  geo->GetOptS(fNam, "cell_mat_name", cell_mat_name);

  //select the cell material
  G4Material *cell_mat=0x0;
  if( cell_mat_name.find("QCal2_") != std::string::npos ) {
    cell_mat = G4Material::GetMaterial(cell_mat_name); // local material
  } else {
    cell_mat = G4NistManager::Instance()->FindOrBuildMaterial(cell_mat_name); // nist material
  }

  G4cout << "  " << fNam << ", cell material: " << cell_mat->GetName() << G4endl;

  //cell volume
  G4Box *cell_shape = new G4Box(fNam+"_cell", cell_xy/2., cell_xy/2., cell_z/2.);
  G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, cell_mat, cell_shape->GetName());

  //cell visibility
  ColorDecoder cell_vis("1:0:0:2");
  cell_vol->SetVisAttributes(cell_vis.MakeVis(geo, fNam, "cell_vis"));



  return cell_vol;

}//MakeCell





















