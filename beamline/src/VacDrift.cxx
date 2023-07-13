
//_____________________________________________________________________________
//
// Vacuum drift section
//
//_____________________________________________________________________________

//Boost
#include <boost/tokenizer.hpp>

//ROOT
#include "TMath.h"

//Geant
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4GenericTrap.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//local classes
#include "VacDrift.h"
#include "GeoParser.h"
#include "ColorDecoder.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
VacDrift::VacDrift(G4String nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "VacDrift: " << fNam << G4endl;

  // full size in y, mm
  G4double ysiz = geo->GetD(fNam, "ysiz")*mm;

  //at lower z
  G4double z0T = geo->GetD(fNam, "z0T")*mm;
  G4double x0T = geo->GetD(fNam, "x0T")*mm;
  G4double z0B = z0T;
  geo->GetOptD(fNam, "z0B", z0B, GeoParser::Unit(mm));
  G4double x0B = geo->GetD(fNam, "x0B")*mm;

  //at larger z
  G4double z1T = geo->GetD(fNam, "z1T")*mm;
  G4double x1T = geo->GetD(fNam, "x1T")*mm;
  G4double z1B = z1T;
  geo->GetOptD(fNam, "z1B", z1B, GeoParser::Unit(mm));
  G4double x1B = geo->GetD(fNam, "x1B")*mm;

  //vessel shape
  G4GenericTrap *vessel_shape = MakeGT(z0T, x0T, z0B, x0B, z1T, x1T, z1B, x1B, ysiz, fNam+"_vessel_shape");

  //material for the vessel
  G4String mat_name = "G4_Galactic";
  geo->GetOptS(nam, "mat_name", mat_name);
  G4cout << "  " << fNam << ", mat_name: " << mat_name << G4endl;

  //vessel logical volume
  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial(mat_name);
  G4LogicalVolume *vol = new G4LogicalVolume(vessel_shape, mat, fNam);

  //vessel visibility
  ColorDecoder vessel_vis("0:0:1:2"); //red:green:blue:alpha
  vol->SetVisAttributes(vessel_vis.MakeVis(geo, fNam, "vis"));

  //vertical center position (along y), mm
  G4double ypos=0;
  geo->GetOptD(fNam, "ypos", ypos, GeoParser::Unit(mm));

  //placement in top
  G4ThreeVector pos(0, ypos, 0);
  G4RotationMatrix rot(G4ThreeVector(1, 0, 0), TMath::Pi()/2); //CLHEP::HepRotation
  G4Transform3D transform(rot, pos); //HepGeom::Transform3D

  new G4PVPlacement(transform, vol, fNam, top, false, 0);

}//VacDrift

//_____________________________________________________________________________
G4GenericTrap *VacDrift::MakeGT(
    G4double z0T, G4double x0T, G4double z0B, G4double x0B, G4double z1T, G4double x1T, G4double z1B, G4double x1B,
    G4double ysiz, G4String nam) {

  //G4GenericTrap or TGeoArb8
  //generic trapezoid native coordinates: 4 xy points plane at -dz, 4 xy points plane at +dz, both clockwise
  //rotation by +pi/2 about x from generic trapezoid coordinates to detector frame: y -> z,  z -> y

  //vertices for the trapezoid
  vector<G4TwoVector> ver(8);

  ver[0].set(x0B, z0B); // point #1

  ver[1].set(x1B, z1B); // point #2

  ver[2].set(x1T, z1T); // point #3

  ver[3].set(x0T, z0T); // point #4

  //plane at lower y
  for(int i=4; i<8; i++) {
    ver[i].set(ver[i-4].x(), ver[i-4].y());
  }

  return new G4GenericTrap(nam, ysiz/2, ver);

}//MakeGT




