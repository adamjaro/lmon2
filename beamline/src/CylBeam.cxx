
//_____________________________________________________________________________
//
// Cylindrical beampipe section
//
//_____________________________________________________________________________

//ROOT
#include "TMath.h"

//Geant
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4CutTubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//local classes
#include "CylBeam.h"
#include "GeoParser.h"
#include "ColorDecoder.h"

using namespace std;

//_____________________________________________________________________________
CylBeam::CylBeam(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "CylBeam: " << fNam << G4endl;

  //center position along z, mm
  G4double zpos = 0;
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //full length along z, mm
  G4double dz = geo->GetD(fNam, "dz")*mm;

  //inner and outer radius
  G4double rmin = geo->GetD(fNam, "rmin")*mm;
  G4double rmax = geo->GetD(fNam, "rmax")*mm;

  //angle along y at lower z (0) and higher z (1)
  G4double theta0 = 0., theta1 = 0.;
  geo->GetOptD(fNam, "theta0", theta0, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "theta1", theta1, GeoParser::Unit(mm));

  //outer shape with vacuum
  G4CutTubs *shape_outer = new G4CutTubs(fNam, 0, rmax, dz/2., 0, 360*deg,
    G4ThreeVector(-TMath::Sin(theta0), 0, -TMath::Cos(theta0)), G4ThreeVector(TMath::Sin(theta1), 0, TMath::Cos(theta1)));

  G4Material *vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume *vol_outer = new G4LogicalVolume(shape_outer, vac, fNam);
  vol_outer->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //wall of a given material
  G4CutTubs *shape_wall = new G4CutTubs(fNam+"_wall", rmin, rmax, dz/2., 0, 360*deg,
    G4ThreeVector(-TMath::Sin(theta0), 0, -TMath::Cos(theta0)), G4ThreeVector(TMath::Sin(theta1), 0, TMath::Cos(theta1)));

  //material for the wall
  G4String wall_mat_name = "G4_Al";
  geo->GetOptS(nam, "wall_mat_name", wall_mat_name);

  //wall in outer volume
  G4Material *mat_wall = G4NistManager::Instance()->FindOrBuildMaterial(wall_mat_name);
  G4LogicalVolume *vol_wall = new G4LogicalVolume(shape_wall, mat_wall, fNam+"_wall");
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), vol_wall, fNam+"_wall", vol_outer, false, 0);

  //wall visibility
  ColorDecoder dec("0.5:0.5:0.5:1");
  vol_wall->SetVisAttributes( dec.MakeVis(geo, fNam, "vis") );

  //outer volume in the top
  new G4PVPlacement(0, G4ThreeVector(0, 0, zpos), vol_outer, fNam, top, false, 0);

}//CylBeam

















