
//_____________________________________________________________________________
//
// Conical beampipe section
//
//_____________________________________________________________________________

//Geant
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Cons.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//local classes
#include "ConeBeam.h"
#include "GeoParser.h"
#include "ColorDecoder.h"

using namespace std;

//_____________________________________________________________________________
ConeBeam::ConeBeam(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "ConeBeam: " << fNam << G4endl;

  //center position along z, mm
  G4double zpos = 0;
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //outer radius (rmax) and wall inner radius (rmin) at lower z (0) and higher z (1), mm
  G4double rmin0 = geo->GetD(fNam, "rmin0")*mm;
  G4double rmax0 = geo->GetD(fNam, "rmax0")*mm;
  G4double rmin1 = geo->GetD(fNam, "rmin1")*mm;
  G4double rmax1 = geo->GetD(fNam, "rmax1")*mm;

  //full length along z, mm
  G4double dz = geo->GetD(fNam, "dz")*mm;

  //outer shape with vacuum
  G4Cons *shape_outer = new G4Cons(fNam, 0, rmax0, 0, rmax1, dz/2., 0, 360*deg);

  G4Material *vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume *vol_outer = new G4LogicalVolume(shape_outer, vac, fNam);
  vol_outer->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //wall of a given material
  G4String nam_wall = fNam+"_wall";
  G4Cons *shape_wall = new G4Cons(nam_wall, rmin0, rmax0, rmin1, rmax1, dz/2., 0, 360*deg);

  //material for the wall
  G4String wall_mat_name = "G4_Al";
  geo->GetOptS(nam, "wall_mat_name", wall_mat_name);

  G4Material *mat_wall = G4NistManager::Instance()->FindOrBuildMaterial(wall_mat_name);
  G4LogicalVolume *vol_wall = new G4LogicalVolume(shape_wall, mat_wall, nam_wall);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), vol_wall, nam_wall, vol_outer, false, 0);

  //wall visibility
  ColorDecoder dec("0.5:0.5:0.5:1");
  vol_wall->SetVisAttributes( dec.MakeVis(geo, fNam, "vis") );

  new G4PVPlacement(0, G4ThreeVector(0, 0, zpos), vol_outer, fNam, top, false, 0);

}//ConeBeam











