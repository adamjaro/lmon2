
//_____________________________________________________________________________
//
// Cylindrical construction segment
//
//_____________________________________________________________________________

//ROOT
#include "TMath.h"

//Geant
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//local classes
#include "CylSegment.h"
#include "GeoParser.h"
#include "ColorDecoder.h"

using namespace std;

//_____________________________________________________________________________
CylSegment::CylSegment(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "CylSegment: " << fNam << G4endl;

  //center position along z, mm
  G4double zpos = 0;
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //full length along z, mm
  G4double dz = geo->GetD(fNam, "dz")*mm;

  //inner and outer radius
  G4double rmin = geo->GetD(fNam, "rmin")*mm;
  G4double rmax = geo->GetD(fNam, "rmax")*mm;

  G4Tubs *shape = new G4Tubs(fNam, rmin, rmax, dz/2., 0., 360.*deg);

  //logical volume
  G4String mat_name = "G4_Al";
  geo->GetOptS(nam, "mat_name", mat_name);
  G4cout << "  " << fNam << ", mat_name: " << mat_name << G4endl;
  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial(mat_name);
  G4LogicalVolume *vol = new G4LogicalVolume(shape, mat, fNam);

  //visibility
  ColorDecoder dec("0:0:1:2"); //red:green:blue:alpha
  vol->SetVisAttributes( dec.MakeVis(geo, fNam, "vis") );

  //segment in the top volume
  new G4PVPlacement(0, G4ThreeVector(0, 0, zpos), vol, fNam, top, false, 0);

}//CylSegment

























