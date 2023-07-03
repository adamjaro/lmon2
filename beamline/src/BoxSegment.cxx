
//_____________________________________________________________________________
//
// Construction box segment
//
//_____________________________________________________________________________

//Geant
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//local classes
#include "BoxSegment.h"
#include "GeoParser.h"
#include "ColorDecoder.h"

using namespace std;

//_____________________________________________________________________________
BoxSegment::BoxSegment(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "BoxSegment: " << fNam << G4endl;

  //full size in x, y and z, mm
  G4double dx = geo->GetD(fNam, "dx")*mm;
  G4double dy = geo->GetD(fNam, "dy")*mm;
  G4double dz = geo->GetD(fNam, "dz")*mm;
  G4Box *shape = new G4Box(fNam, dx/2., dy/2., dz/2.);

  //logical volume
  G4String mat_name = "G4_Galactic";
  geo->GetOptS(nam, "mat_name", mat_name);
  G4cout << "  " << fNam << ", mat_name: " << mat_name << G4endl;
  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial(mat_name);
  G4LogicalVolume *vol = new G4LogicalVolume(shape, mat, fNam);

  //visibility
  ColorDecoder vis("0:0:1:2");
  vol->SetVisAttributes( vis.MakeVis(geo, fNam, "vis") );

  //center position, mm
  G4double xpos=0*mm, ypos=0*mm, zpos=0*mm;
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "ypos", ypos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  G4cout << "  " << fNam << ", xpos: " << xpos << G4endl;
  G4cout << "  " << fNam << ", ypos: " << ypos << G4endl;
  G4cout << "  " << fNam << ", zpos: " << zpos << G4endl;

  //polar angle along the y axis
  G4double theta = 0;
  geo->GetOptD(fNam, "theta", theta, GeoParser::Unit(rad));

  //select the rotation axis
  G4ThreeVector rot_axis(0, 1, 0); // y
  G4bool rotate_x = false;
  geo->GetOptB(fNam, "rotate_x", rotate_x);
  if(rotate_x) {
    rot_axis = G4ThreeVector(1, 0, 0); // x
  }

  G4RotationMatrix rot(rot_axis, theta); //CLHEP::HepRotation
  G4cout << "  " << fNam << ", theta: " << theta << G4endl;

  //placement in mother volume
  G4ThreeVector pos(xpos, ypos, zpos);
  G4Transform3D transform(rot, pos); //HepGeom::Transform3D

  //get the mother volume
  G4LogicalVolume *mother_vol = top;
  G4String mother_nam;
  if( geo->GetOptS(fNam, "place_into", mother_nam) ) {
    mother_vol = GetMotherVolume(mother_nam, top);
  }

  //make the placement
  new G4PVPlacement(transform, vol, fNam, mother_vol, false, 0);

}//BoxSegment

//_____________________________________________________________________________
G4LogicalVolume* BoxSegment::GetMotherVolume(G4String mother_nam, G4LogicalVolume *top) {

  for(size_t i=0; i<top->GetNoDaughters(); i++) {

    G4LogicalVolume *dv = top->GetDaughter(i)->GetLogicalVolume();

    if( dv->GetName() == mother_nam ) {
      return dv;
    }
  }

  return 0x0;

}//GetMotherVolume












