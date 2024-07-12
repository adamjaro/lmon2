
//_____________________________________________________________________________
//
// Universal quadrupole - dipole beam magnet, rectangular outer vessel,
// circular inner magnet volume
//
//_____________________________________________________________________________

//ROOT

//Geant
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleTable.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"

//local classes
#include "GeoParser.h"
#include "ColorDecoder.h"
#include "QDMagnet.h"

using namespace std;

//_____________________________________________________________________________
QDMagnet::QDMagnet(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "QDMagnet: " << fNam << G4endl;

  //magnet center along x and z, mm
  G4double xpos = 0, zpos = 0;
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //magnet length along z, mm
  G4double length = geo->GetD(fNam, "length")*mm;

  G4cout << fNam << ", xpos, zpos, length: " << xpos << " " << zpos << " " << length << G4endl;

  //inner radius, mm
  G4double dr = geo->GetD(fNam, "dr")*mm;

  //vessel size in x and y
  G4double dx = geo->GetD(fNam, "dx")*mm;
  G4double dy = geo->GetD(fNam, "dy")*mm;

  //rotation angle along y axis
  G4double theta = 0;
  geo->GetOptD(fNam, "theta", theta, GeoParser::Unit(rad));

  //beam gamma
  fGamma = geo->GetD(fNam, "gamma");

  //dipole bending field
  G4double magnet_angle = 0;
  G4bool is_dipole = geo->GetOptD(fNam, "magnet_angle", magnet_angle, GeoParser::Unit(rad));

  //quadrupole magnet
  G4double K1L = 0;
  G4bool is_quadrupole = geo->GetOptD(fNam, "K1L", K1L);

  G4cout << "is_dipole: " << is_dipole << G4endl;
  G4cout << "is_quadrupole: " << is_quadrupole << G4endl;

  //mass for beam particle
  fMass = G4ParticleTable::GetParticleTable()->FindParticle("e-")->GetPDGMass();

  //outer vessel
  G4Box *shape_vessel = new G4Box(fNam+"_vessel", dx/2, dy/2, length/2);
  G4Material *mat_vessel = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
  G4LogicalVolume *vol_vessel = new G4LogicalVolume(shape_vessel, mat_vessel, shape_vessel->GetName());

  //vessel visibility
  ColorDecoder vessel_vis("1:0:1:2"); //red:green:blue:alpha  magenta: 1:0:1
  vol_vessel->SetVisAttributes(vessel_vis.MakeVis(geo, fNam, "vessel_vis"));

  //put the vessel to the mother volume
  G4ThreeVector pos(xpos, 0, zpos);
  G4RotationMatrix rot(G4ThreeVector(0, 1, 0), theta); //CLHEP::HepRotation
  G4Transform3D transform(rot, pos); //HepGeom::Transform3D

  G4LogicalVolume *mother_vol = top;
  new G4PVPlacement(transform, vol_vessel, vol_vessel->GetName(), mother_vol, false, 0);

  //inner magnet vacuum volume
  G4Tubs *shape_inner = new G4Tubs(fNam+"_inner", 0, dr, length/2, 0, 2*M_PI*rad);
  G4Material *mat_vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume *vol_inner = new G4LogicalVolume(shape_inner, mat_vac, shape_inner->GetName());

  //visibility for inner volume
  ColorDecoder inner_vis("0:0:1:2"); //red:green:blue:alpha
  vol_inner->SetVisAttributes(inner_vis.MakeVis(geo, fNam, "inner_vis"));

  //inner vacuum volume in outer magnet vessel
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), vol_inner, vol_inner->GetName(), vol_vessel, false, 0);

  //magnetic field for inner vacuum volume
  G4MagneticField *field = nullptr;

  if( is_dipole ) field = MakeDipoleField(length, magnet_angle);
  if( is_quadrupole ) field = MakeQuadrupoleField(length, K1L);

  if(field) {

    G4FieldManager *fieldman = new G4FieldManager(field);
    fieldman->SetDetectorField(field);
    fieldman->CreateChordFinder(field);
    vol_inner->SetFieldManager(fieldman, true);

  } else {
    G4cout << "No field applied for " << fNam << G4endl;
  }

}//QDMagnet

//_____________________________________________________________________________
G4MagneticField* QDMagnet::MakeDipoleField(G4double length, G4double angle) {

  // By[T]*rho[m] = 3.33564*pc[GeV]
  // rho[m] = l[m]/(2*sin(theta[rad]/2)) (RBEND)
  // pc = mc^2*sqrt(gamma^2 - 1)
  G4double By = ((3.33564*(fMass/GeV)*sqrt(pow(fGamma,2)-1.0))/
    ((length/m)/(2.*sin(angle/2.))))*tesla;

  G4cout << "QDMagnet::MakeDipoleField for " << fNam << ", ";
  G4cout << length << " " << angle << " " << fMass << " " << fGamma << " " << By << G4endl;;

  return new G4UniformMagField( G4ThreeVector(0, By, 0) );

}//MakeDipoleField

//_____________________________________________________________________________
G4MagneticField* QDMagnet::MakeQuadrupoleField(G4double length, G4double K1L) {

  // k1[1/m^2] = 0.2998*grad[T/m]/pc[GeV], k1[1/m^2] = k1l[1/m]/l[m]
  G4double grad = (3.33564*(K1L/(1/m))/(length/m)*
    (fMass/GeV)*sqrt(pow(fGamma,2)-1.0))*(tesla/m);

  G4cout << "QDMagnet::MakeQuadrupoleField for " << fNam << ", ";
  G4cout << length << " " << K1L << " " << fMass << " " << fGamma << " " << grad << G4endl;;

  return new G4QuadrupoleMagField(grad);

}//MakeQuadrupoleField



























