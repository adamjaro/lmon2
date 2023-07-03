
//_____________________________________________________________________________
//
// Beam drift section between B2 and Q3 magnets
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
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//local classes
#include "BeamDrift.h"
#include "GeoParser.h"
#include "ColorDecoder.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
BeamDrift::BeamDrift(G4String nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "BeamDrift: " << fNam << G4endl;

  // full inner size in y, mm
  G4double ysiz = geo->GetD(fNam, "ysiz")*mm;

  //wall thickness, for vertical size
  G4double delta = geo->GetD(fNam, "delta")*mm;

  //at lower z
  //outer
  G4double z0TO = geo->GetD(fNam, "z0TO")*mm;
  G4double x0TO = geo->GetD(fNam, "x0TO")*mm;
  G4double z0BO = z0TO;
  geo->GetOptD(fNam, "z0BO", z0BO, GeoParser::Unit(mm));
  G4double x0BO = geo->GetD(fNam, "x0BO")*mm;
  //inner
  G4double z0TI = z0TO;
  geo->GetOptD(fNam, "z0TI", z0TI, GeoParser::Unit(mm));
  G4double x0TI = geo->GetD(fNam, "x0TI")*mm;
  G4double z0BI = z0TO;
  geo->GetOptD(fNam, "z0BI", z0BI, GeoParser::Unit(mm));
  G4double x0BI = geo->GetD(fNam, "x0BI")*mm;

  //at larger z
  //outer
  G4double z1TO = geo->GetD(fNam, "z1TO")*mm;
  G4double x1TO = geo->GetD(fNam, "x1TO")*mm;
  G4double z1BO = z1TO;
  geo->GetOptD(fNam, "z1BO", z1BO, GeoParser::Unit(mm));
  G4double x1BO = geo->GetD(fNam, "x1BO")*mm;
  //inner
  G4double z1TI = z1TO;
  geo->GetOptD(fNam, "z1TI", z1TI, GeoParser::Unit(mm));
  G4double x1TI = geo->GetD(fNam, "x1TI")*mm;
  G4double z1BI = z1TO;
  geo->GetOptD(fNam, "z1BI", z1BI, GeoParser::Unit(mm));
  G4double x1BI = geo->GetD(fNam, "x1BI")*mm;

  //outer and inner vessel
  G4GenericTrap *vessel_outer = MakeGT(z0TO, x0TO, z0BO, x0BO, z1TO, x1TO, z1BO, x1BO, ysiz+2*delta, fNam+"_vessel_outer");
  G4GenericTrap *vessel_inner = MakeGT(z0TI, x0TI, z0BI, x0BI, z1TI, x1TI, z1BI, x1BI, ysiz, fNam+"_vessel_inner");

  //vessel shape
  G4SubtractionSolid *vessel_shape = new G4SubtractionSolid(fNam+"_vessel_shape", vessel_outer, vessel_inner);

  //vessel logical volume
  G4Material *vessel_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4LogicalVolume *vessel_vol = new G4LogicalVolume(vessel_shape, vessel_mat, fNam+"_vessel_shape");

  //vessel visibility
  ColorDecoder vessel_dec("0:0:1:2"); //red:green:blue:alpha
  vessel_vol->SetVisAttributes(vessel_dec.MakeVis(geo, fNam, "vis"));

  //outer logical volume
  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume *vol = new G4LogicalVolume(vessel_outer, mat, fNam);

  //visibility for inner and outer volume
  ColorDecoder inout_dec("0:0:0:3");

  vol->SetVisAttributes( inout_dec.MakeVis(geo, fNam, "vis_inout") );

  //vessel volume in outer shape
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), vessel_vol, fNam+"_vessel_shape", vol, false, 0);

  //inner volume
  G4LogicalVolume *inner_vol = new G4LogicalVolume(vessel_inner, mat, fNam+"_vessel_inner");
  inner_vol->SetVisAttributes( inout_dec.MakeVis(geo, fNam, "vis_inout") );

  //inner volume in outer shape
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), inner_vol, fNam+"_vessel_inner", vol, false, 0);

  //placement in top
  G4ThreeVector pos(0, 0, 0);
  G4RotationMatrix rot(G4ThreeVector(1, 0, 0), TMath::Pi()/2); //CLHEP::HepRotation
  G4Transform3D transform(rot, pos); //HepGeom::Transform3D

  new G4PVPlacement(transform, vol, fNam, top, false, 0);

}//BeamDrift

//_____________________________________________________________________________
G4GenericTrap *BeamDrift::MakeGT(
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















