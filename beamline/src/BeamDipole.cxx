
//_____________________________________________________________________________
//
// beamline dipole magnet
//
//_____________________________________________________________________________

//C++
#include <typeinfo>

//Boost
#include <boost/tokenizer.hpp>

//Geant headers
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Cons.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4VIntegrationDriver.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4IntegrationDriver.hh"
#include "G4Mag_UsualEqRhs.hh"

//local headers
#include "BeamDipole.h"
#include "GeoParser.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
BeamDipole::BeamDipole(G4String nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "BeamDipole: " << fNam << G4endl;

  //center position along x and z, mm
  G4double zpos = 0, xpos = 0;
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));

  //polar angle along y axis
  G4double theta = 0;
  geo->GetOptD(fNam, "theta", theta, GeoParser::Unit(rad));

  //total length in z, mm
  G4double length = geo->GetD(fNam, "length")*mm;

  //conical inner core, entrance radius (r1) and exit radius (r2)
  G4double r1 = 10*mm, r2 = 10*mm;
  geo->GetOptD(fNam, "r1", r1, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "r2", r2, GeoParser::Unit(mm));

  G4String nam_inner = fNam+"_inner";
  G4Cons *shape_inner = new G4Cons(nam_inner, 0, r2, 0, r1, length/2, 0, 360*deg);

  G4Material *mat_inner = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume *vol_inner = new G4LogicalVolume(shape_inner, mat_inner, nam_inner);
  vol_inner->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //magnetic field inside the inner core
  G4double dipole_field = geo->GetD(fNam, "field") * tesla; // magnet field
  G4UniformMagField *field = new G4UniformMagField(G4ThreeVector(0, dipole_field, 0));
  G4FieldManager *fman = new G4FieldManager();
  fman->SetDetectorField(field);

  //alternative stepper
  //G4ClassicalRK4 *stepper = new G4ClassicalRK4(new G4Mag_UsualEqRhs(field));
  //G4VIntegrationDriver *driver = new G4IntegrationDriver<G4ClassicalRK4>(5e-05*mm, stepper);
  //G4ChordFinder *finder = new G4ChordFinder(driver);
  //fman->SetChordFinder(finder);
  //fman->SetMinimumEpsilonStep(5e-4);
  //fman->SetMaximumEpsilonStep(0.01);

  fman->CreateChordFinder(field);
  //PrintField(fman);

  vol_inner->SetFieldManager(fman, true);

  //cylindrical outer shape
  G4double r3 = r2*2;
  geo->GetOptD(fNam, "r3", r3, GeoParser::Unit(mm)); // vessel outer radius
  G4Tubs *shape_outer = new G4Tubs(fNam+"_outer", 0., r3, length/2-1e-4*meter, 0., 360.*deg);

  //main magnet volume
  G4LogicalVolume *main_vol = new G4LogicalVolume(shape_outer, mat_inner, fNam+"_main");
  main_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //magnet vessel around the inner magnetic core
  G4SubtractionSolid *shape_vessel = new G4SubtractionSolid(fNam, shape_outer, shape_inner);

  G4Material *mat_outer = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
  G4LogicalVolume *vol_vessel = new G4LogicalVolume(shape_vessel, mat_outer, fNam);

  //vessel visibility
  vol_vessel->SetVisAttributes(ColorDecoder(geo));

  //put the magnet vessel to the main magnet volume
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), vol_vessel, fNam, main_vol, false, 0);

  //put the inner core to the main magnet volume
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), vol_inner, nam_inner, main_vol, false, 0);

  //put main magnet volume to the top
  G4RotationMatrix main_rot(G4ThreeVector(0, 1, 0), theta); //CLHEP::HepRotation
  G4ThreeVector main_pos(xpos, 0, zpos);
  G4Transform3D main_trans(main_rot, main_pos); //HepGeom::Transform3D
  new G4PVPlacement(main_trans, main_vol, fNam, top, false, 0);

}//BeamDipole

//_____________________________________________________________________________
void BeamDipole::PrintField(G4FieldManager *fman) {

  G4ChordFinder *finder = fman->GetChordFinder();

  G4cout << "BeamDipole:" << G4endl;
  G4cout << "eps_min: " << fman->GetMinimumEpsilonStep()/mm << G4endl;
  G4cout << "eps_max: " << fman->GetMaximumEpsilonStep()/mm << G4endl;
  //G4cout << "min_chord_step: " << fman-> << G4endl;
  G4cout << "delta_chord: " << finder->GetDeltaChord()/mm << G4endl;
  G4cout << "delta_intersection: " << fman->GetDeltaIntersection()/mm << G4endl;
  G4cout << "delta_one_step: " << fman->GetDeltaOneStep()/mm << G4endl;
  //G4cout << "largest_step: " << fman-> << G4endl;

  G4VIntegrationDriver *driver = fman->GetChordFinder()->GetIntegrationDriver();
  driver->SetVerboseLevel(3);

  G4MagIntegratorStepper *stepper = driver->GetStepper();
  G4cout << "stepper: " << typeid(*stepper).name() << G4endl;

  G4EquationOfMotion *eq = stepper->GetEquationOfMotion();
  G4cout << "equation: " << typeid(*eq).name() << G4endl;

}//PrintField

//_____________________________________________________________________________
G4VisAttributes *BeamDipole::ColorDecoder(GeoParser *geo) {

  G4String col("0:1:0:0.7"); // red:green:blue:alpha 
  geo->GetOptS(fNam, "vis", col);

  char_separator<char> sep(":");
  tokenizer< char_separator<char> > clin(col, sep);
  tokenizer< char_separator<char> >::iterator it = clin.begin();

  stringstream st;
  for(int i=0; i<4; i++) {
    st << *(it++) << " ";
  }

  G4double red=0, green=0, blue=0, alpha=0;
  st >> red >> green >> blue >> alpha;

  G4VisAttributes *vis = new G4VisAttributes();
  if(alpha < 1.1) {
    vis->SetColor(red, green, blue, alpha);
    vis->SetForceSolid(true);
  } else {
    vis->SetColor(red, green, blue);
    vis->SetForceAuxEdgeVisible(true);
  }

  return vis;

}//ColorDecoder
































