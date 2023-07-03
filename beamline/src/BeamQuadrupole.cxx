
//_____________________________________________________________________________
//
// beam quadrupole magnet, its vessel absorbes all particles, no secondaries
//
//_____________________________________________________________________________

//C++
#include "math.h"

//Boost
#include <boost/tokenizer.hpp>

//Geant
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4QuadrupoleMagField.hh"
#include "G4FieldManager.hh"
#include "G4RotationMatrix.hh"
#include "G4VIntegrationDriver.hh"
#include "G4ChordFinder.hh"
#include "G4DormandPrinceRK78.hh"
#include "G4DormandPrince745.hh"
#include "G4DormandPrinceRK56.hh"
#include "G4MagErrorStepper.hh"
#include "G4ClassicalRK4.hh"
//#include "G4TDormandPrince45.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4IntegrationDriver.hh"

//local classes
#include "GeoParser.h"
#include "BeamQuadrupole.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
BeamQuadrupole::BeamQuadrupole(G4String nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "BeamQuadrupole: " << fNam << G4endl;

  //magnet center along x and z, mm
  G4double xpos = 0, zpos = 0;
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //length, mm
  G4double length = geo->GetD(fNam, "length")*mm;

  //inner radii closer to the IP (r1) and further from the IP (r2), m
  G4double r1 = 10*mm, r2 = 10*mm;
  geo->GetOptD(fNam, "r1", r1, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "r2", r2, GeoParser::Unit(mm));

  //radial thickness, mm
  G4double dr = r2;
  geo->GetOptD(fNam, "dr", dr, GeoParser::Unit(mm));

  //magnet top volume holding the core and vessel
  G4Cons *shape_top = new G4Cons(fNam+"_top", 0, r2+dr, 0, r1+dr, length/2, 0, 360*deg);
  G4Material *mat_core = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume *vol_top = new G4LogicalVolume(shape_top, mat_core, fNam+"_top");
  vol_top->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //polar angle along y axis for the magnet
  G4double theta = 0;
  geo->GetOptD(fNam, "theta", theta, GeoParser::Unit(rad));
  G4RotationMatrix rot_top(G4ThreeVector(0, 1, 0), theta); //CLHEP::HepRotation

  //magnet top volume in global top
  G4ThreeVector pos(xpos, 0, zpos);
  G4Transform3D transform(rot_top, pos); //HepGeom::Transform3D
  new G4PVPlacement(transform, vol_top, fNam+"_top", top, false, 0);

  //magnet inner core
  G4String nam_mag = fNam+"_mag";
  G4Cons *shape_mag = new G4Cons(nam_mag, 0, r2, 0, r1, length/2, 0, 360*deg);

  //core logical volume
  G4LogicalVolume *vol_mag = new G4LogicalVolume(shape_mag, mat_core, nam_mag);
  vol_mag->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //quadrupole field inside the core
  G4double grad = geo->GetD(fNam, "grad") * tesla/meter; // field gradient, T/m
  G4QuadrupoleMagField *field = 0x0;
  G4double angle = 0; // field angle, rad
  if( geo->GetOptD(fNam, "angle", angle, GeoParser::Unit(rad)) ) {
    G4RotationMatrix *rot = new G4RotationMatrix(G4ThreeVector(0, 0, 1), angle); //is typedef to CLHEP::HepRotation
    field = new G4QuadrupoleMagField(grad, G4ThreeVector(0, 0, 0), rot);
  } else {
    field = new G4QuadrupoleMagField(grad);
  }

  G4FieldManager *fman = new G4FieldManager();
  fman->SetDetectorField(field);

  //alternative stepper
  //fman->SetMinimumEpsilonStep(5e-2);
  //fman->SetMaximumEpsilonStep(0.1);
  //G4ClassicalRK4 *stepper = new G4ClassicalRK4(new G4Mag_UsualEqRhs(field));
  //G4VIntegrationDriver *driver = new G4IntegrationDriver<G4ClassicalRK4>(5e-2*mm, stepper);
  //G4ChordFinder *finder = new G4ChordFinder(driver);
  //fman->SetChordFinder(finder);

  fman->CreateChordFinder(field);
  //PrintField(fman);
  vol_mag->SetFieldManager(fman, true);

  //put the inner core to the top magnet volume
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), vol_mag, nam_mag, vol_top, false, 0);

  //vessel shape
  G4Cons *shape = new G4Cons(fNam, r2, r2+dr, r1, r1+dr, length/2, 0, 360*deg);

  //vessel logical volume
  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
  G4LogicalVolume *vol = new G4LogicalVolume(shape, mat, fNam);

  //vessel visibility
  vol->SetVisAttributes(ColorDecoder(geo));

  //put the vessel cone to the top volume
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), vol, fNam, vol_top, false, 0);

}//ConeAperture

//_____________________________________________________________________________
void BeamQuadrupole::PrintField(G4FieldManager *fman) {

  G4ChordFinder *finder = fman->GetChordFinder();

  G4cout << "BeamQuadrupole:" << G4endl;
  G4cout << "eps_min: " << fman->GetMinimumEpsilonStep()/mm << G4endl;
  G4cout << "eps_max: " << fman->GetMaximumEpsilonStep()/mm << G4endl;
  //G4cout << "min_chord_step: " << fman-> << G4endl;
  G4cout << "delta_chord: " << finder->GetDeltaChord()/mm << G4endl;
  G4cout << "delta_intersection: " << fman->GetDeltaIntersection()/mm << G4endl;
  G4cout << "delta_one_step: " << fman->GetDeltaOneStep()/mm << G4endl;
  //G4cout << "largest_step: " << fman-> << G4endl;

  G4cout << "finder: " << typeid(*fman->GetChordFinder()).name() << G4endl;

  G4VIntegrationDriver *driver = fman->GetChordFinder()->GetIntegrationDriver();
  G4cout << "driver: " << typeid(*driver).name() << G4endl;
  //driver->SetVerboseLevel(3);

  G4MagIntegratorStepper *stepper = driver->GetStepper();
  G4cout << "stepper: " << typeid(*stepper).name() << G4endl;
  //G4DormandPrince745
  //G4DormandPrinceRK78
  //G4DormandPrinceRK56
  //G4MagErrorStepper
  //G4ClassicalRK4
  //G4MagIntegratorStepper
  //G4TDormandPrince45
  //G4cout << "cast: " << dynamic_cast<G4TDormandPrince45*>(stepper) << G4endl;
  //G4cout << "cast: " << dynamic_cast<G4TDormandPrince45<G4Mag_UsualEqRhs>*>(stepper) << G4endl;

  G4EquationOfMotion *eq = stepper->GetEquationOfMotion();
  G4cout << "equation: " << typeid(*eq).name() << G4endl;

}//PrintField

//_____________________________________________________________________________
G4VisAttributes *BeamQuadrupole::ColorDecoder(GeoParser *geo) {

  G4String col("0:0:1:0.7"); // red:green:blue:alpha 
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












