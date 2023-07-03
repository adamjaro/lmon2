
//_____________________________________________________________________________
//
// tungsten/scintillator sampling calorimeter following
// ZEUS Beam Pipe Calorimeter (BPC), Nucl.Instrum.Meth.A 565 (2006) 572-588
//
//_____________________________________________________________________________

//C++
#include <map>
#include <vector>
#include <string>

//ROOT
#include "TTree.h"

//Geant
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4VisAttributes.hh"

//local headers
#include "CaloBPC.h"
#include "GeoParser.h"
#include "ColorDecoder.h"

using namespace std;

//_____________________________________________________________________________
CaloBPC::CaloBPC(const G4String& nam, GeoParser *geo, G4LogicalVolume *top): Detector(),
  G4VSensitiveDetector(nam), fNam(nam) {

  G4cout << "CaloBPC: " << fNam << G4endl;

  //module location in top
  G4double mod_zpos = geo->GetD(fNam, "zpos")*mm; // center position in z, mm
  G4double mod_xpos = 0; // center position in x, mm
  geo->GetOptD(fNam, "xpos", mod_xpos, GeoParser::Unit(mm));
  G4double mod_rot_y = 0; // module rotation in x-z plane by rotation along y, rad
  geo->GetOptD(fNam, "rot_y", mod_rot_y, GeoParser::Unit(rad));

  //module internal dimensions
  G4double modxy = 138*mm; // module transverse size
  geo->GetOptD(fNam, "modxy", modxy, GeoParser::Unit(mm));
  G4double abso_z = 3.5*mm; // absorber thickness along z
  G4double scin_z = 2.6*mm; // scintillator thickness along z

  G4int nlay = 26; // number of absorber and scintillator layers in module

  fNscin = 17; // number of scintillator strips, data member for scintillator addressing
  G4double scin_x = 7.9*mm; // scintillator width in x
  G4double scin_spacing = 8*mm; // spacing of scintillator strips
  G4double scin_ofs = 1*mm; // first scintillator offset from the module edge

  geo->GetOptI(fNam, "nscin", fNscin);
  geo->GetOptD(fNam, "scin_x", scin_x, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "scin_spacing", scin_spacing, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "scin_ofs", scin_ofs, GeoParser::Unit(mm));

  //geometry for module holding the absorbers and scintillators
  G4double modz = nlay*(abso_z + scin_z); // module length along z
  G4cout << "  " << fNam << ", modz (mm): " << modz << G4endl;
  G4String modnam = fNam+"_mod"; //module name

  //box shape for the module
  G4Box *mods = new G4Box(modnam, modxy/2, modxy/2, modz/2);

  //module default material
  G4Material *mod_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  //module logical volume
  G4LogicalVolume *modv = new G4LogicalVolume(mods, mod_mat, modnam);

  //module visibility
  ColorDecoder modvis("0:0:1:2");
  modv->SetVisAttributes(modvis.MakeVis(geo, fNam, "mod_vis"));

  //put the module to its mother volume
  G4RotationMatrix modrot(G4ThreeVector(0, 1, 0), mod_rot_y); //is typedef to CLHEP::HepRotation
  G4ThreeVector modpos(mod_xpos, 0, mod_zpos);
  G4Transform3D modtrans(modrot, modpos); // is HepGeom::Transform3D

  //mother volume for the module
  G4LogicalVolume *mother_vol = top;
  G4String mother_nam;
  if( geo->GetOptS(fNam, "place_into", mother_nam) ) {
    mother_vol = GetMotherVolume(mother_nam, top);
  }

  new G4PVPlacement(modtrans, modv, modnam, mother_vol, false, 0);

  //layer holding the absorber and scintillator strips
  G4double lay_z = abso_z+scin_z; // layer thickness along z
  G4String lay_nam = fNam+"_layer"; //layer name

  //layer shape
  G4Box *lay_shape = new G4Box(lay_nam, modxy/2, modxy/2, lay_z/2);

  //layer logical volume
  G4LogicalVolume *lay_vol = new G4LogicalVolume(lay_shape, mod_mat, lay_nam);
  lay_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //put the layers to the module
  for(G4int ilay=0; ilay<nlay; ilay++) {
    G4double lay_zpos = -ilay*lay_z + modz/2 - lay_z/2;

    //rotation by pi/2 for every second layer
    G4RotationMatrix *lay_rot = 0x0;
    if( ilay%2 != 0 ) {
      lay_rot = new G4RotationMatrix();
      lay_rot->rotateZ(-90*deg);
    }

    new G4PVPlacement(lay_rot, G4ThreeVector(0, 0, lay_zpos), lay_vol, lay_nam, modv, false, ilay);
  }

  //absorber plates
  G4String abso_nam = fNam+"_abso"; //absorber name

  //absorber plate
  G4Box *abso_shape = new G4Box(abso_nam, modxy/2, modxy/2, abso_z/2);

  //absorber material
  G4Material *abso_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

  //absorber logical volume
  G4LogicalVolume *abso_vol = new G4LogicalVolume(abso_shape, abso_mat, abso_nam);

  //absorber visibility
  ColorDecoder abso_vis("1:0:0:0.5");
  abso_vol->SetVisAttributes(abso_vis.MakeVis(geo, fNam, "abso_vis"));

  //put the absorber to the layer
  new G4PVPlacement(0, G4ThreeVector(0, 0, scin_z/2.), abso_vol, abso_nam, lay_vol, false, 0);

  //scintillator strips
  G4Box *scin_shape = new G4Box(fNam, scin_x/2, modxy/2, scin_z/2);

  //polystyrene for SCSN-38
  G4Material *scin_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE");

  //scintillator logical volume
  G4LogicalVolume *scin_vol = new G4LogicalVolume(scin_shape, scin_mat, fNam);

  //scintillator visibility
  ColorDecoder scin_vis("1:1:0:0.5");
  scin_vol->SetVisAttributes(scin_vis.MakeVis(geo, fNam, "scin_vis"));

  //put scintillator strips to the layer
  for(G4int ix=0; ix<fNscin; ix++) {
    //vertical strips along x
    G4double xscin_pos = ix*scin_spacing - modxy/2 + scin_ofs + scin_spacing/2;
    new G4PVPlacement(0, G4ThreeVector(xscin_pos, 0, -scin_z/2.), scin_vol, fNam, lay_vol, false, ix);

    //scintillator positions for hits
    //fHits.SetScinPos(ix, xscin_pos);
  }

}//CaloBPC

//_____________________________________________________________________________
G4bool CaloBPC::ProcessHits(G4Step *step, G4TouchableHistory*) {

  //deposited energy
  G4double edep_in_step = step->GetTotalEnergyDeposit()/GeV;

  //consider only steps with energy deposit
  //if( edep_in_step < 1e-12 ) return true;

  //strip location
  const G4TouchableHandle& hnd = step->GetPreStepPoint()->GetTouchableHandle();

  G4int istrip = hnd->GetCopyNumber(); // strip index on the layer
  G4int ilay = hnd->GetCopyNumber(1); // layer index in the module

  //global strip position
  G4ThreeVector origin(0, 0, 0);
  G4ThreeVector gpos = hnd->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);
  G4double x = gpos.x()/mm;
  G4double y = gpos.y()/mm;
  G4double z = gpos.z()/mm;

  //hit at given strip location
  CaloBPCHits::Hit& hit = fHits.ConstructedAt(make_pair(istrip, ilay), CaloBPCHits::Hit(istrip, ilay, x, y, z));

  //add deposited energy in the hit
  hit.en += edep_in_step;

  return true;

}//ProcessHits

//_____________________________________________________________________________
void CaloBPC::FinishEvent() {

  fHits.FinishEvent();

}//FinishEvent

//_____________________________________________________________________________
void CaloBPC::ClearEvent() {

  fHits.ClearEvent();

}//ClearEvent

//_____________________________________________________________________________
void CaloBPC::CreateOutput(TTree *tree) {

  fHits.CreateOutput(fNam, tree);

}//CreateOutput

//_____________________________________________________________________________
G4LogicalVolume* CaloBPC::GetMotherVolume(G4String mother_nam, G4LogicalVolume *top) {

  for(size_t i=0; i<top->GetNoDaughters(); i++) {

    G4LogicalVolume *dv = top->GetDaughter(i)->GetLogicalVolume();

    if( dv->GetName() == mother_nam ) {
      return dv;
    }
  }

  return 0x0;

}//GetMotherVolume





























