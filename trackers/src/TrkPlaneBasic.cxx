
//_____________________________________________________________________________
//
// Si tracking layer, basic implementation
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
#include "G4PVReplica.hh"
#include "G4RunManager.hh"

//local headers
#include "TrkPlaneBasic.h"
#include "GeoParser.h"
#include "ColorDecoder.h"
#include "MCParticleAction.h"

using namespace std;

//_____________________________________________________________________________
TrkPlaneBasic::TrkPlaneBasic(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), G4VSensitiveDetector(nam), fNam(nam) {

  G4cout << "TrkPlaneBasic: " << fNam << G4endl;

  G4double dz = geo->GetD(fNam, "dz")*mm; // full layer thickness along z, mm
  G4double dxy = geo->GetD(fNam, "dxy")*mm; // full pixel size in x and y, mm
  G4int nx = geo->GetI(fNam, "nx"); // number of pixels along x
  G4int ny = geo->GetI(fNam, "ny"); // number of pixels along y

  //G4cout << nx << " " << ny << " " << geo->GetS(fNam, "nx") << G4endl;

  //tracker layer
  G4double layx = nx*dxy;
  G4double layy = ny*dxy;
  G4cout << "  " << nam << ": layx, layy, dxy (mm): " << layx << ", " << layy << ", " << dxy << G4endl;
  G4Box *lays = new G4Box(nam+"_lay", layx/2., layy/2., dz/2.);

  //layer volume
  G4Material *laym = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4LogicalVolume *layv = new G4LogicalVolume(lays, laym, lays->GetName());

  //layer visibility
  ColorDecoder lay_vis("1:0:0:0.3");
  layv->SetVisAttributes(lay_vis.MakeVis(geo, fNam, "lay_vis"));

  //mother volume for the layer
  G4LogicalVolume *mother_vol = top;
  G4String mother_nam;
  if( geo->GetOptS(fNam, "place_into", mother_nam) ) {
    mother_vol = GetMotherVolume(mother_nam, top);
  }

  //center position for the layer in x, y and z, mm
  G4double xpos=0, ypos=0, zpos=0;
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "ypos", ypos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));
  G4ThreeVector lay_pos(xpos, ypos, zpos);

  //layer rotation along y, rad
  G4double rotate_y = 0; // rotation along y axis
  geo->GetOptD(fNam, "rotate_y", rotate_y, GeoParser::Unit(rad));
  G4RotationMatrix lay_rot(G4ThreeVector(0, 1, 0), rotate_y); //CLHEP::HepRotation

  //layer in its mother volume
  G4Transform3D lay_transform(lay_rot, lay_pos); //HepGeom::Transform3D
  new G4PVPlacement(lay_transform, layv, lays->GetName(), mother_vol, false, 0);

  //row holding individual pixels
  G4Box *row_shape = new G4Box(nam+"_row", layx/2., dxy/2., dz/2.);
  G4LogicalVolume *row_vol = new G4LogicalVolume(row_shape, laym, row_shape->GetName());
  row_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //rows in tracker layer
  new G4PVReplica(row_shape->GetName(), row_vol, layv, kYAxis, ny, dxy);

  //sensitive pixels in each row
  G4Box *pix_shape = new G4Box(nam, dxy/2., dxy/2., dz/2.);

  //silicon material for pixels
  G4Material *pix_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  //pixel volume
  G4LogicalVolume *pix_vol = new G4LogicalVolume(pix_shape, pix_mat, nam);
  ColorDecoder pix_vis("1:0:0:3");
  pix_vol->SetVisAttributes(pix_vis.MakeVis(geo, fNam, "pix_vis"));

  //pixels in row
  new G4PVReplica(nam, pix_vol, row_vol, kXAxis, nx, dxy);

}//TrkPlaneBasic

//_____________________________________________________________________________
G4bool TrkPlaneBasic::ProcessHits(G4Step *step, G4TouchableHistory*) {

  //pixel location
  const G4TouchableHandle& hnd = step->GetPreStepPoint()->GetTouchableHandle();
  G4int ipix = hnd->GetCopyNumber(); // pixel index in the row
  G4int irow = hnd->GetCopyNumber(1); // row index in the layer

  //global pixel position
  G4ThreeVector origin(0, 0, 0);
  G4ThreeVector gpos = hnd->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);
  G4double x = gpos.x()/mm;
  G4double y = gpos.y()/mm;
  G4double z = gpos.z()/mm;

  //deposited energy in step
  G4double en = step->GetTotalEnergyDeposit()/keV;

  //track ID and PDG in the step
  G4Track *track = step->GetTrack();
  G4int itrk = track->GetTrackID();
  G4int pdg = track->GetParticleDefinition()->GetPDGEncoding();

  //track by primary particle
  G4bool is_prim = track->GetParentID() > 0 ? false : true;

  //ID of primary particle for the hit
  Int_t prim_id = fStack->GetPrimaryID(itrk);

  //G4cout << "TrkPlaneBasic::ProcessHits: " << track->GetTrackID() << " " << track->GetParentID() << " " << is_prim << endl;

  //hit at the given pixel location
  TrkPlaneBasicHits::Hit& hit = fHits.ConstructedAt( make_pair(ipix, irow),
    TrkPlaneBasicHits::Hit(ipix, irow, x, y, z, itrk, pdg, is_prim, prim_id) );

  //add deposited energy in the hit
  hit.en += en;

  //lowest track ID associated with the hit and its PDG code
  if( itrk < hit.itrk ) {
    hit.itrk = itrk;
    hit.pdg = pdg;
  }

  //update primary flag if not set
  if( hit.is_prim == kFALSE ) {
    hit.is_prim = is_prim;
  }

  //lowest primary ID for the hit
  if( prim_id < hit.prim_id ) {

    hit.prim_id = prim_id;
  }

  return true;

}//ProcessHits

//_____________________________________________________________________________
void TrkPlaneBasic::CreateOutput(TTree *tree) {

  fHits.CreateOutput(fNam, tree);

  fStack = static_cast<const MCParticleAction*>( G4RunManager::GetRunManager()->GetUserTrackingAction() );

}//CreateOutput

//_____________________________________________________________________________
void TrkPlaneBasic::ClearEvent() {

  fHits.ClearEvent();

}//ClearEvent

//_____________________________________________________________________________
void TrkPlaneBasic::FinishEvent() {

  //G4cout << "FinishEvent for " << fNam << G4endl;

  fHits.FinishEvent();

}//FinishEvent

//_____________________________________________________________________________
G4LogicalVolume* TrkPlaneBasic::GetMotherVolume(G4String mother_nam, G4LogicalVolume *top) {

  for(size_t i=0; i<top->GetNoDaughters(); i++) {

    G4LogicalVolume *dv = top->GetDaughter(i)->GetLogicalVolume();

    if( dv->GetName() == mother_nam ) {
      return dv;
    }
  }

  return 0x0;

}//GetMotherVolume
































