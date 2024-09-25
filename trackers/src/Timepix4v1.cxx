
//_____________________________________________________________________________
//
// Timepix4, implementation version 1
//
// Geometry parameters:
//  dxy (mm): pixel size in x and y
//  nx, ny: number of pixels along x and y
//
// Optional geometry parameters:
//  place_into: name of volume to place the Timepix4v1 to
//  xpos, ypos, zpos (mm): position in the volume where the detector is placed
//  lay_Cu_dz, lay_PCB_circuit_dz, lay_PCB_TSV_dz
//  lay_ASIC_dz, lay_Si_dz (mm): thickness of individual layers
//  lay_vis, pix_vis: visibility of the entire layer and for individual pixels
//  Cu_mat_name, Cu_vis: material name and visibility for Cu back layer
//  PCB_circuit_mat_name, PCB_circuit_vis: material and visibility for PCB
//  TSV_mat_name, TSV_vis: material and visibility for TSV PCB layer
//  ASIC_mat_name, ASIC_vis: material and visibility for ASIC layer
//
//_____________________________________________________________________________

//C++

//ROOT

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

//local classes
#include "GeoParser.h"
#include "ColorDecoder.h"
#include "MCParticleAction.h"
#include "Timepix4v1.h"
#include "LogicVolFinder.h"

//_____________________________________________________________________________
Timepix4v1::Timepix4v1(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
   Detector(), G4VSensitiveDetector(nam), fNam(nam), fGeo(geo), fStack(0x0) {

  G4cout << "Timepix4v1: " << fNam << G4endl;

  //transverse dimensions
  G4double dxy = geo->GetD(fNam, "dxy")*mm; // full pixel size in x and y, mm
  G4int nx = geo->GetI(fNam, "nx"); // number of pixels along x
  G4int ny = geo->GetI(fNam, "ny"); // number of pixels along y

  //layers from the back (lower z)
  G4double lay_Cu_dz = 0.1*mm; // backplane Copper layer
  G4double lay_PCB_circuit_dz = 0.2*mm; // circuit PCB
  G4double lay_PCB_TSV_dz = 0.1*mm; // Through-silicon-Vias
  G4double lay_ASIC_dz = 0.1*mm; // Timepix4 ASIC
  G4double lay_Si_dz = 0.06*mm; // Silicon sensor
  geo->GetOptD(fNam, "lay_Cu_dz", lay_Cu_dz, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "lay_PCB_circuit_dz", lay_PCB_circuit_dz, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "lay_PCB_TSV_dz", lay_PCB_TSV_dz, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "lay_ASIC_dz", lay_ASIC_dz, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "lay_Si_dz", lay_Si_dz, GeoParser::Unit(mm));

  //total layer thickness and transverse size
  G4double lay_dz = lay_Cu_dz + lay_PCB_circuit_dz + lay_PCB_TSV_dz + lay_ASIC_dz + lay_Si_dz;
  G4double lay_dx = nx*dxy;
  G4double lay_dy = ny*dxy;

  G4cout << "  " << fNam <<", lay_dx, dy, dz (mm): " << lay_dx << " " << lay_dy << " " << lay_dz << G4endl;

  //main layer volume to contain all individual layers
  G4Box *main_shape = new G4Box(nam+"_main", lay_dx/2., lay_dy/2., lay_dz/2.);
  G4Material *main_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *main_vol = new G4LogicalVolume(main_shape, main_mat, main_shape->GetName());

  //main layer visibility
  ColorDecoder lay_vis("0:0:1:2");
  main_vol->SetVisAttributes(lay_vis.MakeVis(geo, fNam, "lay_vis"));

  //optional mother volume for the layer (top assigned as default)
  G4LogicalVolume *mother_vol = LogicVolFinder::GetMotherVolume("place_into", top, geo, fNam);
  //case of two levels of volume to place, volume 'place_into2' inside 'place_into1'
  if( geo->HasParameter(fNam, "place_into2") ) {
    mother_vol = LogicVolFinder::GetMotherVolume2("place_into1", "place_into2", top, geo, fNam);
  }

  //center position for the main layer in x, y and z, mm
  G4double xpos=0, ypos=0, zpos=0;
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "ypos", ypos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //layer main volume in its mother volume
  new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos), main_vol, main_vol->GetName(), mother_vol, false, 0);

  //tracker sensor layer (sensitive pixels)
  G4LogicalVolume *lay_sens = MakeSensorLayer(nx, ny, dxy, lay_Si_dz, lay_dx, lay_dy);

  //tracker layer in main volume
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*(lay_dz-lay_Si_dz)), lay_sens, lay_sens->GetName(), main_vol, false, 0);

  //PCB material as FR4, made of quartz and epoxy, based on ExN03

  //Epoxy for FR4
  if( !G4Material::GetMaterial("Timepix4v1_Epoxy", false) ) {
    G4Material *epoxy = new G4Material("Timepix4v1_Epoxy", 1.3*g/cm3, 3);
    epoxy->AddElement(G4NistManager::Instance()->FindOrBuildElement("H"), 44);
    epoxy->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 15);
    epoxy->AddElement(G4NistManager::Instance()->FindOrBuildElement("O"), 7);
    // https://allpix-squared-forum.web.cern.ch/t/new-passive-materials/144
  }

  //FR4 material
  if( !G4Material::GetMaterial("Timepix4v1_FR4", false) ) {
    G4Material *fr4 = new G4Material("Timepix4v1_FR4", 1.86*g/cm3, 2);
    fr4->AddMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE"), 0.528); // SiO_2
    fr4->AddMaterial(G4Material::GetMaterial("Timepix4v1_Epoxy"), 0.472);
    // https://www.phenix.bnl.gov/~suhanov/ncc/geant/rad-source/src/ExN03DetectorConstruction.cc
  }

  //MakeMaterialLayer(G4double layx, G4double layy, G4double dz,
  //G4String mat_name, G4String mat_cmd, G4String vis, G4String vis_cmd) {

  //Copper back layer
  G4LogicalVolume *cu_lay = MakeMaterialLayer(lay_dx, lay_dy, lay_Cu_dz, "G4_Cu", "Cu_mat_name", "1:0:0:3", "Cu_vis");
  new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*lay_dz+0.5*lay_Cu_dz), cu_lay, cu_lay->GetName(), main_vol, false, 0);

  //PCB circuit layer
  G4LogicalVolume *pcb_lay = MakeMaterialLayer(lay_dx, lay_dy, lay_PCB_circuit_dz,
    "Timepix4v1_FR4", "PCB_circuit_mat_name", "0:1:0:3", "PCB_circuit_vis");
  G4double pcb_lay_zpos = -0.5*lay_dz + lay_Cu_dz + 0.5*lay_PCB_circuit_dz;
  new G4PVPlacement(0, G4ThreeVector(0, 0, pcb_lay_zpos), pcb_lay, pcb_lay->GetName(), main_vol, false, 0);

  //PCB TSV layer
  G4LogicalVolume *tsv_lay = MakeMaterialLayer(lay_dx, lay_dy, lay_PCB_TSV_dz,
    "Timepix4v1_FR4", "TSV_mat_name", "1:1:0:3", "TSV_vis");
  G4double tsv_lay_zpos = -0.5*lay_dz + lay_Cu_dz + lay_PCB_circuit_dz + 0.5*lay_PCB_TSV_dz;
  new G4PVPlacement(0, G4ThreeVector(0, 0, tsv_lay_zpos), tsv_lay, tsv_lay->GetName(), main_vol, false, 0);

  //ASIC layer
  G4LogicalVolume *asic_lay = MakeMaterialLayer(lay_dx, lay_dy, lay_ASIC_dz,
    "G4_Si", "ASIC_mat_name", "0:1:1:3", "ASIC_vis");
  G4double asic_lay_zpos = -0.5*lay_dz + lay_Cu_dz + lay_PCB_circuit_dz + lay_PCB_TSV_dz + 0.5*lay_ASIC_dz;
  new G4PVPlacement(0, G4ThreeVector(0, 0, asic_lay_zpos), asic_lay, asic_lay->GetName(), main_vol, false, 0);

}//Timepix4v1

//_____________________________________________________________________________
G4LogicalVolume* Timepix4v1::MakeSensorLayer(G4int nx, G4int ny, G4double dxy, G4double dz, G4double layx, G4double layy) {

  //tracker sensor layer
  G4Box *lays = new G4Box(fNam+"_Si_lay", layx/2., layy/2., dz/2.);
  G4Material *laym = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *layv = new G4LogicalVolume(lays, laym, lays->GetName());
  layv->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //row holding individual pixels
  G4Box *row_shape = new G4Box(fNam+"_Si_row", layx/2., dxy/2., dz/2.);
  G4LogicalVolume *row_vol = new G4LogicalVolume(row_shape, laym, row_shape->GetName());
  row_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //rows in tracker layer
  new G4PVReplica(row_shape->GetName(), row_vol, layv, kYAxis, ny, dxy);

  //sensitive pixels in each row
  G4Box *pix_shape = new G4Box(fNam, dxy/2., dxy/2., dz/2.);

  //silicon material for pixels
  G4Material *pix_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  //pixel volume
  G4LogicalVolume *pix_vol = new G4LogicalVolume(pix_shape, pix_mat, fNam);
  ColorDecoder pix_vis("1:0:0:3");
  pix_vol->SetVisAttributes(pix_vis.MakeVis(fGeo, fNam, "pix_vis"));

  //pixels in row
  new G4PVReplica(fNam, pix_vol, row_vol, kXAxis, nx, dxy);

  return layv;

}//MakeSensorLayer

//_____________________________________________________________________________
G4LogicalVolume* Timepix4v1::MakeMaterialLayer(G4double layx, G4double layy, G4double dz,
   G4String mat_name, G4String mat_cmd, G4String vis, G4String vis_cmd) {

  //construction material layer

  //layer material
  G4String lay_mat_name = mat_name;
  fGeo->GetOptS(fNam, mat_cmd, lay_mat_name);

  //select the layer material
  G4Material *lay_mat=0x0;
  if( lay_mat_name.find("Timepix4v1_") != std::string::npos ) {
    lay_mat = G4Material::GetMaterial(lay_mat_name); // local material
  } else {
    lay_mat = G4NistManager::Instance()->FindOrBuildMaterial(lay_mat_name); // nist material
  }

  G4cout << "  " << fNam << ", layer material: " << lay_mat->GetName() << G4endl;

  //material layer volume
  G4Box *lays = new G4Box(fNam+"_"+lay_mat_name, layx/2., layy/2., dz/2.);
  G4LogicalVolume *layv = new G4LogicalVolume(lays, lay_mat, lays->GetName());

  //layer visibility
  ColorDecoder lay_vis(vis);
  layv->SetVisAttributes(lay_vis.MakeVis(fGeo, fNam, vis_cmd));

  return layv;

}//MakeMaterialLayer

//_____________________________________________________________________________
G4bool Timepix4v1::ProcessHits(G4Step *step, G4TouchableHistory*) {

  //G4cout << "Timepix4v1::ProcessHits" << G4endl;

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

  //G4cout << "ipix, irow, x, y, z (mm): " << ipix << " " << irow << " " << x << " " << y << " " << z << G4endl;

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

  //hit at the given pixel location
  TrkPlaneBasicHits::Hit& hit = fHits.ConstructedAt( std::make_pair(ipix, irow),
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
void Timepix4v1::CreateOutput(TTree *tree) {

  fHits.CreateOutput(fNam, tree);

  fStack = static_cast<const MCParticleAction*>( G4RunManager::GetRunManager()->GetUserTrackingAction() );

}//CreateOutput

//_____________________________________________________________________________
void Timepix4v1::ClearEvent() {

  fHits.ClearEvent();

}//ClearEvent

//_____________________________________________________________________________
void Timepix4v1::FinishEvent() {

  fHits.FinishEvent();

}//FinishEvent




















