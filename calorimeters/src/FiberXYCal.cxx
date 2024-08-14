
//_____________________________________________________________________________
//
// W Scifi calorimeter with fibers perpendicular to shower axis
// in alternating layers, geometry by Aranya Giri and Simon Gardner
//
// Optical properties follow Nucl.Instrum.Meth.A 386 (1997) 397-408
// (H1 lead-fiber calorimeter).
//
// Preliminary scintillation properties are taken from example
// examples/extended/optical/wls/
//
// mod_xy (mm): module transverze size (and cell longitudinal size), must equal to cell_xy * ncells
// mod_z (mm): module length along z, must equal to nlay * cell_xy
// nlay: number of layers in calorimeter, must be even
// cell_xy (mm): cell size along x and y
// ncells: number of cells in layer
// fiber_core_r (mm): radius for fiber core (cladding and groove thickness defaults in the code)
// fiber_dx (mm): fiber spacing (their mutual distances)
//
//_____________________________________________________________________________

//C++
#include <vector>
#include <math.h>

//ROOT
#include "TTree.h"

//Geant
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalMaterialProperties.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4PVReplica.hh"

//local classes
#include "GeoParser.h"
#include "ColorDecoder.h"
#include "OpSiDet.h"
#include "FiberXYCal.h"

//_____________________________________________________________________________
FiberXYCal::FiberXYCal(const G4String& nam, GeoParser *geo, G4LogicalVolume *top): Detector(),
   G4VSensitiveDetector(nam), fNam(nam) {

  G4cout << "FiberXYCal: " << fNam << G4endl;

  //layer holding individual cells
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm; // cell size
  G4int ncells = geo->GetI(fNam, "ncells"); // number of cells in layer
  G4double mod_xy = geo->GetD(fNam, "mod_xy")*mm; // must be equal to cell_xy * ncells

  //layer volume
  G4Box *layer_shape = new G4Box(fNam+"_lay", mod_xy/2, cell_xy/2, mod_xy/2);
  G4Material *air_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *layer_vol = new G4LogicalVolume(layer_shape, air_mat, layer_shape->GetName());

  //layer visibility
  ColorDecoder layer_vis("1:0:0:3");
  layer_vol->SetVisAttributes(layer_vis.MakeVis(geo, fNam, "layer_vis"));

  //cell volume holding all structure of absorber and fibers
  G4LogicalVolume *cell_vol = MakeCell(geo);

  //cells in layer
  new G4PVReplica(cell_vol->GetName(), cell_vol, layer_vol, kXAxis, ncells, cell_xy);

  //calorimeter module
  G4int nlay = geo->GetI(fNam, "nlay"); // number of layers in calorimeter module (even number)
  G4double mod_z = geo->GetD(fNam, "mod_z")*mm; // module size along z

  //module volume
  G4Box *mod_shape = new G4Box(fNam+"_mod", mod_xy/2, mod_xy/2, mod_z/2);
  G4LogicalVolume *mod_vol = new G4LogicalVolume(mod_shape, air_mat, mod_shape->GetName());

  //module visibility
  ColorDecoder mod_vis("1:0:0:3");
  mod_vol->SetVisAttributes(mod_vis.MakeVis(geo, fNam, "mod_vis"));

  //alternating rotation for layers in the module
  G4RotationMatrix *rotx = new G4RotationMatrix();
  rotx->rotateX(-90*deg);
  G4RotationMatrix *rotzx = new G4RotationMatrix();
  rotzx->rotateZ(-90*deg);
  rotzx->rotateX(90*deg);

  //layers by two in alternating rotation in calorimeter module
  G4int layer_cnt = 0;
  //layer loop
  for(G4int i=0; i<nlay/2; i++) {

    //alternating rotation
    G4RotationMatrix *rot = rotx;
    if( i%2 != 0 ) rot = rotzx;

    //two layers for a given rotation
    new G4PVPlacement(rot, G4ThreeVector(0, 0, -0.5*mod_z+(0.5+layer_cnt)*cell_xy),
      layer_vol, layer_vol->GetName(), mod_vol, false, layer_cnt++);

    new G4PVPlacement(rot, G4ThreeVector(0, 0, -0.5*mod_z+(0.5+layer_cnt)*cell_xy),
      layer_vol, layer_vol->GetName(), mod_vol, false, layer_cnt++);

  }//layer loop

  //module center position
  G4double xpos=0, ypos=0, zpos=0;
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "ypos", ypos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //calorimeter module in its mother volume
  G4LogicalVolume *mother_vol = top;
  G4String mother_nam;
  if( geo->GetOptS(fNam, "place_into", mother_nam) ) {
    mother_vol = GetMotherVolume(mother_nam, top);
  }
  new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos), mod_vol, mod_vol->GetName(), mother_vol, false, 0);

}//FiberXYCal

//_____________________________________________________________________________
G4LogicalVolume* FiberXYCal::MakeCell(GeoParser *geo) {

  //cell size in x-y
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm;

  //cell size in z, same as module in x-y
  G4double cell_z = geo->GetD(fNam, "mod_xy")*mm;

  G4cout << "FiberXYCal::MakeCell, cell_xy, z (mm): " << cell_xy << ", " << cell_z << G4endl;

  //thickness for optical detector to be attached to fiber cores
  G4double opdet_dz = 0.3*mm; // thickness in z
  geo->GetOptD(fNam, "opdet_dz", opdet_dz, GeoParser::Unit(mm));

  //local materials for the cell
  BuildCellMaterials();

  //main volume for the cell
  G4Box *cell_shape = new G4Box(fNam+"_cell", cell_xy/2., cell_xy/2., cell_z/2.);
  G4Material *air_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, air_mat, cell_shape->GetName());

  //cell visibility
  ColorDecoder cell_vis("0:0:1:2");
  cell_vol->SetVisAttributes(cell_vis.MakeVis(geo, fNam, "cell_vis"));

  //absorber size along z, same as cell size, different variable for potential space at sides
  G4double abso_z = cell_z;

  //absorber material   FiberXYCal_WPowderplusEpoxy   G4_W
  G4Material *abso_mat = G4NistManager::Instance()->FindOrBuildMaterial("FiberXYCal_WPowderplusEpoxy");

  //absorber volume
  G4Box *abso_shape = new G4Box(fNam+"_abso", cell_xy/2., cell_xy/2., abso_z/2.);
  G4LogicalVolume *abso_vol = new G4LogicalVolume(abso_shape, abso_mat, abso_shape->GetName());

  //absorber visibility
  ColorDecoder abso_vis("0:0:1:3");
  abso_vol->SetVisAttributes(abso_vis.MakeVis(geo, fNam, "abso_vis"));

  //absorber in cell volume
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), abso_vol, abso_vol->GetName(), cell_vol, false, 0);

  //fiber dimensions
  G4double fiber_core_r = geo->GetD(fNam, "fiber_core_r")*mm; // core radius
  G4double fiber_clad_dr = 0.02*mm; // cladding thickness
  G4double fiber_air_dr = 0.01*mm; // air thickness in groove for the fiber
  geo->GetOptD(fNam, "fiber_clad_dr", fiber_clad_dr, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "fiber_air_dr", fiber_air_dr, GeoParser::Unit(mm));

  G4cout << "FiberXYCal::MakeCell, fiber_core_r, clad_dr, air_dr (mm): ";
  G4cout << fiber_core_r << ", " << fiber_clad_dr << ", " << fiber_air_dr << G4endl;

  //fiber groove filled with air
  G4double fiber_groove_r = fiber_core_r + fiber_clad_dr + fiber_air_dr;
  G4Tubs *groove_shape = new G4Tubs(fNam+"_groove", 0, fiber_groove_r, abso_z/2., 0, 360*deg);

  //groove volume
  G4Material *groove_mat = G4NistManager::Instance()->FindOrBuildMaterial("FiberXYCal_Air");
  G4LogicalVolume *groove_vol = new G4LogicalVolume(groove_shape, groove_mat, groove_shape->GetName());
  groove_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //fiber cladding
  G4Tubs *clad_shape = new G4Tubs(fNam+"_clad", 0, fiber_core_r+fiber_clad_dr, abso_z/2., 0, 360*deg);

  //cladding volume
  G4Material *clad_mat = G4NistManager::Instance()->FindOrBuildMaterial("FiberXYCal_Pethylene");
  G4LogicalVolume *clad_vol = new G4LogicalVolume(clad_shape, clad_mat, clad_shape->GetName());

  //cladding visibility
  ColorDecoder clad_vis("1:0:0:3");
  clad_vol->SetVisAttributes(clad_vis.MakeVis(geo, fNam, "clad_vis"));

  //cladding in the groove
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), clad_vol, clad_vol->GetName(), groove_vol, false, 0);

  //fiber core, sensitive volume, shorter to allow for optical detector
  G4Tubs *core_shape = new G4Tubs(fNam, 0, fiber_core_r, (abso_z/2.)-opdet_dz, 0, 360*deg);

  //core material
  G4Material *core_mat = G4NistManager::Instance()->FindOrBuildMaterial("FiberXYCal_Polystyrene");

  //core volume, sensitive
  G4LogicalVolume *core_vol = new G4LogicalVolume(core_shape, core_mat, fNam);

  //core visibility
  ColorDecoder core_vis("1:0:0:3"); // 0.4
  core_vol->SetVisAttributes(core_vis.MakeVis(geo, fNam, "core_vis"));

  //core in the cladding
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), core_vol, fNam, clad_vol, false, 0);

  //optical photon detector
  OpSiDet *opdet = new OpSiDet(fNam+"_opdet");
  fOpDet = dynamic_cast<Detector*>( opdet );

  //volume for optical detector
  ColorDecoder opdet_vis("1:1:0:3"); // 1
  G4LogicalVolume *opdet_vol = opdet->CreateGeometry(fiber_core_r, opdet_dz,
    opdet_vis.MakeVis(geo, fNam, "opdet_vis"));

  //put optical detector at the end of the fiber, negative local z
  new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*abso_z+0.5*opdet_dz), opdet_vol, opdet_vol->GetName(), clad_vol, false, 0);

  //mirror for reflection at opposite end to photon detector, positive local z
  G4Tubs *mirror_shape = new G4Tubs(fNam+"_mirror", 0, fiber_core_r, opdet_dz/2., 0, 360*deg);

  //mirror volume
  G4Material *mirror_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4LogicalVolume *mirror_vol = new G4LogicalVolume(mirror_shape, mirror_mat, mirror_shape->GetName());

  //mirror visibility
  ColorDecoder mirror_vis("0:1:0:3"); // 0.4
  mirror_vol->SetVisAttributes(mirror_vis.MakeVis(geo, fNam, "mirror_vis"));

  //put the mirror at the end of the fiber, positive local z
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*abso_z-0.5*opdet_dz), mirror_vol, mirror_vol->GetName(), clad_vol, false, 0);

  //mirror surface
  G4OpticalSurface *mirror_surf = new G4OpticalSurface("mirror_surf", glisur, ground, dielectric_metal, 1);

  G4MaterialPropertiesTable *mirror_tab = new G4MaterialPropertiesTable();
  std::vector<G4double> wlen = { 7, 0.1 }; // wavelength in micro meters
  G4OpticalMaterialProperties::ConvertToEnergy(wlen);
  std::vector<G4double> reflectivity = { 1, 1 };
  std::vector<G4double> efficiency = { 0, 0 };

  mirror_tab->AddProperty("REFLECTIVITY", wlen, reflectivity);
  mirror_tab->AddProperty("EFFICIENCY", wlen, efficiency);

  mirror_surf->SetMaterialPropertiesTable(mirror_tab);

  new G4LogicalSkinSurface("mirror_surf", mirror_vol, mirror_surf);

  //fiber spacing in module
  G4double fiber_dx = geo->GetD(fNam, "fiber_dx")*mm;

  G4cout << "FiberXYCal::MakeCell, fiber_dx (mm): " << fiber_dx << G4endl;

  G4int nfib = std::floor(cell_xy/fiber_dx); // number of fibers along x and y
  G4double ofs = (cell_xy - (nfib-1)*fiber_dx)/2.; // offset for the first fiber

  //fibers in the cell, by the groove, holding the clad fiber
  G4int fiber_cnt = 0; // fiber count in cell
  for(G4int ix=0; ix<nfib; ix++) {
    for(G4int iy=0; iy<nfib; iy++) {

      //fiber position in x and y
      G4double fib_x = -cell_xy/2. + ofs + ix*fiber_dx;
      G4double fib_y = -cell_xy/2. + ofs + iy*fiber_dx;

      //test for fiber at the edge
      if( fib_x-fiber_groove_r < -0.5*cell_xy or fib_x+fiber_groove_r > 0.5*cell_xy ) continue;
      if( fib_y-fiber_groove_r < -0.5*cell_xy or fib_y+fiber_groove_r > 0.5*cell_xy ) continue;

      //fiber (by the groove) in the absorber
      new G4PVPlacement(0, G4ThreeVector(fib_x, fib_y, 0), groove_vol, groove_vol->GetName(), abso_vol, false, fiber_cnt++);

    }//iy
  }//ix

  return cell_vol;

}//MakeCell

//_____________________________________________________________________________
void FiberXYCal::BuildCellMaterials() {

  //local materials for the cell

  G4NistManager *nist = G4NistManager::Instance();

  //polystyrene for fiber core
  if( !G4Material::GetMaterial("FiberXYCal_Polystyrene", false) ) {
    std::vector<G4String> elements = {"C", "H"};
    std::vector<G4int> natoms = {8, 8};
    G4Material *poly = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
    G4Material *mat = nist->ConstructNewMaterial("FiberXYCal_Polystyrene", elements, natoms, poly->GetDensity());

    //optical properties for polystyrene
    std::vector<G4double> rindex = { 1.6, 1.6 };
    std::vector<G4double> abslen = { 3*m, 3*m };
    std::vector<G4double> wlen = { 7, 0.1 }; // wavelength in micro meters
    G4OpticalMaterialProperties::ConvertToEnergy(wlen);

    //scintillation properties, preliminary according to examples/extended/optical/wls/
    std::vector<G4double> scin_wlen = { 0.46, 0.4 }; // micro meter
    G4OpticalMaterialProperties::ConvertToEnergy(scin_wlen);
    std::vector<G4double> scin_fast = { 1, 1 };

    G4MaterialPropertiesTable *tab = new G4MaterialPropertiesTable();
    tab->AddProperty("RINDEX", wlen, rindex);
    tab->AddProperty("ABSLENGTH", wlen, abslen);
    tab->AddProperty("SCINTILLATIONCOMPONENT1", scin_wlen, scin_fast);
    tab->AddConstProperty("SCINTILLATIONYIELD", 10. / keV);
    tab->AddConstProperty("RESOLUTIONSCALE", 1.0);
    tab->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 10. * ns);
    mat->SetMaterialPropertiesTable(tab);
  }

  //polyethylene for cladding
  if( !G4Material::GetMaterial("FiberXYCal_Pethylene", false) ) {
    std::vector<G4String> elements = {"C", "H"};
    std::vector<G4int> natoms = {1, 2};
    G4Material *pethy = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
    G4Material *mat = nist->ConstructNewMaterial("FiberXYCal_Pethylene", elements, natoms, pethy->GetDensity());

    //optical properties for polyethylene
    std::vector<G4double> rindex = { 1.49, 1.49 };
    std::vector<G4double> abslen = { 20*m, 20*m };
    std::vector<G4double> wlen = { 7, 0.1 }; // wavelength in micro meters
    G4OpticalMaterialProperties::ConvertToEnergy(wlen);

    G4MaterialPropertiesTable *tab = new G4MaterialPropertiesTable();
    tab->AddProperty("RINDEX", wlen, rindex);
    tab->AddProperty("ABSLENGTH", wlen, abslen);
    mat->SetMaterialPropertiesTable(tab);
  }

  //air for groove
  if( !G4Material::GetMaterial("FiberXYCal_Air", false) ) {
    G4Material *air = nist->FindOrBuildMaterial("G4_AIR");
    G4Material *mat = new G4Material("FiberXYCal_Air", air->GetDensity(), air);

    //optical properties for groove air
    G4MaterialPropertiesTable *tab = new G4MaterialPropertiesTable();
    tab->AddProperty("RINDEX", "Air");
    mat->SetMaterialPropertiesTable(tab);
  }

  //tungsten powder and epoxy for absorber
  if( !G4Material::GetMaterial("FiberXYCal_WPowderplusEpoxy", false) ) {
    std::vector<G4String> epoxy_elements = {"H", "C", "O"};
    std::vector<G4int> epoxy_natoms = {44, 15, 7};
    G4Material *epoxy = nist->ConstructNewMaterial("FiberXYCal_Epoxy", epoxy_elements, epoxy_natoms, 1.3*g/cm3);

    G4Material *tungsten = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");

    G4Material *mat = new G4Material("FiberXYCal_WPowderplusEpoxy", 10.95*g/cm3, 2);
    mat->AddMaterial(tungsten, 0.97);
    mat->AddMaterial(epoxy, 0.03);
  }

}//BuildCellMaterials

//_____________________________________________________________________________
G4bool FiberXYCal::ProcessHits(G4Step *step, G4TouchableHistory*) {

  //G4cout << "FiberXYCal::ProcessHits" << G4endl;

  //cell location in the module for the given fiber core
  const G4TouchableHandle& hnd = step->GetPreStepPoint()->GetTouchableHandle();

  //G4cout << hnd->GetVolume()->GetName() << G4endl;
  //G4cout << hnd->GetVolume()->GetName() << " " << hnd->GetVolume(1)->GetName();// << G4endl;
  //G4cout << " " << hnd->GetVolume(2)->GetName() << " " << hnd->GetVolume(3)->GetName();// << G4endl;
  //G4cout << " " << hnd->GetVolume(4)->GetName() << " " << hnd->GetVolume(5)->GetName() << G4endl;

  //volume hierarchy: core -> cladding -> groove -> absorber -> cell
  hnd->MoveUpHistory(4); // #4 to get from core to cell

  //global cell position
  G4ThreeVector origin(0, 0, 0);
  G4ThreeVector gpos = hnd->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);

  //ID for the cell
  Int_t cell_id = hnd->GetCopyNumber(0); // cell volume after MoveUpHistory(34)

  //G4cout << cell_id << " " << gpos.x() << " " << gpos.y() << " " << gpos.z() << G4endl;

  //make the hit for the cell
  CalPWOHits::Hit& hit = fHits.ConstructedAt(cell_id, CalPWOHits::Hit(cell_id, gpos.x()/mm, gpos.y()/mm, gpos.z()/mm));

  //increment energy deposit in the hit
  hit.en += step->GetTotalEnergyDeposit()/GeV;

  return true;

}//ProcessHits

//_____________________________________________________________________________
void FiberXYCal::Add(std::vector<Detector*> *vec) {

  //add this detector and its optical detector to sensitive detectors

  vec->push_back(this);
  fOpDet->Add(vec);

}//Add

//_____________________________________________________________________________
G4LogicalVolume* FiberXYCal::GetMotherVolume(G4String mother_nam, G4LogicalVolume *top) {

  for(size_t i=0; i<top->GetNoDaughters(); i++) {

    G4LogicalVolume *dv = top->GetDaughter(i)->GetLogicalVolume();

    if( dv->GetName() == mother_nam ) {
      return dv;
    }
  }

  return 0x0;

}//GetMotherVolume























