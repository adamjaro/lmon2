
//_____________________________________________________________________________
//
// Cherenkov fiber calorimeter with fibers along its axis (0 degrees),
// acording to ALICE ZDC in NIMA 456 (2001) 248-258
// and also NIMA 361 (1995) 161-179
//
// Geometry parameters:
// modx, mody, modz
// nx, ny
// cell_xy, cell_z
// cell_vis
// fiber_clad_D
// fiber_core_D
// fiber_dx
//
// optional parameters:
// xpos, ypos, zpos
// opdet_dz
// cell_mat_name
// core_vis
// mod_vis, opdet_vis
//
//_____________________________________________________________________________

//C++
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

//local classes
#include "GeoParser.h"
#include "ColorDecoder.h"
#include "OpSiDet.h"
#include "QCal1.h"

//_____________________________________________________________________________
QCal1::QCal1(const G4String& nam, GeoParser *geo, G4LogicalVolume *top): Detector(),
   G4VSensitiveDetector(nam), fNam(nam), fOpDet(0x0) {

  G4cout << "QCal1: " << fNam << G4endl;

  //module size, must accommodate for all cells and photon detector, explained below
  G4double modx = geo->GetD(fNam, "modx")*mm; // full size in x, mm
  G4double mody = geo->GetD(fNam, "mody")*mm; // full size in y, mm
  G4double modz = geo->GetD(fNam, "modz")*mm; // full size in z, mm

  //calorimeter module
  G4Box *mod_shape = new G4Box(fNam+"_mod", modx/2., mody/2., modz/2.);
  G4Material *mod_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *mod_vol = new G4LogicalVolume(mod_shape, mod_mat, mod_shape->GetName());

  //module visibility
  ColorDecoder mod_vis("1:0:0:3");
  mod_vol->SetVisAttributes(mod_vis.MakeVis(geo, fNam, "mod_vis"));

  //number of cells in module along x and y
  G4int nx = geo->GetI(fNam, "nx"); // num cells in x
  G4int ny = geo->GetI(fNam, "ny"); // num cells in y

  //cell volume
  G4LogicalVolume *cell_vol = MakeCell(geo);

  //cells in module, modx = nx*cell_xy must hold and the same for y
  G4int cell_cnt = 0; // cell count in module
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm; // cell transverse size
  for(G4int ix=0; ix<nx; ix++) {
    for(G4int iy=0; iy<ny; iy++) {

      //cell position in x and y
      G4double cell_posx = -modx/2. + cell_xy/2. + ix*cell_xy;
      G4double cell_posy = -mody/2. + cell_xy/2. + iy*cell_xy;

      //put the cell in the module
      new G4PVPlacement(0, G4ThreeVector(cell_posx, cell_posy, 0), cell_vol, cell_vol->GetName(), mod_vol, false, cell_cnt++);
    }//iy
  }//ix

  //module center position
  G4double xpos=0, ypos=0, zpos=0;
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "ypos", ypos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //module in top volume
  new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos), mod_vol, mod_vol->GetName(), top, false, 0);

}//QCal1

//_____________________________________________________________________________
G4LogicalVolume* QCal1::MakeCell(GeoParser *geo) {

  //cell (absorber) size
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm;
  G4double cell_z = geo->GetD(fNam, "modz")*mm; // initial value as modz
  geo->GetOptD(fNam, "cell_z", cell_z, GeoParser::Unit(mm));

  //main volume for the cell holding absorber with fibers and photon detector
  G4double modz = geo->GetD(fNam, "modz")*mm;
  G4Box *main_shape = new G4Box(fNam+"_main", cell_xy/2., cell_xy/2., modz/2.);
  G4Material *main_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *main_vol = new G4LogicalVolume(main_shape, main_mat, main_shape->GetName());
  main_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //local materials for the cell
  if( !G4Material::GetMaterial("QCal1_Copper", false) ) {
    new G4Material("QCal1_Copper", 29., 63.55*g/mole, 8.960*g/cm3); //z, a, density
  }

  //name for cell material
  G4String cell_mat_name = "G4_W";
  geo->GetOptS(fNam, "cell_mat_name", cell_mat_name);

  //select the cell material
  G4Material *cell_mat=0x0;
  if( cell_mat_name.find("QCal1_") != std::string::npos ) {
    cell_mat = G4Material::GetMaterial(cell_mat_name); // local material
  } else {
    cell_mat = G4NistManager::Instance()->FindOrBuildMaterial(cell_mat_name); // nist material
  }

  G4cout << "  " << fNam << ", cell material: " << cell_mat->GetName() << G4endl;

  //cell volume
  G4Box *cell_shape = new G4Box(fNam+"_cell", cell_xy/2., cell_xy/2., cell_z/2.);
  G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, cell_mat, cell_shape->GetName());

  //cell visibility
  ColorDecoder cell_vis("0:0:1:2");
  cell_vol->SetVisAttributes(cell_vis.MakeVis(geo, fNam, "cell_vis"));

  //fiber cladding
  G4double fiber_clad_D = geo->GetD(fNam, "fiber_clad_D")*mm; // cladding diameter
  G4Tubs *clad_shape = new G4Tubs(fNam+"_clad", 0, fiber_clad_D/2., cell_z/2., 0, 360*deg);

  //PMMA material
  G4Material *pmma_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS");

  //optics for PMMA
  G4MaterialPropertiesTable *pmma_tab = new G4MaterialPropertiesTable();
  pmma_tab->AddProperty("RINDEX", "PMMA");
  pmma_mat->SetMaterialPropertiesTable(pmma_tab);

  //fiber cladding volume
  G4LogicalVolume *clad_vol = new G4LogicalVolume(clad_shape, pmma_mat, clad_shape->GetName());
  clad_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //SiO2 material
  G4Material *siO2_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  //optics for SiO2 as Fused silica
  G4MaterialPropertiesTable *siO2_tab = new G4MaterialPropertiesTable();
  siO2_tab->AddProperty("RINDEX", "Fused Silica");
  siO2_mat->SetMaterialPropertiesTable(siO2_tab);

  //optical detector to be attached to fiber cores
  G4double opdet_dz = 0.3*mm; // thickness in z
  geo->GetOptD(fNam, "opdet_dz", opdet_dz, GeoParser::Unit(mm));

  //fiber core, sensitive volume
  G4double fiber_core_D = geo->GetD(fNam, "fiber_core_D")*mm; // core diameter
  //G4Tubs *core_shape = new G4Tubs(fNam, 0, fiber_core_D/2., cell_z/2., 0, 360*deg);
  G4Tubs *core_shape = new G4Tubs(fNam, 0, fiber_core_D/2., (cell_z-opdet_dz)/2., 0, 360*deg);

  //fiber core volume
  G4LogicalVolume *core_vol = new G4LogicalVolume(core_shape, siO2_mat, fNam);

  //core visibility
  ColorDecoder core_vis("1:0:0:3");
  core_vol->SetVisAttributes(core_vis.MakeVis(geo, fNam, "core_vis"));

  //core in the cladding
  //G4VPhysicalVolume *core_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), core_vol, fNam, clad_vol, false, 0);
  G4VPhysicalVolume *core_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*opdet_dz), core_vol, fNam, clad_vol, false, 0);

  //make the optical detector
  OpSiDet *opdet = new OpSiDet(fNam+"_opdet");
  fOpDet = dynamic_cast<Detector*>( opdet );
  //G4LogicalVolume *opdet_vol = opdet->CreateGeometry(cell_xy, cell_xy, opdet_dz,
  //  opdet_vis.MakeVis(geo, fNam, "opdet_vis"));
  ColorDecoder opdet_vis("1:1:0:3");
  G4LogicalVolume *opdet_vol = opdet->CreateGeometry(fiber_core_D/2, opdet_dz,
    opdet_vis.MakeVis(geo, fNam, "opdet_vis"));

  //optical detector at the end of the fiber
  G4VPhysicalVolume *opdet_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*cell_z+0.5*opdet_dz),
   opdet_vol, opdet_vol->GetName(), clad_vol, false, 0);

  //optical detector in the cell main volume
  //G4VPhysicalVolume *opdet_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*modz-cell_z-0.5*opdet_dz),
   //opdet_vol, opdet_vol->GetName(), main_vol, false, 0);

  //optical coupling from fiber core to optical detector
  //opdet->MakeBoundary(core_phys, opdet_phys);

  //fiber spacing in module
  G4double fiber_dx = geo->GetD(fNam, "fiber_dx")*mm;
  G4int nfib = std::floor(cell_xy/fiber_dx); // number of fibers along x and y
  G4double ofs = (cell_xy - (nfib-1)*fiber_dx)/2.; // offset for the first fiber

  //fibers in the cell, by the cladding holding the core
  G4int fiber_cnt = 0; // fiber count in cell
  for(G4int ix=0; ix<nfib; ix++) {
    for(G4int iy=0; iy<nfib; iy++) {

      //fiber position in x and y
      G4double fib_x = -cell_xy/2. + ofs + ix*fiber_dx;
      G4double fib_y = -cell_xy/2. + ofs + iy*fiber_dx;

      //test for fiber at the edge
      if( fib_x-0.5*fiber_clad_D < -0.5*cell_xy or fib_x+0.5*fiber_clad_D > 0.5*cell_xy ) continue;
      if( fib_y-0.5*fiber_clad_D < -0.5*cell_xy or fib_y+0.5*fiber_clad_D > 0.5*cell_xy ) continue;

      //fiber (by the cladding) in the cell
      new G4PVPlacement(0, G4ThreeVector(fib_x, fib_y, 0), clad_vol, clad_vol->GetName(), cell_vol, false, fiber_cnt++);

    }//iy
  }//ix

  //absorber cell in its main volume
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*(modz-cell_z)), cell_vol, cell_vol->GetName(), main_vol, false, 0);

  return main_vol;

}//MakeCell

//_____________________________________________________________________________
G4bool QCal1::ProcessHits(G4Step *step, G4TouchableHistory*) {

  //G4cout << "QCal1::ProcessHits" << G4endl;

  //cell location in the module
  const G4TouchableHandle& hnd = step->GetPreStepPoint()->GetTouchableHandle();

  //G4cout << hnd->GetVolume()->GetName() << " " << hnd->GetVolume(1)->GetName();// << G4endl;
  //G4cout << " " << hnd->GetVolume(2)->GetName() << " " << hnd->GetVolume(3)->GetName();// << G4endl;
  //G4cout << " " << hnd->GetVolume(4)->GetName() << " " << hnd->GetVolume(5)->GetName() << G4endl;

  //fiber_cladding (#0) -> cell (#1) -> cell_main (#2) -> module (#3) -> top (4)
  //G4cout << hnd->GetCopyNumber(0) << " " << hnd->GetCopyNumber(1);// << " " << hnd->GetCopyNumber(2);
  //G4cout << " " << hnd->GetCopyNumber(3) << " " << hnd->GetCopyNumber(4) << " " << hnd->GetCopyNumber(5) << G4endl;

  //get cell_main in the history: qcal qcal_clad qcal_cell qcal_main qcal_mod topv_p
  hnd->MoveUpHistory(3); // #3 from fiber core (qcal) to the cell main (qcal_main)

  //global cell position
  G4ThreeVector origin(0, 0, 0);
  G4ThreeVector gpos = hnd->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);

  //make the hit for the cell
  Int_t cell_id = hnd->GetCopyNumber(0); // cell_main after MoveUpHistory(3)
  CalPWOHits::Hit& hit = fHits.ConstructedAt(cell_id, CalPWOHits::Hit(cell_id, gpos.x()/mm, gpos.y()/mm, gpos.z()/mm));

  //deposited energy in step, GeV
  G4double edep = step->GetTotalEnergyDeposit()/GeV;

  //increment energy deposit in the hit
  hit.en += edep;

  return true;

}//ProcessHits

//_____________________________________________________________________________
void QCal1::CreateOutput(TTree *tree) {

  //create output for the hits
  fHits.CreateOutput(fNam, tree);

}//CreateOutput

//_____________________________________________________________________________
void QCal1::ClearEvent() {

  fHits.ClearEvent();

}//ClearEvent

//_____________________________________________________________________________
void QCal1::FinishEvent() {

  fHits.FinishEvent();

}//FinishEvent

//_____________________________________________________________________________
void QCal1::Add(std::vector<Detector*> *vec) {

  //add this detector and its optical detector to sensitive detectors

  vec->push_back(this);
  fOpDet->Add(vec);

}//Add


























