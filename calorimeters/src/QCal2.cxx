
//_____________________________________________________________________________
//
// Cherenkov fiber calorimeter with fibers at an angle to its axis
//
//_____________________________________________________________________________

//C++

//ROOT
#include "TTree.h"
#include "TMath.h"

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
#include "LogicVolFinder.h"
#include "OpSiDet.h"
#include "QCal2.h"

//_____________________________________________________________________________
QCal2::QCal2(const G4String& nam, GeoParser *geo, G4LogicalVolume *top): Detector(),
   G4VSensitiveDetector(nam), fNam(nam), fOpDet(NULL) {

  G4cout << "QCal2: " << fNam << G4endl;

  //module size
  G4double modx = geo->GetD(fNam, "modx")*mm; // full size in x, mm
  G4double mody = geo->GetD(fNam, "mody")*mm; // full size in y, mm
  G4double modz = geo->GetD(fNam, "modz")*mm; // full size in z, mm

  //calorimeter module
  G4Box *mod_shape = new G4Box(fNam+"_mod", modx/2., mody/2., modz/2.);
  G4Material *mod_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *mod_vol = new G4LogicalVolume(mod_shape, mod_mat, mod_shape->GetName());

  //module visibility
  ColorDecoder mod_vis("0:0:1:2");
  mod_vol->SetVisAttributes(mod_vis.MakeVis(geo, fNam, "mod_vis"));

  //cell volume
  G4LogicalVolume *cell_vol = MakeCell(geo);

  //cells in module
  G4RotationMatrix *cell_rot = new G4RotationMatrix(); // is HepRotation
  cell_rot->rotateX((90+geo->GetD(fNam, "cell_phi"))*deg);

  G4int cell_cnt = 0; // cell count in module

  //initial position for the first cell along z
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm;
  G4double cell_z = geo->GetD(fNam, "cell_z")*mm;
  G4double cell_phi = geo->GetD(fNam, "cell_phi")*deg; // the 'cell_phi' value is in rad after geo-> call is *deg
  G4double cell_posz_init = 0.5*(cell_xy*TMath::Cos(cell_phi) + cell_z*TMath::Sin(cell_phi));

  //number of cells in module along z and x
  G4int nz = geo->GetI(fNam, "nz"); // num cells in z
  G4int nx = geo->GetI(fNam, "nx"); // num cells in x

  G4double cell_posz = 0.5*modz - cell_posz_init;

  //z-loop
  for(G4int iz=0; iz<nz; iz++) {
    //x-loop
    for(G4int ix=0; ix<nx; ix++) {

      G4double cell_posx = -0.5*modx + 0.5*cell_xy + ix*cell_xy;
      new G4PVPlacement(cell_rot, G4ThreeVector(cell_posx, 0, cell_posz), cell_vol, cell_vol->GetName(), mod_vol, false, cell_cnt++);

    }//x-loop
    cell_posz -= cell_xy/TMath::Cos(cell_phi);

  }//z-loop




  //module center position
  G4double xpos=0, ypos=0, zpos=0;

  //module in its mother volume
  G4LogicalVolume *mother_vol = LogicVolFinder::GetMotherVolume("place_into", top, geo, fNam);
  //case of two levels of volume to place, volume 'place_into2' inside 'place_into1'
  if( geo->HasParameter(fNam, "place_into2") ) {
    mother_vol = LogicVolFinder::GetMotherVolume2("place_into1", "place_into2", top, geo, fNam);
  }

  new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos), mod_vol, mod_vol->GetName(), mother_vol, false, 0);

}//QCal2

//_____________________________________________________________________________
G4LogicalVolume* QCal2::MakeCell(GeoParser *geo) {

  //cell main volume size
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm;
  G4double cell_z = geo->GetD(fNam, "cell_z")*mm;

  //name for absorber material
  G4String cell_mat_name = "G4_W";
  geo->GetOptS(fNam, "cell_mat_name", cell_mat_name);

  //local materials for the absorber
  if( !G4Material::GetMaterial("QCal2_Copper", false) ) {
    new G4Material("QCal2_Copper", 29., 63.55*g/mole, 8.960*g/cm3); //z, a, density
  }

  //select the absorber material, local or NIST
  G4Material *cell_mat=0x0;
  if( cell_mat_name.find("QCal2_") != std::string::npos ) {
    cell_mat = G4Material::GetMaterial(cell_mat_name); // local material
  } else {
    cell_mat = G4NistManager::Instance()->FindOrBuildMaterial(cell_mat_name); // nist material
  }

  G4cout << "  " << fNam << ", absorber material: " << cell_mat->GetName() << G4endl;

  //cell main (absorber) volume
  G4Box *cell_shape = new G4Box(fNam+"_cell", cell_xy/2., cell_xy/2., cell_z/2.);
  G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, cell_mat, cell_shape->GetName());

  //cell visibility
  ColorDecoder cell_vis("1:0:0:2");
  cell_vol->SetVisAttributes(cell_vis.MakeVis(geo, fNam, "cell_vis"));

  //fiber cladding
  G4double fiber_clad_D = geo->GetD(fNam, "fiber_clad_D")*mm; // cladding diameter
  G4Tubs *clad_shape = new G4Tubs(fNam+"_clad", 0, fiber_clad_D/2., cell_z/2., 0, 360*deg);

  //PMMA material for the cladding
  G4Material *pmma_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS");

  //optics for PMMA
  G4MaterialPropertiesTable *pmma_tab = new G4MaterialPropertiesTable();
  pmma_tab->AddProperty("RINDEX", "PMMA");
  pmma_mat->SetMaterialPropertiesTable(pmma_tab);

  //fiber cladding volume
  G4LogicalVolume *clad_vol = new G4LogicalVolume(clad_shape, pmma_mat, clad_shape->GetName());
  ColorDecoder clad_vis("1:0:0:2");
  clad_vol->SetVisAttributes(clad_vis.MakeVis(geo, fNam, "clad_vis"));

  //fiber spacing in module
  G4double fiber_dx = geo->GetD(fNam, "fiber_dx")*mm;
  G4int nfib = std::floor(cell_xy/fiber_dx); // number of fibers along x and y
  G4double ofs = (cell_xy - (nfib-1)*fiber_dx)/2.; // offset for the first fiber

  //fibers (by the cladding) in the cell
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

  //optical detector to be attached to fiber cores
  OpSiDet *opdet = new OpSiDet(fNam+"_opdet");
  fOpDet = dynamic_cast<Detector*>( opdet );

  //dimensions for optical detector
  G4double fiber_core_D = geo->GetD(fNam, "fiber_core_D")*mm; // core diameter
  G4double opdet_dz = 1*mm; // thickness in z (0.3)
  geo->GetOptD(fNam, "opdet_dz", opdet_dz, GeoParser::Unit(mm));

  //volume for optical detector
  ColorDecoder opdet_vis("1:1:0:2");
  G4LogicalVolume *opdet_vol = opdet->CreateGeometry(fiber_core_D/2, opdet_dz, opdet_vis.MakeVis(geo, fNam, "opdet_vis"));

  //optical detector in the cladding, at the end of the fiber
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*cell_z-0.5*opdet_dz), opdet_vol, opdet_vol->GetName(), clad_vol, false, 0);

  //fiber core, sensitive volume
  G4Tubs *core_shape = new G4Tubs(fNam, 0, fiber_core_D/2., (cell_z-opdet_dz)/2., 0, 360*deg);

  //SiO2 material for the core
  G4Material *siO2_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  //optics for SiO2 as Fused silica
  G4MaterialPropertiesTable *siO2_tab = new G4MaterialPropertiesTable();
  siO2_tab->AddProperty("RINDEX", "Fused Silica");
  siO2_mat->SetMaterialPropertiesTable(siO2_tab);

  //fiber core volume
  G4LogicalVolume *core_vol = new G4LogicalVolume(core_shape, siO2_mat, fNam);
  ColorDecoder core_vis("0:1:1:1");
  core_vol->SetVisAttributes(core_vis.MakeVis(geo, fNam, "core_vis"));

  //core in the cladding
  G4VPhysicalVolume *core_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*opdet_dz), core_vol, fNam, clad_vol, false, 0);


  return cell_vol;

}//MakeCell

//_____________________________________________________________________________
void QCal2::Add(std::vector<Detector*> *vec) {

  //add this detector and its optical detector to sensitive detectors

  vec->push_back(this);
  fOpDet->Add(vec);

}//Add



















