
//_____________________________________________________________________________
//
// Cherenkov fiber calorimeter, cells at an angle to its axis, bent fibers
// to SiPM sensors
//
// Geometry parameters:
//
// modx, mody, modz
// xpos, ypos, zpos
// cell_xy, cell_z
// cell_phi
// nz, nx
// cell_mat_name
// fiber_clad_D
// fiber_core_D
// fiber_dx
// opdet_dz
// place_into, place_into1, place_into2
//
// mod_vis, cell_vis, clad_vis, core_vis
//_____________________________________________________________________________

//C++
#include <array>

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
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"

//local classes
#include "GeoParser.h"
#include "ColorDecoder.h"
#include "LogicVolFinder.h"
#include "OpSiDet.h"
#include "FiberYZ.h"
#include "QCal2Fibers.h"

//_____________________________________________________________________________
QCal2Fibers::QCal2Fibers(const G4String& nam, GeoParser *geo, G4LogicalVolume *top): Detector(),
   G4VSensitiveDetector(nam), fNam(nam), fOpDet(NULL) {

  G4cout << "QCal2Fibers: " << fNam << G4endl;

  //module size
  //G4double modx = geo->GetD(fNam, "modx")*mm; // full size in x, mm
  //G4double mody = geo->GetD(fNam, "mody")*mm; // full size in y, mm
  //G4double modz = geo->GetD(fNam, "modz")*mm; // full size in z, mm
  G4double modx = 35*mm; // full size in x, mm
  G4double mody = 35*mm; // full size in y, mm
  G4double modz = 40*mm; // full size in z, mm

  //calorimeter module
  G4Box *mod_shape = new G4Box(fNam+"_mod", modx/2., mody/2., modz/2.);
  G4Material *mod_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *mod_vol = new G4LogicalVolume(mod_shape, mod_mat, mod_shape->GetName());

  //module visibility
  ColorDecoder mod_vis("0:0:1:2");
  mod_vol->SetVisAttributes(mod_vis.MakeVis(geo, fNam, "mod_vis"));

  //cell volume
  G4LogicalVolume *cell_vol = MakeCell(geo);

  //FIXME single cell
  //G4RotationMatrix *cell_rot = new G4RotationMatrix(); // is HepRotation
  //cell_rot->rotateX((90+geo->GetD(fNam, "cell_phi"))*deg);
  //cell_rot->rotateX(10*deg);
  //new G4PVPlacement(cell_rot, G4ThreeVector(0, 0, 0), cell_vol, cell_vol->GetName(), mod_vol, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), cell_vol, cell_vol->GetName(), mod_vol, false, 0);

/*

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

  */

  //module center position
  G4double xpos=0, ypos=0, zpos=0;
  geo->GetOptD(fNam, "xpos", xpos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "ypos", ypos, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "zpos", zpos, GeoParser::Unit(mm));

  //module in its mother volume
  G4LogicalVolume *mother_vol = LogicVolFinder::GetMotherVolume("place_into", top, geo, fNam);
  //case of two levels of volume to place, volume 'place_into2' inside 'place_into1'
  if( geo->HasParameter(fNam, "place_into2") ) {
    mother_vol = LogicVolFinder::GetMotherVolume2("place_into1", "place_into2", top, geo, fNam);
  }

  new G4PVPlacement(0, G4ThreeVector(xpos, ypos, zpos), mod_vol, mod_vol->GetName(), mother_vol, false, 0);

}//QCal2Fibers

//_____________________________________________________________________________
G4LogicalVolume* QCal2Fibers::MakeCell(GeoParser *geo) {

  //cell main volume size
  G4double cell_xy = geo->GetD(fNam, "cell_xy")*mm;
  //G4double cell_z = geo->GetD(fNam, "cell_z")*mm;
  //G4double cell_xy = 32*mm;
  G4double cell_z = 37*mm;

  //cell main  volume
  //G4Box *cell_shape = new G4Box(fNam+"_cell", cell_xy/2., cell_xy/2., cell_z/2.); // FIXME fibers to be sensitive
  G4Box *cell_shape = new G4Box(fNam, cell_xy/2., cell_xy/2., cell_z/2.);
  G4Material *cell_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, cell_mat, cell_shape->GetName());
  //G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, pmma_mat, cell_shape->GetName());

  //cell visibility
  ColorDecoder cell_vis("1:0:0:2");
  cell_vol->SetVisAttributes(cell_vis.MakeVis(geo, fNam, "cell_vis"));

  //length in z for individual sections
  G4double opdet_dz = 1*mm; // thickness in z for optical detector
  geo->GetOptD(fNam, "opdet_dz", opdet_dz, GeoParser::Unit(mm));

  //length for lightguide section
  G4double lguide_z = geo->GetD(fNam, "lguide_z")*mm;

  //length for in z for absorber and Cherenkov fibers
  G4double abso_z = cell_z - lguide_z - opdet_dz;

  //diameter for Cherenkov fiber cladding and core
  G4double fiber_clad_D = geo->GetD(fNam, "fiber_clad_D")*mm; // cladding diameter
  G4double fiber_core_D = geo->GetD(fNam, "fiber_core_D")*mm; // core diameter

  //PMMA material for fiber cladding
  G4Material *pmma_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS");

  //optics for PMMA
  G4MaterialPropertiesTable *pmma_tab = new G4MaterialPropertiesTable();
  pmma_tab->AddProperty("RINDEX", "PMMA");
  pmma_mat->SetMaterialPropertiesTable(pmma_tab);

  //SiO2 material for fiber core
  G4Material *siO2_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  //optics for SiO2 as Fused silica
  G4MaterialPropertiesTable *siO2_tab = new G4MaterialPropertiesTable();
  siO2_tab->AddProperty("RINDEX", "Fused Silica");
  siO2_mat->SetMaterialPropertiesTable(siO2_tab);

  //Cherenkov radiator fiber in cladding, FIXME fiber core to be sensitive volume
  G4LogicalVolume *clad_vol = MakeStraightFib(fiber_clad_D, fiber_core_D, abso_z, pmma_mat, siO2_mat,
    fNam+"_clad", fNam+"_core", geo);


  //fiber (by cladding) in the cell
  new G4PVPlacement(0, G4ThreeVector(-6,0,0.5*abso_z-0.5*cell_z), clad_vol, clad_vol->GetName(), cell_vol, false, 0);
  new G4PVPlacement(0, G4ThreeVector(6,0,0.5*abso_z-0.5*cell_z), clad_vol, clad_vol->GetName(), cell_vol, false, 1);




  //optical detector at positive z end
  OpSiDet *opdet = new OpSiDet(fNam+"_opdet");
  fOpDet = dynamic_cast<Detector*>( opdet );

  //placement for the optical detector
  ColorDecoder opdet_vis("1:1:0:1");
  G4LogicalVolume *opdet_vol = opdet->CreateGeometry(6*mm, 6*mm, opdet_dz, opdet_vis.MakeVis(geo, fNam, "opdet_vis"));
  new G4PVPlacement(0, G4ThreeVector(0,0,0.5*cell_z-0.5*opdet_dz), opdet_vol, opdet_vol->GetName(), cell_vol, false, 0);

  //straight section at the end of bent fibers
  G4double fib_lguide_end_dz = 5*mm;
  geo->GetOptD(fNam, "fib_lguide_end_dz", fib_lguide_end_dz, GeoParser::Unit(mm));

  G4LogicalVolume *lguide_end_vol = MakeStraightFib(fiber_clad_D, fiber_core_D, fib_lguide_end_dz, pmma_mat, siO2_mat,
    fNam+"_lguide_end_clad", fNam+"_lguide_end_core", geo);
  new G4PVPlacement(0, G4ThreeVector(-2,0,0.5*cell_z-opdet_dz-0.5*fib_lguide_end_dz), lguide_end_vol,
    lguide_end_vol->GetName(), cell_vol, false, 0); // negative x
  new G4PVPlacement(0, G4ThreeVector(2,0,0.5*cell_z-opdet_dz-0.5*fib_lguide_end_dz), lguide_end_vol,
    lguide_end_vol->GetName(), cell_vol, false, 1); // positive x

  //bent fiber at negative x
  G4LogicalVolume *fibYZ_neg_x_clad = MakeFiberYZ(lguide_z-fib_lguide_end_dz, 4, fiber_clad_D/2,
    "fibYZ_neg_x_clad", pmma_mat, TMath::Pi()/2);
  G4LogicalVolume *fibYZ_neg_x_core = MakeFiberYZ(lguide_z-fib_lguide_end_dz, 4, fiber_core_D/2,
    "fibYZ_neg_x_core", siO2_mat, TMath::Pi()/2);
  new G4PVPlacement(0, G4ThreeVector(0,0,0), fibYZ_neg_x_core, fibYZ_neg_x_core->GetName(), fibYZ_neg_x_clad, false, 0);

  new G4PVPlacement(0, G4ThreeVector(-2,0,-0.5*cell_z+abso_z), fibYZ_neg_x_clad, fibYZ_neg_x_clad->GetName(), cell_vol, false, 0);

  //bent fiber at positive x
  G4LogicalVolume *fibYZ_pos_x_clad = MakeFiberYZ(lguide_z-fib_lguide_end_dz, 4, fiber_clad_D/2,
    "fibYZ_pos_x_clad", pmma_mat, -TMath::Pi()/2);
  G4LogicalVolume *fibYZ_pos_x_core = MakeFiberYZ(lguide_z-fib_lguide_end_dz, 4, fiber_core_D/2,
    "fibYZ_pos_x_core", siO2_mat, -TMath::Pi()/2);
  new G4PVPlacement(0, G4ThreeVector(0,0,0), fibYZ_pos_x_core, fibYZ_pos_x_core->GetName(), fibYZ_pos_x_clad, false, 0);

  new G4PVPlacement(0, G4ThreeVector(2,0,-0.5*cell_z+abso_z), fibYZ_pos_x_clad, fibYZ_pos_x_clad->GetName(), cell_vol, false, 0);

  //visibility for bent fibers
  ColorDecoder fibYZ_clad_vis("0:1:1:2");
  fibYZ_neg_x_clad->SetVisAttributes(fibYZ_clad_vis.MakeVis(geo, fNam, "fibYZ_clad_vis"));
  fibYZ_pos_x_clad->SetVisAttributes(fibYZ_clad_vis.MakeVis(geo, fNam, "fibYZ_clad_vis"));
  ColorDecoder fibYZ_vis("0:1:1:0.3");
  fibYZ_neg_x_core->SetVisAttributes(fibYZ_vis.MakeVis(geo, fNam, "fibYZ_vis"));
  fibYZ_pos_x_core->SetVisAttributes(fibYZ_vis.MakeVis(geo, fNam, "fibYZ_vis"));



/*
  G4double ypos = 4;


  //bent fiber, cladding
  G4LogicalVolume *fiberYZ_clad_vol = MakeFiberYZ(10, ypos, 1.4, "fiberYZ_clad", pmma_mat);
  ColorDecoder fibYZ_clad_vis("0:1:1:2");
  fiberYZ_clad_vol->SetVisAttributes(fibYZ_clad_vis.MakeVis(geo, fNam, "fibYZ_clad_vis"));

  new G4PVPlacement(0, G4ThreeVector(0,0,0), fiberYZ_clad_vol, fiberYZ_clad_vol->GetName(), cell_vol, false, 0);

  //bent fiber, core
  G4LogicalVolume *fiberYZ_vol = MakeFiberYZ(10, ypos, 1.1, "fiberYZ", siO2_mat);
  ColorDecoder fibYZ_vis("0:1:1:0.3");
  fiberYZ_vol->SetVisAttributes(fibYZ_vis.MakeVis(geo, fNam, "fibYZ_vis"));
  new G4PVPlacement(0, G4ThreeVector(0,0,0), fiberYZ_vol, fiberYZ_vol->GetName(), fiberYZ_clad_vol, false, 0);


  //fiber cladding
  G4Tubs *clad_shape = new G4Tubs(fNam+"_clad", 0, 1.4, 4/2, 0, 360*deg);

  //fiber cladding volume
  G4LogicalVolume *clad_vol = new G4LogicalVolume(clad_shape, pmma_mat, clad_shape->GetName());
  ColorDecoder clad_vis("1:0:0:2");
  clad_vol->SetVisAttributes(clad_vis.MakeVis(geo, fNam, "clad_vis"));

  //upper fiber by cladding
  //new G4PVPlacement(0, G4ThreeVector(0,ypos,-2), clad_vol, clad_vol->GetName(), cell_vol, false, 0);
  new G4PVPlacement(0, G4ThreeVector(ypos,0,-2), clad_vol, clad_vol->GetName(), cell_vol, false, 0);


  //fiber core, FIXME to be sensitive volume
  G4Tubs *core_shape = new G4Tubs(fNam+"_core", 0, 1.1, 4/2, 0, 360*deg);
  G4LogicalVolume *core_vol = new G4LogicalVolume(core_shape, siO2_mat, core_shape->GetName());
  ColorDecoder core_vis("0:1:1:0.3");
  core_vol->SetVisAttributes(core_vis.MakeVis(geo, fNam, "core_vis"));

  //fiber in cladding
  new G4PVPlacement(0, G4ThreeVector(0,0.0), core_vol, core_vol->GetName(), clad_vol, false, 0);



  //new G4PVPlacement(0, G4ThreeVector(0,8,12), core_vol, core_vol->GetName(), cell_vol, false, 0);
  //new G4PVPlacement(0, G4ThreeVector(0,8,8), core_vol, core_vol->GetName(), cell_vol, false, 1);

  //lower fiber
  //new G4PVPlacement(0, G4ThreeVector(0,0,-2), core_vol, core_vol->GetName(), cell_vol, false, 0);

  //optical detector to be attached to fiber cores
  OpSiDet *opdet = new OpSiDet(fNam+"_opdet");
  fOpDet = dynamic_cast<Detector*>( opdet );

  ColorDecoder opdet_vis("1:1:0:1");
  //G4LogicalVolume *opdet_vol = opdet->CreateGeometry(1.1, 4, opdet_vis.MakeVis(geo, fNam, "opdet_vis"));
  //new G4PVPlacement(0, G4ThreeVector(0,8,7.8), opdet_vol, opdet_vol->GetName(), cell_vol, false, 0);
  //new G4PVPlacement(0, G4ThreeVector(0,0,-2), opdet_vol, opdet_vol->GetName(), cell_vol, false, 0);
  //new G4PVPlacement(0, G4ThreeVector(0,0,-6), opdet_vol, opdet_vol->GetName(), cell_vol, false, 0);

  G4LogicalVolume *opdet_vol = opdet->CreateGeometry(4, 4, 1, opdet_vis.MakeVis(geo, fNam, "opdet_vis"));
  new G4PVPlacement(0, G4ThreeVector(0,0,10.5), opdet_vol, opdet_vol->GetName(), cell_vol, false, 0);
*/


  return cell_vol;

}//MakeCell

//_____________________________________________________________________________
G4LogicalVolume* QCal2Fibers::MakeFiberYZ(Double_t L, Double_t yL, Double_t r, const G4String& nam, G4Material *mat, Double_t theta) {

  FiberYZ fib(L, yL, r, 0.1); // 0.9  0.1  0.05
  fib.InvertZ();
  fib.RotateXY(theta);

  //G4cout << "fib: " << fib.GetNFacets() << G4endl;

  G4TessellatedSolid *fiberYZ_shape = new G4TessellatedSolid(nam);

  //facet loop
  for(size_t i=0; i<fib.GetNFacets(); i++) {

    //get the facet
    const FiberYZ::facet& fct = fib.GetFacet(i);

    //facet points
    //const std::array<Double_t, 3>& p0 = fct.p0;
    //const std::array<Double_t, 3>& p1 = fct.p1;
    //const std::array<Double_t, 3>& p2 = fct.p2;

    //reversed order with inverted z
    const std::array<Double_t, 3>& p0 = fct.p2;
    const std::array<Double_t, 3>& p1 = fct.p1;
    const std::array<Double_t, 3>& p2 = fct.p0;

    G4TriangularFacet *gf = new
    G4TriangularFacet(G4ThreeVector(p0[0],p0[1],p0[2]), G4ThreeVector(p1[0],p1[1],p1[2]), G4ThreeVector(p2[0],p2[1],p2[2]), ABSOLUTE);

    fiberYZ_shape->AddFacet(gf);

  }//facet loop

  fiberYZ_shape->SetSolidClosed(true);

  G4LogicalVolume *fiberYZ_vol = new G4LogicalVolume(fiberYZ_shape, mat, fiberYZ_shape->GetName());

  return fiberYZ_vol;

}//MakeFiberYZ

//_____________________________________________________________________________
G4LogicalVolume* QCal2Fibers::MakeStraightFib(G4double cladD, G4double coreD, G4double Lz,
  G4Material *clad_mat, G4Material *core_mat, const G4String& clad_nam, const G4String& core_nam, GeoParser *geo) {

  //straight fiber

  //fiber cladding
  G4Tubs *clad_shape = new G4Tubs(clad_nam, 0, cladD/2., Lz/2., 0, 360*deg);

  G4LogicalVolume *clad_vol = new G4LogicalVolume(clad_shape, clad_mat, clad_shape->GetName());
  ColorDecoder clad_vis("1:0:0:2");
  clad_vol->SetVisAttributes(clad_vis.MakeVis(geo, fNam, "clad_vis"));

  G4Tubs *core_shape = new G4Tubs(core_nam, 0, coreD/2., Lz/2., 0, 360*deg);

  G4LogicalVolume *core_vol = new G4LogicalVolume(core_shape, core_mat, core_shape->GetName());
  ColorDecoder core_vis("0:1:1:0.3");
  core_vol->SetVisAttributes(core_vis.MakeVis(geo, fNam, "core_vis"));

  //fiber core in cladding
  new G4PVPlacement(0, G4ThreeVector(0,0,0), core_vol, core_vol->GetName(), clad_vol, false, 0);

  return clad_vol;

}//MakeStraightFib

//_____________________________________________________________________________
void QCal2Fibers::Add(std::vector<Detector*> *vec) {

  //add this detector and its optical detector to sensitive detectors

  vec->push_back(this);
  fOpDet->Add(vec);

}//Add

//_____________________________________________________________________________
G4bool QCal2Fibers::ProcessHits(G4Step *step, G4TouchableHistory*) {

  return true;

}//ProcessHits













