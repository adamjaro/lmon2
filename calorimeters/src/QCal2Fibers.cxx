
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
// mod_vis, cell_vis, clad_vis, core_vis, fibYZ_vis
//_____________________________________________________________________________

//C++
#include <array>
#include <vector>

//ROOT
#include "TTree.h"
#include "TMath.h"

//Geant
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4MultiUnion.hh"
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
  G4Box *cell_shape = new G4Box(fNam+"_cell", cell_xy/2., cell_xy/2., cell_z/2.);
  G4Material *cell_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, cell_mat, cell_shape->GetName());
  //G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, pmma_mat, cell_shape->GetName());

  //cell visibility
  ColorDecoder cell_vis("1:0:0:2");
  cell_vol->SetVisAttributes(cell_vis.MakeVis(geo, fNam, "cell_vis"));

  //optical detector at positive z end
  OpSiDet *opdet = new OpSiDet(fNam+"_opdet");
  opdet->SetUpHistory(1); // from opdet to cell in its hits
  fOpDet = dynamic_cast<Detector*>( opdet );

  //placement for the optical detector
  G4double opdet_dz=1*mm, opdet_dxy=6*mm; // thickness in z and size in xy for optical detector
  geo->GetOptD(fNam, "opdet_dz", opdet_dz, GeoParser::Unit(mm));
  geo->GetOptD(fNam, "opdet_dxy", opdet_dxy, GeoParser::Unit(mm));
  ColorDecoder opdet_vis("1:1:0:1");
  G4LogicalVolume *opdet_vol = opdet->CreateGeometry(opdet_dxy, opdet_dxy, opdet_dz, opdet_vis.MakeVis(geo, fNam, "opdet_vis"));
  new G4PVPlacement(0, G4ThreeVector(0,0,0.5*cell_z-0.5*opdet_dz), opdet_vol, opdet_vol->GetName(), cell_vol, false, 0);

  //length for lightguide section
  G4double lguide_z = geo->GetD(fNam, "lguide_z")*mm;

  //length in z for absorber and Cherenkov fibers
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

  //Cherenkov radiator fiber in cladding, same lenght as absorber, fiber core as sensitive volume
  G4LogicalVolume *clad_vol = MakeStraightFib(fiber_clad_D, fiber_core_D, abso_z, pmma_mat, siO2_mat,
    fNam+"_clad", fNam, geo);

  //fiber spacing in cell
  const G4double fiber_dx = geo->GetD(fNam, "fiber_dx")*mm;

  //number of fibers along x and y in the cell
  const G4int nfib = std::floor((cell_xy-fiber_dx)/fiber_dx)+1;

  //offset for the first fiber (number of intervals between fibers is nfib-1)
  const G4double ofs = 0.5*(cell_xy - (nfib-1)*fiber_dx);

  //spacing for bent lightguide fibers at optical detector
  G4double odet_fib_dxy = opdet_dxy/nfib;

  //fiber positions for lightguide
  std::vector<std::array<G4double, 4>> fib_lguide_pos;

  //fibers (by the cladding) in the cell
  G4int fiber_cnt = 0; // fiber count in cell
  for(G4int ix=0; ix<nfib; ix++) {
    for(G4int iy=0; iy<nfib; iy++) {

      //fiber position in x and y
      G4double fib_x = -0.5*cell_xy + ofs + ix*fiber_dx;
      G4double fib_y = -0.5*cell_xy + ofs + iy*fiber_dx;

      //fiber (by the cladding) in the cell, starting from negative end in z
      new G4PVPlacement(0, G4ThreeVector(fib_x, fib_y, -0.5*cell_z+0.5*abso_z), clad_vol, clad_vol->GetName(),
        cell_vol, false, fiber_cnt++);

      //corresponding lightguide fiber at SiPM
      G4double odet_fib_x = -0.5*opdet_dxy + 0.5*odet_fib_dxy + ix*odet_fib_dxy;
      G4double odet_fib_y = -0.5*opdet_dxy + 0.5*odet_fib_dxy + iy*odet_fib_dxy;

      fib_lguide_pos.push_back({fib_x, fib_y, odet_fib_x, odet_fib_y}); // indices: 0, 1, 2, 3

      //G4cout << fib_x << " " << fib_y << " " << odet_fib_x << " " << odet_fib_y << G4endl; 

    }//iy
  }//ix

  //lightguide straight section at the end of bent fibers
  G4double fib_lguide_end_dz = 5*mm;
  geo->GetOptD(fNam, "fib_lguide_end_dz", fib_lguide_end_dz, GeoParser::Unit(mm));

  //volume for straight end fibers
  G4LogicalVolume *lguide_end_vol = MakeStraightFib(fiber_clad_D, fiber_core_D, fib_lguide_end_dz, pmma_mat, siO2_mat,
    fNam+"_lguide_end_clad", fNam+"_lguide_end_core", geo);

  //placement for straight end fibers
  G4int fib_end_cnt=0;
  for(std::array<G4double, 4>& pos: fib_lguide_pos) {

    //x and y at the sensor at indices 2 and 3
    new G4PVPlacement(0, G4ThreeVector(pos[2],pos[3],0.5*cell_z-opdet_dz-0.5*fib_lguide_end_dz), lguide_end_vol,
      lguide_end_vol->GetName(), cell_vol, false, fib_end_cnt++);

  }//lguide pos loop

  //lightguide by bent fibers
  fib_end_cnt = 0;
  for(std::array<G4double, 4>& pos: fib_lguide_pos) {

    //distance in x and y to connect by the fiber
    G4double dx = pos[0]-pos[2]; // fib_x - odet_fib_x
    G4double dy = pos[1]-pos[3]; // fib_y - odet_fib_y
    G4double Lz = TMath::Sqrt(dx*dx+dy*dy);

    //angle for the fiber
    G4double theta = -0.5*TMath::Pi(); // bent fiber natively along y axis

    //positive y
    if(pos[1] > 0 and Lz > 1e-9) {

      theta += TMath::ACos(dx/Lz);

      //G4cout << dx << " " << dy << " " << Lz << " " << TMath::ACos(dx/Lz) << G4endl;
    } else if(Lz > 1e-9) {

      //negative y
      theta += TMath::Pi()+TMath::ACos(-dx/Lz);

      //if(pos[1] < 0 and Lz > 1e-9) {
      //G4cout << dx << " " << dy << " " << Lz << " " << TMath::Pi()+TMath::ACos(-dx/Lz) << G4endl;
    }

    //construct the fiber
    G4LogicalVolume *fibYZ_clad = MakeFiberYZ(lguide_z-fib_lguide_end_dz, Lz, fiber_clad_D/2,
      "fibYZ_clad_"+std::to_string(fib_end_cnt), pmma_mat, theta);
    G4LogicalVolume *fibYZ_core = MakeFiberYZ(lguide_z-fib_lguide_end_dz, Lz, fiber_core_D/2,
      "fibYZ_core_"+std::to_string(fib_end_cnt), siO2_mat, theta);

    //fiber visibility
    ColorDecoder fibYZ_clad_vis("0:1:1:3");
    fibYZ_clad->SetVisAttributes(fibYZ_clad_vis.MakeVis(geo, fNam, "fibYZ_clad_vis"));
    ColorDecoder fibYZ_vis("0:1:1:0.3");
    fibYZ_core->SetVisAttributes(fibYZ_vis.MakeVis(geo, fNam, "fibYZ_vis"));

    //core in the cladding
    new G4PVPlacement(0, G4ThreeVector(0,0,0), fibYZ_core, fibYZ_core->GetName(), fibYZ_clad, false, 0);

    //bent fiber in the cell
    new G4PVPlacement(0, G4ThreeVector(pos[2],pos[3],-0.5*cell_z+abso_z), fibYZ_clad, fibYZ_clad->GetName(), cell_vol, false, 0);

    fib_end_cnt++;

  }//lguide pos loop

  //name for absorber material
  G4String abso_mat_name = "G4_W";
  geo->GetOptS(fNam, "abso_mat_name", abso_mat_name);

  //local materials for the absorber
  if( !G4Material::GetMaterial("QCal2Fibers_Copper", false) ) {
    new G4Material("QCal2Fibers_Copper", 29., 63.55*g/mole, 8.960*g/cm3); //z, a, density
  }

  //select the absorber material, local or NIST
  G4Material *abso_mat=0x0;
  if( abso_mat_name.find("QCal2Fibers_") != std::string::npos ) {
    abso_mat = G4Material::GetMaterial(abso_mat_name); // local material
  } else {
    abso_mat = G4NistManager::Instance()->FindOrBuildMaterial(abso_mat_name); // nist material
  }

  G4cout << "  " << fNam << ", absorber material: " << abso_mat->GetName() << G4endl;

  //absorber opening for fibers
  G4double fiber_abso_delt = 0.02*mm; // add space between fibers and absorber
  geo->GetOptD(fNam, "fiber_abso_delt", fiber_abso_delt, GeoParser::Unit(mm));
  G4double fiber_abso_xy = fiber_clad_D+fiber_abso_delt;

  //absorber block
  G4MultiUnion *abso_shape = new G4MultiUnion(fNam+"_abso_shape");

  //upper and lower (in y) absorber base plates
  G4Box *abso_base_shape = new G4Box(fNam+"_abso_base_shape", cell_xy/2., (ofs-0.5*fiber_abso_xy)/2., abso_z/2.);
  abso_shape->AddNode(abso_base_shape, G4Transform3D(G4RotationMatrix(),
    G4ThreeVector(0, -0.5*cell_xy+abso_base_shape->GetYHalfLength(), 0))); // lower
  abso_shape->AddNode(abso_base_shape, G4Transform3D(G4RotationMatrix(),
    G4ThreeVector(0, 0.5*cell_xy-abso_base_shape->GetYHalfLength(), 0))); // upper

  //plates between layers of fibers (in y)
  G4Box *abso_layer_shape = new G4Box(fNam+"_abso_layer_shape", cell_xy/2., (fiber_dx-fiber_abso_xy)/2., abso_z/2.);
  for(G4int iy=0; iy<nfib-1; iy++) {
    G4double posY = -0.5*cell_xy + ofs + 0.5*fiber_dx + iy*fiber_dx; // position in y between the fibers
    abso_shape->AddNode(abso_layer_shape, G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, posY, 0)));
  }//iy

  //edge part in layer of fibers
  G4Box *abso_edge_shape = new G4Box(fNam+"_abso_edge_shape", (ofs-0.5*fiber_abso_xy)/2., fiber_abso_xy/2., abso_z/2.);
  for(G4int iy=0; iy<nfib; iy++) {
    G4double posY = -0.5*cell_xy + ofs + iy*fiber_dx; // fiber position in y
    G4double posX1 = 0.5*cell_xy - abso_edge_shape->GetXHalfLength(); // positive x
    G4double posX2 = -0.5*cell_xy + abso_edge_shape->GetXHalfLength(); // negative x

    abso_shape->AddNode(abso_edge_shape, G4Transform3D(G4RotationMatrix(), G4ThreeVector(posX1, posY, 0)));
    abso_shape->AddNode(abso_edge_shape, G4Transform3D(G4RotationMatrix(), G4ThreeVector(posX2, posY, 0)));
  }//iy

  //close the absorber block
  abso_shape->Voxelize();

  //absorber volume
  G4LogicalVolume *abso_vol = new G4LogicalVolume(abso_shape, abso_mat, abso_shape->GetName());
  ColorDecoder abso_vis("1:0:0:2");
  abso_vol->SetVisAttributes(abso_vis.MakeVis(geo, fNam, "abso_vis"));

  //absorber in the cell
  new G4PVPlacement(0, G4ThreeVector(0,0,-0.5*cell_z+0.5*abso_z), abso_vol, abso_vol->GetName(), cell_vol, false, 0);

  //absorber middle parts in layer of fibers (outside of G4MultiUnion because visualization would not work)
  G4Box *abso_mid_shape = new G4Box(fNam+"_abso_mid_shape", (fiber_dx-fiber_abso_xy)/2., fiber_abso_xy/2., abso_z/2.);
  G4LogicalVolume *abso_mid_vol = new G4LogicalVolume(abso_mid_shape, abso_mat, abso_mid_shape->GetName());
  abso_mid_vol->SetVisAttributes(abso_vis.MakeVis(geo, fNam, "abso_vis"));

  G4int abso_mid_shape_cnt = 0;
  for(G4int iy=0; iy<nfib; iy++) {
    G4double posY = -0.5*cell_xy + ofs + iy*fiber_dx; // fiber position in y

    for(G4int ix=0; ix<nfib-1; ix++) {
      G4double posX = -0.5*cell_xy + ofs + 0.5*fiber_dx + ix*fiber_dx; // position in x between the fibers

      new G4PVPlacement(0, G4ThreeVector(posX,posY,-0.5*cell_z+0.5*abso_z),
        abso_mid_vol, abso_mid_vol->GetName(), cell_vol, false, abso_mid_shape_cnt++);
    }//ix
  }//iy









/*
  //fiber (by cladding) in the cell
  new G4PVPlacement(0, G4ThreeVector(-6,0,0.5*abso_z-0.5*cell_z), clad_vol, clad_vol->GetName(), cell_vol, false, 0);
  new G4PVPlacement(0, G4ThreeVector(6,0,0.5*abso_z-0.5*cell_z), clad_vol, clad_vol->GetName(), cell_vol, false, 1);

  //optical detector at positive z end
  OpSiDet *opdet = new OpSiDet(fNam+"_opdet");
  opdet->SetUpHistory(1); // from opdet to cell
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
  ColorDecoder clad_vis("1:0:0:3");
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













