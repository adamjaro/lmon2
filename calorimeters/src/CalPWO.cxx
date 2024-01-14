

//_____________________________________________________________________________
//
// PbWO4 calorimeter, module consists of individual crystal cells
//
//_____________________________________________________________________________

//C++
#include <vector>

//ROOT
#include "TTree.h"

//Geant
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4RunManager.hh"

//local classes
#include "CalPWO.h"
#include "PhotoCathPMT.h"
#include "GeoParser.h"
#include "ColorDecoder.h"
#include "MCParticleAction.h"

//_____________________________________________________________________________
CalPWO::CalPWO(const G4String& nam, GeoParser *geo, G4LogicalVolume *top) : Detector(),
  G4VSensitiveDetector(nam), fNam(nam), fStack(0) {

  G4cout << "CalPWO: " << fNam << G4endl;

  //module size, must accommodate for wrapping, crystal, photon detector and number of cells, explained below
  G4double modx = geo->GetD(fNam, "modx")*mm; // full size in x, mm
  G4double mody = geo->GetD(fNam, "mody")*mm; // full size in y, mm
  G4double modz = geo->GetD(fNam, "modz")*mm; // full size in z, mm

  G4cout << "  " << fNam << ", modx, y, z (mm): " << modx << " " << mody << " " << modz << " " << G4endl;

  //module volume
  G4Box *mod_shape = new G4Box(fNam+"_mod", modx/2., mody/2., modz/2.);
  G4Material *mod_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"); // default module material
  G4LogicalVolume *mod_vol = new G4LogicalVolume(mod_shape, mod_mat, mod_shape->GetName());

  //module visibility
  ColorDecoder mod_vis("0:0:1:2");
  mod_vol->SetVisAttributes(mod_vis.MakeVis(geo, fNam, "mod_vis"));

  //module center position
  G4double mod_xpos=0., mod_ypos=0., mod_zpos=0.;
  geo->GetOptD(fNam, "xpos", mod_xpos, GeoParser::Unit(mm)); // center position in x, mm
  geo->GetOptD(fNam, "ypos", mod_ypos, GeoParser::Unit(mm)); // center position in y, mm
  geo->GetOptD(fNam, "zpos", mod_zpos, GeoParser::Unit(mm)); // center position in z, mm

  //module in its mother volume
  G4LogicalVolume *mother_vol = top;
  G4String mother_nam;
  if( geo->GetOptS(fNam, "place_into", mother_nam) ) {
    mother_vol = GetMotherVolume(mother_nam, top);
  }
  new G4PVPlacement(0, G4ThreeVector(mod_xpos, mod_ypos, mod_zpos), mod_vol, mod_vol->GetName(), mother_vol, false, 0);

  //crystal size
  G4double crystal_xy = geo->GetD(fNam, "crystal_xy")*mm; // crystal full size in x and y, mm
  G4double crystal_z = geo->GetD(fNam, "crystal_z")*mm; // crystal full length in z, mm

  //wrapping thickness
  G4double wrapping_thickness = geo->GetD(fNam, "wrapping_thickness")*mm; // wrapping thickness, mm

  //cell holding the wrapped crystal and photon detector
  G4double cell_xy = crystal_xy + 2.*wrapping_thickness; // factor of 2 for wrapping on all sides

  //cell volume, modz = wrapping_thickness + crystal_z + phot_z must hold
  G4Box *cell_shape = new G4Box(fNam+"_cell", cell_xy/2., cell_xy/2., modz/2.);
  G4LogicalVolume *cell_vol = new G4LogicalVolume(cell_shape, mod_mat, cell_shape->GetName()); // default volume material
  cell_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //number of cells in module along x and y
  G4int nx = geo->GetI(fNam, "nx"); // num cells in x
  G4int ny = geo->GetI(fNam, "ny"); // num cells in y

  //cells in module, modx = nx*cell_xy must hold and the same for y
  G4int cell_cnt = 0; // cell count in module
  for(G4int ix=0; ix<nx; ix++) {
    for(G4int iy=0; iy<ny; iy++) {

      //cell position in x and y
      G4double cell_posx = -modx/2. + cell_xy/2. + ix*cell_xy;
      G4double cell_posy = -mody/2. + cell_xy/2. + iy*cell_xy;

      //put the cell in the module
      new G4PVPlacement(0, G4ThreeVector(cell_posx, cell_posy, 0), cell_vol, cell_vol->GetName(), mod_vol, false, cell_cnt++);
    }//iy
  }//ix

  //crystal volume, sensitive volume for energy deposition in the crystal
  G4Box *crystal_shape = new G4Box(fNam, crystal_xy/2., crystal_xy/2., crystal_z/2.);
  G4Material *crystal_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_PbWO4");
  G4LogicalVolume *crystal_vol = new G4LogicalVolume(crystal_shape, crystal_mat, fNam);

  //optical properties for the crystal material and surface
  SetCrystalOptics( crystal_mat );
  SetCrystalSurface( crystal_vol );

  //crystal visibility
  ColorDecoder crystal_vis("1:0:0:2");
  crystal_vol->SetVisAttributes(crystal_vis.MakeVis(geo, fNam, "crystal_vis"));

  //crystal in cell, at the front z of the cell with space for wrapping
  G4VPhysicalVolume *crystal_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*(modz-crystal_z)-wrapping_thickness),
    crystal_vol, fNam, cell_vol, false, 0);

  //wrapping volume, along the sides and the front (larger z) of the crystal
  G4Box *wrapping_outer = new G4Box(fNam+"_wrapping_outer", cell_xy/2., cell_xy/2., (crystal_z+wrapping_thickness)/2.);
  G4SubtractionSolid *wrapping_shape
    = new G4SubtractionSolid(fNam+"_wrapping", wrapping_outer, crystal_shape, 0, G4ThreeVector(0, 0, -wrapping_thickness/2.));

  //wrapping material and logical volume
  G4Material *wrapping_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_TEFLON");
  G4LogicalVolume *wrapping_vol = new G4LogicalVolume(wrapping_shape, wrapping_mat, wrapping_shape->GetName());

  //wrapping visibility
  ColorDecoder wrapping_vis("0:1:0:2");
  wrapping_vol->SetVisAttributes(wrapping_vis.MakeVis(geo, fNam, "wrapping_vis"));

  //wrapping volume in the cell, directly at the front (larger z) of the cell
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*(modz-crystal_z-wrapping_thickness)), wrapping_vol, fNam, cell_vol, false, 0);

  //PMT photocathode radius along xy plane
  G4double radius = geo->GetD(fNam, "radius")*mm; // photocathode radius along xy, mm

  //PMT photocathode thickness along z (mm), accounts both for glass plate and photocathode at the end of lower z
  G4double dz = geo->GetD(fNam, "dz")*mm; // photocathode thickness in z, mm

  //photocathode visibility
  ColorDecoder cath_main_vis("0:0:1:2"); // main photocathode volume, glass and cathode layer
  ColorDecoder cath_lay_vis("1:0:0:2"); // photocathode layer itself

  //create the PMT photocathode
  PhotoCathPMT *cath = new PhotoCathPMT(fNam+"_cath", 0x0);
  fPMT = dynamic_cast<Detector*>( cath );
  G4LogicalVolume *cath_vol = cath->CreateGeometry(radius, dz,
    cath_main_vis.MakeVis(geo, fNam, "cath_main_vis"), cath_lay_vis.MakeVis(geo, fNam, "cath_lay_vis"));

  //PMT photocathode in the cell, placed at the back (lower z)
  new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*modz+0.5*dz), cath_vol, fNam+"_cath", cell_vol, false, 0);

  //optical boundary from the crystal to the PMT glass layer
  SetCrystalBoundary(crystal_phys, cath->GetGlassPhysVol());

}//CalPWO

//_____________________________________________________________________________
void CalPWO::Add(std::vector<Detector*> *vec) {

  //add this detector and its PMT to sensitive detectors

  vec->push_back(this);
  fPMT->Add(vec);

}//Add

//_____________________________________________________________________________
void CalPWO::SetCrystalOptics(G4Material *mat) {

  //scintillation and optical properties for the crystal

  std::vector<G4double> scin_lam = {428, 414}; // nm, scintillation wavelength
  std::vector<G4double> scin_fast = {1., 1.};

  G4MaterialPropertiesTable *tab = new G4MaterialPropertiesTable();
  tab->AddProperty("SCINTILLATIONCOMPONENT1", LambdaNMtoEV(scin_lam), scin_fast);
  tab->AddProperty("SCINTILLATIONCOMPONENT2", LambdaNMtoEV(scin_lam), scin_fast);
  //tab->AddConstProperty("FASTTIMECONSTANT", 1*ps);
  //tab->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1.67*ns);
  //tab->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 6.6*ns);
  //tab->AddConstProperty("SCINTILLATIONYIELD", 300/MeV);
  //tab->AddConstProperty("SCINTILLATIONYIELD1", 0.5);
  //tab->AddConstProperty("SCINTILLATIONYIELD2", 0.5);
  tab->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 0.64*ns); // PWO-UF
  tab->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 3.2*ns);
  tab->AddConstProperty("SCINTILLATIONYIELD", 150/MeV);
  tab->AddConstProperty("SCINTILLATIONYIELD1", 0.72);
  tab->AddConstProperty("SCINTILLATIONYIELD2", 0.28);
  tab->AddConstProperty("RESOLUTIONSCALE", 1.);

  //uniform optical properties
  std::vector<G4double> opt_lam = {800, 350}; // nm
  std::vector<G4double> opt_r = {2.2, 2.2};
  std::vector<G4double> opt_abs = {200*cm, 200*cm};

  tab->AddProperty("RINDEX", LambdaNMtoEV(opt_lam), opt_r);
  tab->AddProperty("ABSLENGTH", LambdaNMtoEV(opt_lam), opt_abs);

  mat->SetMaterialPropertiesTable(tab);

}//SetCrystalOptics

//_____________________________________________________________________________
void CalPWO::SetCrystalSurface(G4LogicalVolume *vol) {

  //crystal optical surface

  G4OpticalSurface *surface = new G4OpticalSurface("CrystalSurface");
  surface->SetType(dielectric_metal);
  //surface->SetFinish(polished);
  surface->SetFinish(ground);
  //surface->SetSigmaAlpha(0.1);
  surface->SetModel(unified);

  new G4LogicalSkinSurface("CrystalSurfaceL", vol, surface);

  //surface material
  std::vector<G4double> opt_lam = {800, 350}; // nm
  std::vector<G4double> reflectivity = {0.9, 0.9};
  std::vector<G4double> efficiency = {1., 1.};
  G4MaterialPropertiesTable *surfmat = new G4MaterialPropertiesTable();
  surfmat->AddProperty("REFLECTIVITY", LambdaNMtoEV(opt_lam), reflectivity);
  surfmat->AddProperty("EFFICIENCY", LambdaNMtoEV(opt_lam), efficiency);
  surface->SetMaterialPropertiesTable(surfmat);
  //csurf->DumpInfo();

}//SetCrystalSurface

//_____________________________________________________________________________
void CalPWO::SetCrystalBoundary(G4VPhysicalVolume *crystal, G4VPhysicalVolume *glass) {

  //optical boundary from the crystal to PMT glass

  G4OpticalSurface *surf = new G4OpticalSurface("CrytalGlassS");
  surf->SetType(dielectric_dielectric); // photons go to the detector, must have rindex defined
  surf->SetFinish(polished);
  //surf->SetModel(unified);
  surf->SetModel(glisur);

  new G4LogicalBorderSurface("CrytalGlassB", crystal, glass, surf);

  std::vector<G4double> opt_lam = {800, 350}; // nm
  std::vector<G4double> reflectivity = {1., 1.};
  std::vector<G4double> efficiency = {1., 1.};

  G4MaterialPropertiesTable *surfmat = new G4MaterialPropertiesTable();
  surfmat->AddProperty("REFLECTIVITY", LambdaNMtoEV(opt_lam), reflectivity);
  surfmat->AddProperty("EFFICIENCY", LambdaNMtoEV(opt_lam), efficiency);
  surf->SetMaterialPropertiesTable(surfmat);

}//SetCrystalBoundary

//_____________________________________________________________________________
std::vector<G4double> CalPWO::LambdaNMtoEV(const std::vector<G4double>& lambda) {

  //converting wavelength (lambda) in nm to energy (en) in eV
  /*

    E = h nu

          nu = c/lambda

    E = h c / lambda

        h = 6.625 x 10^-34 J s

        c = 3 x 10^8 m/s

        1 J = 1.602 x 10^-19 eV

        1 m = 10^9 nm

        h c = ((6.625 x 3)/ 1.602) x 10^-34 x 10^8 x 10^19 x 10^9 = 10^2 x ((6.625 x 3)/ 1.602) = 1240.637 eV nm 

    E (eV) = 1240.637 / lambda (nm)

  */

  //energy in eV
  std::vector<G4double> en;

  //wavelength (nm) loop
  for(auto i: lambda) {

    en.push_back( (1240.637/i)*eV );

  }//wavelength (nm) loop

  return en;

}//LambdaNMtoEV

//_____________________________________________________________________________
G4bool CalPWO::ProcessHits(G4Step *step, G4TouchableHistory*) {

  //cell location in the module
  const G4TouchableHandle& hnd = step->GetPreStepPoint()->GetTouchableHandle();

  //global cell position
  G4ThreeVector origin(0, 0, 0);
  G4ThreeVector gpos = hnd->GetHistory()->GetTopTransform().Inverse().TransformPoint(origin);

  //make the hit for the cell, indexing as copy #1 -> cell, copy #0 -> module
  Int_t cell_id = hnd->GetCopyNumber(1);
  CalPWOHits::Hit& hit = fHits.ConstructedAt(cell_id, CalPWOHits::Hit(cell_id, gpos.x()/mm, gpos.y()/mm, gpos.z()/mm));

  //deposited energy in step, GeV
  G4double edep = step->GetTotalEnergyDeposit()/GeV;

  //increment energy deposit in the hit
  hit.en += edep;

  //primary ID for the track in step
  Int_t prim_id = fStack->GetPrimaryID( step->GetTrack()->GetTrackID() );

  //energy deposit for the given primary ID
  if( hit.prim_energy.find(prim_id) == hit.prim_energy.end() ) {
    //initialize with zero
    hit.prim_energy[prim_id] = 0;
  }
  hit.prim_energy[prim_id] += edep; // increment energy deposit for the primary ID

  return true;

}//ProcessHits

//_____________________________________________________________________________
void CalPWO::CreateOutput(TTree *tree) {

  //create output for the hits
  fHits.CreateOutput(fNam, tree);

  //load the tracking action for primary particle IDs
  fStack = static_cast<const MCParticleAction*>( G4RunManager::GetRunManager()->GetUserTrackingAction() );

  //G4cout << "CalPWO::CreateOutput " << fStack << G4endl;

}//CreateOutput

//_____________________________________________________________________________
void CalPWO::ClearEvent() {

  fHits.ClearEvent();

}//ClearEvent

//_____________________________________________________________________________
void CalPWO::FinishEvent() {

  fHits.FinishEvent();

}//FinishEvent

//_____________________________________________________________________________
G4LogicalVolume* CalPWO::GetMotherVolume(G4String mother_nam, G4LogicalVolume *top) {

  for(size_t i=0; i<top->GetNoDaughters(); i++) {

    G4LogicalVolume *dv = top->GetDaughter(i)->GetLogicalVolume();

    if( dv->GetName() == mother_nam ) {
      return dv;
    }
  }

  return 0x0;

}//GetMotherVolume
























