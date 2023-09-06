
//_____________________________________________________________________________
//
// Simple optical Si detector imitating a PIN photodiode
//
//_____________________________________________________________________________

//C++

//ROOT
#include "TTree.h"

//Geant
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"

//local classes
#include "GeoParser.h"
#include "ColorDecoder.h"
#include "OpSiDet.h"

//_____________________________________________________________________________
OpSiDet::OpSiDet(const G4String& nam, GeoParser *geo, G4LogicalVolume *top): Detector(),
   G4VSensitiveDetector(nam), fNam(nam) {

  G4cout << "OpSiDet: " << fNam << G4endl;

}//OpSiDet

//_____________________________________________________________________________
G4LogicalVolume* OpSiDet::CreateGeometry(G4double dx, G4double dy, G4double dz, G4VisAttributes *vis) {

  //detector shape
  G4Box *shape = new G4Box(fNam, dx/2., dy/2., dz/2.);
  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  //optical properties for detector material
  G4MaterialPropertiesTable *tab = new G4MaterialPropertiesTable();
  tab->AddProperty("RINDEX", "Fused Silica");
  mat->SetMaterialPropertiesTable(tab);

  G4LogicalVolume *vol = new G4LogicalVolume(shape, mat, fNam);

  vol->SetVisAttributes( vis );

  return vol;

}//CreateGeometry

//_____________________________________________________________________________
void OpSiDet::MakeBoundary(G4VPhysicalVolume *src_phys, G4VPhysicalVolume *opdet_phys) {

  //optical boundary from photon source volume to the detector volume

  G4OpticalSurface *surf = new G4OpticalSurface("OpSiDet_surface");
  surf->SetType(dielectric_dielectric); // photons go to the detector, must have rindex defined
  surf->SetFinish(polished);
  //surf->SetModel(unified);
  surf->SetModel(glisur);

  new G4LogicalBorderSurface("OpSiDet_surface", src_phys, opdet_phys, surf);

  std::vector<G4double> opt_lam = {800, 350}; // nm
  std::vector<G4double> reflectivity = {1, 1};
  std::vector<G4double> efficiency = {1., 1.};

  G4MaterialPropertiesTable *surfmat = new G4MaterialPropertiesTable();
  surfmat->AddProperty("REFLECTIVITY", LambdaNMtoEV(opt_lam), reflectivity);
  surfmat->AddProperty("EFFICIENCY", LambdaNMtoEV(opt_lam), efficiency);
  surf->SetMaterialPropertiesTable(surfmat);

}//MakeBoundary

//_____________________________________________________________________________
G4bool OpSiDet::ProcessHits(G4Step *step, G4TouchableHistory*) {

  G4cout << "OpSiDet::ProcessHits" << G4endl;

  return true;

}//ProcessHits

//_____________________________________________________________________________
std::vector<G4double> OpSiDet::LambdaNMtoEV(const std::vector<G4double>& lambda) {

  //converting wavelength (lambda) in nm to energy (en) in eV
  /*

    E = h nu,  nu = c/lambda

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




















