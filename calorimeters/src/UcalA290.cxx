
//_____________________________________________________________________________
//
// Uranium-scintillator calorimeter for ZEUS according to NIM A290 (1990) 95-108
//
//_____________________________________________________________________________

//C++
#include "math.h"

//ROOT
#include "TTree.h"
#include "TH1F.h"

//Geant
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4VisAttributes.hh"
#include "G4PVReplica.hh"

//local classes
#include "UcalA290.h"
#include "GeoParser.h"

//_____________________________________________________________________________
UcalA290::UcalA290(const G4String& nam, GeoParser *geo, G4LogicalVolume *top) : Detector(),
  G4VSensitiveDetector(nam), fNam(nam), fGeo(geo) {

  G4cout << "  UcalA290: " << fNam << G4endl;

  DefineMaterials();

  fModXY = 800.; // module transverse size (800), mm

  G4double al_z = 15; // aluminum front plate, thickness in z, mm

  //scintillator layer
  G4double scin_z = 2.6; // scintillator thickness in z, mm
  G4double spacer_z = 4.; // spacer for scintillators, thickness in z, mm

  //absorber layers
  G4double abso_z = 3.3; // uranium thickness, mm
  geo->GetOptD(nam, "abso_z", abso_z);
  G4cout << "    abso_z: " << abso_z << G4endl;
  G4double clad_emc_z = 0.2; // cladding in EMC layers, mm
  G4double clad_hac_z = 0.4; // cladding in HAC layers, mm

  G4double gap_emc = 10; // gap for semiconductor detectors in EMC section, mm

  G4int nLayEMC = 25; // number of EMC layers, after the front Al plate with its scintillator layer
  geo->GetOptI(nam, "n_emc", nLayEMC);

  G4int nLayHAC = 80; // number of layers in each HAC1 and HAC2 sections
  geo->GetOptI(nam, "n_hac", nLayHAC);

  G4bool use_clad = true; // use cladding for absorber layers
  geo->GetOptB(nam, "use_clad", use_clad);
  G4cout << "    use_clad: " << use_clad << G4endl;

  //Birks correction
  fUseBirksCorrection = true;
  fBirksCoefficient = 0.126*mm/MeV;
  geo->GetOptB(nam, "use_Birks_correction", fUseBirksCorrection);
  geo->GetOptD(nam, "Birks_coefficient", fBirksCoefficient);

  G4cout << "    use_Birks_correction: " << fUseBirksCorrection << G4endl;
  G4cout << "    Birks_coefficient: " << fBirksCoefficient << G4endl;

  fMaxTime = 50.; //maximal time for signal integration, ns
  geo->GetOptD(nam, "max_time", fMaxTime);
  fMaxTime = fMaxTime*ns;
  G4cout << "    max_time: " << fMaxTime << G4endl;

  //module length along z, given by the sum of all layers
  G4double modz = al_z + spacer_z + nLayEMC*(abso_z+2*clad_emc_z+spacer_z) + 2*gap_emc + nLayHAC*2*(abso_z+2*clad_hac_z+spacer_z);

  G4cout << "    Module length (mm): " << modz << G4endl;

  //indices of the first layers in HAC1 and HAC2 sections
  fStartHAC1 = nLayEMC+1; // includes the layer after the front Al plate (26)
  fStartHAC2 = fStartHAC1 + nLayHAC; // (106)
  G4cout << "    HAC1 and HAC2: " << fStartHAC1 << " " << fStartHAC2 << G4endl;

  //positions of scintillator centers along the module
  std::vector<G4double> scin_pos;

  //event outputs
  fELayer = new std::vector<Float_t>(nLayEMC+1+2*nLayHAC); // deposited energy in individual layers
  fTimeStep = new TH1F("time_step", "time_step", 20, 0, 200); // time in step, ns

  //top calorimeter volume, box shape
  G4String modnam = fNam+"_mod"; // module name
  G4Box *mods = new G4Box(modnam, fModXY/2, fModXY/2, modz/2);

  //module logical volume with default material
  G4Material *mod_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *modv = new G4LogicalVolume(mods, mod_mat, modnam);
  new G4PVPlacement(0, G4ThreeVector(0, 0, modz/2), modv, modnam, top, false, 0);
  //new G4PVPlacement(0, G4ThreeVector(0, 0, 0), modv, modnam, top, false, 0);
  //new G4PVPlacement(0, G4ThreeVector(0, 0, -modz/2), modv, modnam, top, false, 0);

  //front aluminum plate
  MakeAlFront(modv, al_z);

  //scintillator layer
  G4LogicalVolume *scinv = MakeScinLayer(mod_mat, scin_z, spacer_z);

  //absorber layers in EMC and HAC sections
  G4LogicalVolume *abso_vol_emc, *abso_vol_hac;
  if(use_clad) {
    abso_vol_emc = MakeAbsoLayer(abso_z, clad_emc_z, "emc");
    abso_vol_hac = MakeAbsoLayer(abso_z, clad_hac_z, "hac");
  } else {
    abso_vol_emc = MakeAbsoNoClad(abso_z, "emc");
    abso_vol_hac = MakeAbsoNoClad(abso_z, "hac");
  }


  //scintillator positioning
  G4double scin_hz = dynamic_cast<G4Box*>(scinv->GetSolid())->GetZHalfLength(); // half of z for scintillator layer
  G4int scin_cnt = 0; // counter for scintillator layers as placed to the module

  //absorber positioning
  G4double abso_emc_hz = dynamic_cast<G4Box*>(abso_vol_emc->GetSolid())->GetZHalfLength(); // half z for EMC layer
  G4double abso_hac_hz = dynamic_cast<G4Box*>(abso_vol_hac->GetSolid())->GetZHalfLength(); // half z for HAC layer
  G4int emc_cnt = 0; // counter for EMC layers
  G4int hac_cnt = 0; // counter for HAC layers

  //individulal components in the module
  G4double current_z = -modz/2 + al_z; // z position along the module, initialized to the point after the front Al plate

  //first scintillator layer after aluminum front plate
  new G4PVPlacement(0, G4ThreeVector(0, 0, current_z+scin_hz), scinv, scinv->GetName(), modv, false, scin_cnt++);
  scin_pos.push_back( current_z+scin_hz+modz/2 );
  current_z += 2*scin_hz; // move the z position to the point after the first scintillator plate

  //G4cout << "EMC layers:" << G4endl;

  //EMC layers loop
  for(G4int i=0; i<nLayEMC; i++) {

    //absorber
    new G4PVPlacement(0, G4ThreeVector(0, 0, current_z+abso_emc_hz), abso_vol_emc, abso_vol_emc->GetName(), modv, false, emc_cnt++);
    current_z += 2*abso_emc_hz;

    //G4cout << scin_cnt << " " << current_z+scin_hz+modz/2 << G4endl;

    //scintillator
    new G4PVPlacement(0, G4ThreeVector(0, 0, current_z+scin_hz), scinv, scinv->GetName(), modv, false, scin_cnt++);
    scin_pos.push_back( current_z+scin_hz+modz/2 );
    current_z += 2*scin_hz;

    //gaps after layers 3 and 6
    if(i == 2) current_z += gap_emc;
    if(i == 5) current_z += gap_emc;

  }//EMC layers loop

  //G4cout << "HAC layers:" << G4endl;

  //HAC layers loop
  for(G4int i=0; i<2*nLayHAC; i++) {

    //absorber
    new G4PVPlacement(0, G4ThreeVector(0, 0, current_z+abso_hac_hz), abso_vol_hac, abso_vol_hac->GetName(), modv, false, hac_cnt++);
    current_z += 2*abso_hac_hz;

    //G4cout << scin_cnt << " " << current_z+scin_hz+modz/2 << G4endl;

    //scintillator
    new G4PVPlacement(0, G4ThreeVector(0, 0, current_z+scin_hz), scinv, scinv->GetName(), modv, false, scin_cnt++);
    scin_pos.push_back( current_z+scin_hz+modz/2 );
    current_z += 2*scin_hz;

  }//HAC layers loop
/*
  //print the scintillator positions
  G4cout << "---- Scintillator positions ----" << G4endl;
  for(unsigned i=0; i<scin_pos.size(); i++) {
    //G4cout << i << " " << scin_pos[i] << G4endl;
    G4cout << scin_pos[i] << G4endl;
  }
  G4cout << "---- Scintillator positions ----" << G4endl;
*/
}//UcalA290

//_____________________________________________________________________________
void UcalA290::MakeAlFront(G4LogicalVolume *modv, G4double al_z) {

  //front aluminum plate

  G4String nam = fNam+"_front";
  G4Box *shape = new G4Box(nam, fModXY/2, fModXY/2, al_z/2);

  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4LogicalVolume *vol = new G4LogicalVolume(shape, mat, nam);

  //visibility for aluminum front plate
  G4VisAttributes *vis = new G4VisAttributes();
  vis->SetColor(0, 1, 0, 0.5);
  vis->SetForceSolid(true);
  vol->SetVisAttributes(vis);

  //put the front plate to the module
  G4double hz = dynamic_cast<G4Box*>(modv->GetSolid())->GetZHalfLength();
  new G4PVPlacement(0, G4ThreeVector(0, 0, (al_z/2)-hz), vol, nam, modv, false, 0);

}//MakeAlFront

//_____________________________________________________________________________
G4LogicalVolume *UcalA290::MakeScinLayer(G4Material *mod_mat, G4double scin_z, G4double spacer_z) {

  //scintillator layer defined by the spacer
  G4String nam = fNam+"_scin_layer";
  G4Box *layer_shape = new G4Box(nam, fModXY/2, fModXY/2, spacer_z/2);
  G4LogicalVolume *layer_vol = new G4LogicalVolume(layer_shape, mod_mat, nam);
  layer_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //scintillator active volume
  G4Box *scin_shape = new G4Box(fNam, fModXY/2, fModXY/2, scin_z/2);
  G4Material *scin_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4LogicalVolume *scin_vol = new G4LogicalVolume(scin_shape, scin_mat, fNam);

  //visibility for active scintillator
  G4VisAttributes *vis = new G4VisAttributes();
  //vis->SetColor(0, 0, 1, 0.5);
  vis->SetColor(0, 0, 1);
  vis->SetForceSolid(true);
  scin_vol->SetVisAttributes(vis);

  //put the scintillator to the middle of spacer layer
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), scin_vol, fNam, layer_vol, false, 0);

  return layer_vol;

}//MakeScinLayer

//_____________________________________________________________________________
G4LogicalVolume *UcalA290::MakeAbsoLayer(G4double abso_z, G4double clad_z, G4String eh) {

  //absorber layer with cladding
  G4String nam = fNam+"_abso_layer_"+eh;
  G4Box *layer_shape = new G4Box(nam, fModXY/2, fModXY/2, (abso_z+2*clad_z)/2);
  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  G4LogicalVolume *layer_vol = new G4LogicalVolume(layer_shape, mat, nam);

  //absorber visibility for the entire layer
  G4VisAttributes *vis = new G4VisAttributes();
  if(eh == "hac") {
    vis->SetColor(1, 0, 0); // RGB
  } else {
    vis->SetColor(1, 1, 0); // yellow
  }
  vis->SetForceSolid(true);
  layer_vol->SetVisAttributes(vis);

  //DU absorber plate
  G4Box *abso_shape = new G4Box(fNam+"_abso_du_plate_"+eh, fModXY/2, fModXY/2, abso_z/2);

  G4String abso_mat_name = "G4_U";
  fGeo->GetOptS(fNam, "abso_mat_name", abso_mat_name);
  G4cout << "    abso_mat_name: " << abso_mat_name << G4endl;
  G4Material *abso_mat = G4NistManager::Instance()->FindOrBuildMaterial(abso_mat_name);

  G4LogicalVolume *abso_vol = new G4LogicalVolume(abso_shape, abso_mat, abso_shape->GetName());
  abso_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //put the DU plate in the middle of absorber layer
  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), abso_vol, abso_vol->GetName(), layer_vol, false, 0);

  //stainless steel cladding
  G4Box *clad_shape = new G4Box(fNam+"_abso_clad_"+eh, fModXY/2, fModXY/2, clad_z/2);
  G4String clad_mat_name = "G4_STAINLESS-STEEL";
  fGeo->GetOptS(fNam, "clad_mat_name", clad_mat_name);
  G4cout << "    clad_mat_name: " << clad_mat_name << G4endl;
  G4Material *clad_mat = G4NistManager::Instance()->FindOrBuildMaterial(clad_mat_name);
  G4LogicalVolume *clad_vol = new G4LogicalVolume(clad_shape, clad_mat, clad_shape->GetName());
  clad_vol->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //put the cladding to both sides of DU plate, front and back of the layer
  G4double pos0 = -layer_shape->GetZHalfLength() + clad_z/2;
  G4double pos1 = layer_shape->GetZHalfLength() - clad_z/2;
  new G4PVPlacement(0, G4ThreeVector(0, 0, pos0), clad_vol, clad_vol->GetName(), layer_vol, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0, 0, pos1), clad_vol, clad_vol->GetName(), layer_vol, false, 1);

  return layer_vol;

}//MakeAbsoLayer

//_____________________________________________________________________________
G4LogicalVolume *UcalA290::MakeAbsoNoClad(G4double abso_z, G4String eh) {

  //uranium absorber layer with no cladding

  G4String nam = fNam+"_abso_layer_"+eh;
  G4Box *layer_shape = new G4Box(nam, fModXY/2, fModXY/2, abso_z/2);

  G4String abso_mat_name = "G4_U";
  fGeo->GetOptS(fNam, "abso_mat_name", abso_mat_name);
  G4cout << "    abso_mat_name: " << abso_mat_name << G4endl;
  G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial(abso_mat_name);

  G4LogicalVolume *layer_vol = new G4LogicalVolume(layer_shape, mat, nam);

  //absorber visibility for the entire layer
  G4VisAttributes *vis = new G4VisAttributes();
  if(eh == "hac") {
    vis->SetColor(1, 1, 0); // yellow
  } else {
    vis->SetColor(1, 0, 0); // RGB
  }
  vis->SetForceSolid(true);
  layer_vol->SetVisAttributes(vis);

  return layer_vol;

}//MakeAbsoNoClad

//_____________________________________________________________________________
G4bool UcalA290::ProcessHits(G4Step *step, G4TouchableHistory*) {

  //step time
  G4double time = step->GetPostStepPoint()->GetGlobalTime()/ns;
  fTimeStep->Fill( time );

  //maximal time for signal integration
  if( time > fMaxTime ) return true;

  //scintillator location
  int lay_id = 1;
  const G4TouchableHandle& hnd = step->GetPreStepPoint()->GetTouchableHandle();

  //energy in individual layers
  fELayer->at(hnd->GetCopyNumber(lay_id)) += BirksCorrectedEnergyDeposit(step);

  //G4cout << hnd->GetCopyNumber(lay_id);// << G4endl;
  //G4cout << " " << step->GetPreStepPoint()->GetPosition().z()  << G4endl;

  return true;

}//ProcessHits

//_____________________________________________________________________________
G4double UcalA290::BirksCorrectedEnergyDeposit(G4Step *step) {

  //deposited energy in step, MeV
  G4double edep_step = step->GetTotalEnergyDeposit()/MeV;

  if(!fUseBirksCorrection) return edep_step/GeV; // to GeV

  //Birks attenuation
  G4double step_length = step->GetStepLength()*mm;
  G4double particle_charge = step->GetTrack()->GetDefinition()->GetPDGCharge();

  //value of fBirksCoefficient is in mm/MeV
  G4double edep = edep_step;
  if( std::abs(edep_step*step_length*particle_charge) > 1e-12 ) {

    edep = edep_step/(1. + fBirksCoefficient*edep_step/step_length);
  }

  return edep/GeV; //to GeV

}//BirksCorrectedEnergyDeposit

//_____________________________________________________________________________
void UcalA290::DefineMaterials() {

  G4String name, symbol;
  G4double a, z, density, abundance, fraction;
  G4int iz, n, ncomponents;

  //depleted uranium, DESY 89-128 1989
  if( !G4Material::GetMaterial("DepletedUranium", false) ) {
    G4Isotope *U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
    G4Isotope *U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

    G4Element *elDU = new G4Element(name="depleted Uranium", symbol="U", ncomponents=2);
    elDU->AddIsotope(U8, abundance = 99.797*perCent);
    elDU->AddIsotope(U5, abundance = 0.203*perCent);

    G4Element *elNb = new G4Element(name="Nb", symbol="Nb", z = 41., a = 92.906*g/mole);

    G4Material *mDU = new G4Material("DepletedUranium", density=18.5*g/cm3, ncomponents=2);
    mDU->AddElement(elDU, fraction=98.6*perCent);
    mDU->AddElement(elNb, fraction=1.4*perCent);

    //G4cout << G4NistManager::Instance()->FindOrBuildMaterial("DepletedUranium") << G4endl;
  }

}//DefineMaterials

//_____________________________________________________________________________
void UcalA290::ClearEvent() {

  fEdepEMC = 0.;
  fEdepHAC1 = 0.;
  fEdepHAC2 = 0.;

  std::fill(fELayer->begin(), fELayer->end(), 0.);

  fTimeStep->Reset();

}//ClearEvent

//_____________________________________________________________________________
void UcalA290::FinishEvent() {

  //G4cout << "Event:" << G4endl;

  //energy in individual layers

  //layer loop
  for(unsigned i=0; i<fELayer->size(); i++) {
    //G4cout << "i: " << i << " " << fELayer->at(i) << G4endl;

    Float_t en = fELayer->at(i);

    if(i < fStartHAC1) {
      fEdepEMC += en;
    } else if(i >= fStartHAC2) {
      fEdepHAC2 += en;
    } else {
      fEdepHAC1 += en;
    }

  }//layer loop

  //G4cout << "Event: " << fEdepEMC << " " << fEdepHAC1 << " " << fEdepHAC2 << G4endl;

}//FinishEvent

//_____________________________________________________________________________
void UcalA290::CreateOutput(TTree *tree) {

  tree->Branch((fNam+"_edep_EMC").c_str(), &fEdepEMC, (fNam+"_edep_EMC/D").c_str());
  tree->Branch((fNam+"_edep_HAC1").c_str(), &fEdepHAC1, (fNam+"_edep_HAC1/D").c_str());
  tree->Branch((fNam+"_edep_HAC2").c_str(), &fEdepHAC2, (fNam+"_edep_HAC2/D").c_str());

  tree->Branch((fNam+"_edep_layers").c_str(), fELayer);

  tree->Branch((fNam+"_time_step").c_str(), fTimeStep);

}//CreateOutput





















