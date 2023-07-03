
//_____________________________________________________________________________
//
// W/ScFi originally created by Zhiwan Xu (UCLA)
//
// Version 3 with tower assembly along z axis
//_____________________________________________________________________________

//C++
#include "math.h"

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
#include "G4PVReplica.hh"
#include "CLHEP/Vector/RotationY.h"

//local classes
#include "WScFiZXv3.h"
#include "GeoParser.h"

//_____________________________________________________________________________
WScFiZXv3::WScFiZXv3(const G4String& nam, GeoParser *geo, G4LogicalVolume *top) : Detector(),
  G4VSensitiveDetector(nam), fNam(nam), fCheckOverlaps(false) {

  G4cout << "  WScFiZXv3: " << fNam << G4endl;

  //Birks correction
  fBirksCoefficient = 0.126*mm/MeV;
  geo->GetOptD(nam, "Birks_coefficient", fBirksCoefficient, GeoParser::Unit(mm/MeV));

  //dimensions for single tower
  G4double towerSizeXY = 25*mm;
  geo->GetOptD(nam, "towerSizeXY", towerSizeXY, GeoParser::Unit(mm));
  G4double towerEMZ = 170*mm;
  geo->GetOptD(nam, "towerEMZ", towerEMZ, GeoParser::Unit(mm));
  G4double zpos = towerEMZ/2;
  geo->GetOptD(nam, "zpos", zpos, GeoParser::Unit(mm));

  G4int nxy = 32;
  geo->GetOptI(nam, "nxy", nxy);

  //module size for tower assembly, increased to allow for tower rotation
  G4double modxy = nxy*towerSizeXY; //  + 40*mm
  G4double modz = towerEMZ; //  + 6*mm

  G4cout << "    modxy: " << modxy << G4endl;
  G4cout << "    modz: " << modz << G4endl;

  //top calorimeter volume
  G4Box *mods = new G4Box(fNam+"_mod", modxy/2, modxy/2, modz/2);
  G4LogicalVolume *modv = new G4LogicalVolume(mods, G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"), fNam+"_mod");
  new G4PVPlacement(0, G4ThreeVector(0, 0, zpos), modv, fNam+"_mod", top, false, 0);
  //modv->SetVisAttributes(new G4VisAttributes(G4Color(0, 0, 1)));
  modv->SetVisAttributes( G4VisAttributes::GetInvisible() );

  G4LogicalVolume *towv = MakeTower(towerSizeXY, towerEMZ);

  G4int tcnt = 0;
  G4double xypos0 = -(nxy*towerSizeXY)/2 + towerSizeXY/2;
  for(G4int ix=0; ix<nxy; ix++) {
    for(G4int iy=0; iy<nxy; iy++) {
      G4double xpos = xypos0 + ix*towerSizeXY;
      G4double ypos = xypos0 + iy*towerSizeXY;

      new G4PVPlacement(0, G4ThreeVector(xpos, ypos, 0), towv, towv->GetName(), modv, false, tcnt++);
    }
  }

}//WScFiZXv3

//_____________________________________________________________________________
G4LogicalVolume* WScFiZXv3::MakeTower(G4double calorSizeXY, G4double calorEMZ) {

  const double offset=0.5;//in mm
  const double dist=1.0;
  const double tot_len=calorSizeXY;
  const double h=0.5*sqrt(3)*dist;

  const int nx1=int((tot_len-2*offset)/(dist/2))+1;
  const int ny1=int((tot_len-offset)/(2*h))+1;
  const int ny2=int((tot_len-offset-h)/(2*h))+1;

  const double x0=-((tot_len/2.0)-offset);
  const double y01=((tot_len/2.0)-offset);
  const double y02=((tot_len/2.0)-offset-h);

  //G4Material* defaultMaterial = G4Material::GetMaterial("Galactic");
  //G4Material* gapMaterial2 =G4Material::GetMaterial("G4_POLYSTYRENE");
  G4Material* defaultMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
  G4Material* gapMaterial2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4Material* EMCal_abs_mat=new G4Material("EMCal_fiber_mat",12.4*g/cm3,2);
  G4double a=183.85*g/mole;
  G4Element* elW=new G4Element("Tungsten","W",74.,a);
  EMCal_abs_mat->AddElement(elW,96.0*perCent);
  EMCal_abs_mat->AddMaterial(gapMaterial2,4.0*perCent);


  //EM
  G4LogicalVolume* calorEM;
  G4VSolid* calorimeterEM = new G4Box("CalorimeterEM_", calorSizeXY/2, calorSizeXY/2, calorEMZ/2);
  calorEM=new G4LogicalVolume(calorimeterEM,defaultMaterial, fNam+"_CalorEM");
    
  //absorber
  G4LogicalVolume* absorberEMLV;
  G4VSolid* absorberEM = new G4Box(fNam+"_AbsoEM_solid",            // its name
                                     calorSizeXY/2, calorSizeXY/2, calorEMZ/2); // its size
  absorberEMLV = new G4LogicalVolume(absorberEM,EMCal_abs_mat, fNam+"_AbsoEM");
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),absorberEMLV, fNam+"_AbsoEM_p",calorEM,false,0,fCheckOverlaps);
    
  G4VSolid* gapEM = new G4Tubs(fNam+"_GapEM",             // its name
                               0.0, 0.235*mm, calorEMZ/2,0.0,360.0 * deg); // its size//0.0 * deg, 360.0 * deg
    
  //Fibers
  G4LogicalVolume* gapEMLV;
  gapEMLV = new G4LogicalVolume(gapEM,gapMaterial2, fNam);
  int copynono=0;
  double step_x=(dist/2.0)*mm;
  double step_y=(2.0*h)*mm;
  //G4cout<<"Nx: "<<nx1<<" Ny1: "<<ny1<<" Ny2: "<<ny2<<" step_x: "<<step_x<<" step_y: "<<step_y<<" x0: "<<x0<<" y01: "<<y01<<" y02: "<<y02<<G4endl;

  for(int i=0;i<nx1;i++){
    G4double pos_x=x0*mm+i*step_x;
    G4double pos_y=0.0;
    //if(i==(nFibAr-1)) continue;
    // if(i%2==0) pos_x=(-29.95+i*(0.05))*cm;
    // if(i%2!=0) pos_x=(-29.9+i*(0.05))*cm;
    // pos_x=(-29.95+i*(0.1))*cm;
    int jend=(i%2==0) ? ny1 : ny2;
    for(int j=0;j<jend;j++){
            
            
      if(i%2==0) pos_y=y01*mm-j*step_y;
      if(i%2!=0) pos_y=y02*mm-j*step_y;
      new G4PVPlacement(0,G4ThreeVector(pos_x,pos_y, 0.),gapEMLV,fNam+"_EMGapPhysical",absorberEMLV,false,copynono,fCheckOverlaps);
      //new G4PVPlacement(0,G4ThreeVector(pos_x,pos_y, 0.),gapEMLV,"EMGapPhysical",absorberEMLV,0,copynono);//first try
      copynono++;
      //G4cout<<"Point # "<<copynono<<" x: "<<pos_x<<" y: "<<pos_y<<G4endl;

    }
  }

  absorberEMLV->SetVisAttributes( G4VisAttributes::GetInvisible() );
  gapEMLV->SetVisAttributes( G4VisAttributes::GetInvisible() );

  G4VisAttributes* calorEMvis= new G4VisAttributes(G4Colour(1,0,1));//magenta calorimeter
  calorEMvis->SetForceAuxEdgeVisible(true);
  calorEM->SetVisAttributes(calorEMvis);

  return calorEM;

}//MakeTower

//_____________________________________________________________________________
G4bool WScFiZXv3::ProcessHits(G4Step *step, G4TouchableHistory*) {

  //deposited energy in step (GeV)
  fEdep += BirksCorrectedEnergyDeposit(step);

  return true;

}//ProcessHits

//_____________________________________________________________________________
G4double WScFiZXv3::BirksCorrectedEnergyDeposit(G4Step *step) {

  //deposited energy in step, MeV
  G4double edep_step = step->GetTotalEnergyDeposit()/MeV;

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
void WScFiZXv3::ClearEvent() {

  fEdep = 0;

}//ClearEvent

//_____________________________________________________________________________
void WScFiZXv3::CreateOutput(TTree *tree) {

  tree->Branch((fNam+"_edep").c_str(), &fEdep, (fNam+"_edep/D").c_str());

}//CreateOutput













