
//_____________________________________________________________________________
//
// generic run action, handles running time
// and output creation in DetectorConstruction
//_____________________________________________________________________________

//C++ headers
#include <time.h>

//Geant headers
#include "G4RunManager.hh"

//local headers
#include "RunAction.h"
#include "DetectorConstruction.h"

//_____________________________________________________________________________
RunAction::RunAction() : G4UserRunAction(), fDet(0) {

  //get the detector
  fDet = static_cast<const DetectorConstruction*>( G4RunManager::GetRunManager()->GetUserDetectorConstruction() );

}//RunAction

//_____________________________________________________________________________
void RunAction::BeginOfRunAction(const G4Run*) {

  fStart = clock();

  fDet->CreateOutput();

}//BeginOfRunAction

//_____________________________________________________________________________
void RunAction::EndOfRunAction(const G4Run*) {

  G4cout << "Running time: " << (clock() - fStart)/CLOCKS_PER_SEC << " sec" << G4endl;  

}//EndOfRunAction
















