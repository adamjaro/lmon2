
//_____________________________________________________________________________
//
// standard event action,
// calls ClearEvent and FinishEvent in DetectorConstruction
//_____________________________________________________________________________

//Geant headers
#include "G4Event.hh"
#include "G4RunManager.hh"

//local headers
#include "EventAction.h"
#include "DetectorConstruction.h"
#include "MCParticleAction.h"

//_____________________________________________________________________________
EventAction::EventAction() : G4UserEventAction(), fDet(0), fStack(0) {

  //get the detector
  fDet = static_cast<const DetectorConstruction*>( G4RunManager::GetRunManager()->GetUserDetectorConstruction() );

  //tracking action maintaining local particle stack
  fStack = static_cast<const MCParticleAction*>( G4RunManager::GetRunManager()->GetUserTrackingAction() );

  G4cout << "EventAction::EventAction " << fDet << " " << fStack << G4endl;

}//EventAction

//_____________________________________________________________________________
void EventAction::BeginOfEventAction(const G4Event *evt) {

  //set MC and clear the detectors
  fDet->BeginEvent(evt);

  //reset the TrackingAction for the event
  fStack->Reset();

}//BeginOfEventAction

//_____________________________________________________________________________
void EventAction::EndOfEventAction(const G4Event*) {

  //finish for MC particles
  fStack->FinishEvent();

  //write detector event output
  fDet->FinishEvent();

}//EndOfEventAction














