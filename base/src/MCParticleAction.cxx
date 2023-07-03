
//ROOT
#include "TTree.h"

//Geant
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

//local classes
#include "MCParticleAction.h"
#include "MCParticles.h"

using namespace std;

//_____________________________________________________________________________
MCParticleAction::MCParticleAction() : G4UserTrackingAction() {

  G4cout << "MCParticleAction::MCParticleAction " << this << G4endl;

  //create the local stack
  fStack = new unordered_map<G4int, G4int>;

  //output on MC particles
  fPart = new MCParticles::Coll();

}//MCParticleAction

//_____________________________________________________________________________
void MCParticleAction::PreUserTrackingAction(const G4Track *track) {

  //current track ID
  G4int trackID = track->GetTrackID();

  //track already seen
  if( fStack->find(trackID) != fStack->end() ) return;

  //parent ID for current track
  G4int parentID = track->GetParentID();

  //G4cout << "MCParticleAction::PreUserTrackingAction, track ID: " << trackID << " " << parentID << G4endl;

  //primary track, add the track ID as mapped value
  if( parentID == 0 ) {

    fStack->emplace(trackID, trackID);
    AddMCParticle(track);
    return;
  }

  //secondary track, assign its primary ID as mapped value
  fStack->emplace(trackID, fStack->at(parentID));

  //G4cout << "MCParticleAction::PreUserTrackingAction, prim ID:  " << fStack->at(trackID) << G4endl;

}//PreUserTrackingAction

//_____________________________________________________________________________
G4int MCParticleAction::GetPrimaryID(G4int id) const {

  //retrieve ID of primary particle belonging to a track at id in the argument

  return fStack->at(id);

}//GetPrimaryID

//_____________________________________________________________________________
void MCParticleAction::Reset() const {

  //G4cout << "MCParticleAction::Reset " << this << G4endl;

  //remove all particles
  fStack->clear();

  //clear the MC particles
  fPart->ClearEvent();

}//Reset

//_____________________________________________________________________________
void MCParticleAction::CreateOutput(TTree *tree) const {

  //output on MC particles
  fPart->CreateOutput("mcp", tree);

}//CreateOutput

//_____________________________________________________________________________
void MCParticleAction::AddMCParticle(const G4Track *track) {

  //add MC particle from the track

  //particle kinematics
  const G4DynamicParticle *part = track->GetDynamicParticle();

  //vertex position
  G4ThreeVector vtx = track->GetVertexPosition();

  //add the MC particle
  MCParticles::Part& mc = fPart->Add( MCParticles::Part() );
  mc.pdg = part->GetPDGcode();
  mc.itrk = track->GetTrackID();
  mc.en = part->Get4Momentum().e()/GeV;
  mc.theta = part->Get4Momentum().theta()/rad;
  mc.phi = part->Get4Momentum().phi()/rad;
  mc.vx = vtx.x()/mm;
  mc.vy = vtx.y()/mm;
  mc.vz = vtx.z()/mm;

  //G4cout << "MCParticleAction::AddMCParticle " << track->GetTrackID() << " " << part->GetPDGcode() << " ";
  //G4cout << vtx.x()/mm << " " << vtx.y()/mm << " " << vtx.z()/mm << " " << part->Get4Momentum().e()/GeV << G4endl;

}//AddMCParticle

//_____________________________________________________________________________
void MCParticleAction::FinishEvent() const {

  fPart->FinishEvent();

}//FinishEvent

































