
#ifndef MCParticleAction_h
#define MCParticleAction_h

#include<unordered_map>

#include "Rtypes.h"

#include "G4UserTrackingAction.hh"

class TTree;

namespace MCParticles {
  class Coll;
}

class MCParticleAction : public G4UserTrackingAction {

  public:

    MCParticleAction();

    void PreUserTrackingAction(const G4Track *) override;

    void Reset() const;

    G4int GetPrimaryID(G4int id) const;

    void CreateOutput(TTree *tree) const;

    void FinishEvent() const;

  private:

    std::unordered_map<G4int, G4int> *fStack; // local stack for primary particle IDs

    void AddMCParticle(const G4Track*);

    //output on MC particles
    MCParticles::Coll *fPart;

};

#endif

