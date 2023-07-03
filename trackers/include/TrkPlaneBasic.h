
#ifndef TrkPlaneBasic_h
#define TrkPlaneBasic_h

// MAPS tracking layer, basic implementation

#include "Detector.h"
#include "G4VSensitiveDetector.hh"
#include "TrkPlaneBasicHits.h"

class GeoParser;
class MCParticleAction;

class TrkPlaneBasic : public Detector, public G4VSensitiveDetector {

  public:

    TrkPlaneBasic(const G4String&, GeoParser*, G4LogicalVolume*);
    virtual ~TrkPlaneBasic() {}

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void CreateOutput(TTree*);
    virtual void ClearEvent();
    virtual void FinishEvent();

  private:

    G4LogicalVolume* GetMotherVolume(G4String mother_nam, G4LogicalVolume *top);

    G4String fNam; // name of detector sensitive logical volume

    const MCParticleAction *fStack; // stack for primary particles

    //hits
    TrkPlaneBasicHits::Coll fHits;

};

#endif



















