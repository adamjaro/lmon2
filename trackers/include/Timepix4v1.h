
#ifndef Timepix4v1_h
#define Timepix4v1_h

#include "Detector.h"
#include "G4VSensitiveDetector.hh"
#include "TrkPlaneBasicHits.h"

class Timepix4v1 : public Detector, public G4VSensitiveDetector {

  public:

    Timepix4v1(const G4String&, GeoParser*, G4LogicalVolume*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void CreateOutput(TTree*);
    virtual void ClearEvent();
    virtual void FinishEvent();

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

  private:

    G4LogicalVolume* MakeSensorLayer(G4int nx, G4int ny, G4double dxy, G4double dz, G4double layx, G4double layy);
    G4LogicalVolume* MakeMaterialLayer(G4double, G4double, G4double, G4String, G4String, G4String, G4String);

    G4String fNam; // name of detector sensitive logical volume
    GeoParser *fGeo; // geometry parser

    const MCParticleAction *fStack; // stack for primary particles

    TrkPlaneBasicHits::Coll fHits; // hits collection

};

#endif

