
#ifndef ParticleCounter_h
#define ParticleCounter_h

// Counter plane for incoming particles

class GeoParser;
class G4Step;

#include "Detector.h"
#include "G4VSensitiveDetector.hh"
#include "ParticleCounterHits.h"

class ParticleCounter: public Detector, public G4VSensitiveDetector {

  public:

    ParticleCounter(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void CreateOutput(TTree *tree);
    virtual void ClearEvent();
    virtual void FinishEvent();

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

  private:

    G4LogicalVolume* GetMotherVolume(G4String mother_nam, G4LogicalVolume *top);

    G4String fNam; //detector name

    G4bool fRemoveTracks; // stop and remove tracks incident on the counter

    //hits
    ParticleCounterHits::Coll fHits;

};

#endif

