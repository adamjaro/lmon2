
#ifndef CaloBPC_h
#define CaloBPC_h

// tungsten/scintillator calorimeter following the ZEUS BPC

#include "Detector.h"
#include "G4VSensitiveDetector.hh"
#include "CaloBPCHits.h"

class GeoParser;

class CaloBPC : public Detector, public G4VSensitiveDetector {

  public:

    CaloBPC(const G4String&, GeoParser *, G4LogicalVolume*);
    virtual ~CaloBPC() {}

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void CreateOutput(TTree *);
    virtual void ClearEvent();
    virtual void FinishEvent();

  private:

    G4LogicalVolume* GetMotherVolume(G4String mother_nam, G4LogicalVolume *top);

    G4String fNam; // name of detector sensitive logical volume
    G4int fNscin; // number of scintillator strips on the layer

    //hits
    CaloBPCHits::Coll fHits;

};

#endif



















