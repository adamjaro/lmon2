
#ifndef WScFiZXv3_h
#define WScFiZXv3_h

// W/ScFi by Zhiwan Xu (UCLA), v3

#include "Detector.h"
#include "G4VSensitiveDetector.hh"

class GeoParser;

class WScFiZXv3 : public Detector, public G4VSensitiveDetector {

  public:

    WScFiZXv3(const G4String&, GeoParser*, G4LogicalVolume *top);
    virtual ~WScFiZXv3() {}

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void ClearEvent();
    virtual void CreateOutput(TTree *tree);

  private:

    G4LogicalVolume* MakeTower(G4double calorSizeXY, G4double calorEMZ);

    G4double BirksCorrectedEnergyDeposit(G4Step *step);

    G4String fNam; // name of detector sensitive logical volume
    G4double fBirksCoefficient; // value of Birks coefficient in mm/MeV

    G4bool fCheckOverlaps; // option to activate checking of volumes overlaps

    Double_t fEdep; // deposited energy, GeV

};

#endif

