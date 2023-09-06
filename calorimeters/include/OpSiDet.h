
#ifndef OpSiDet_h
#define OpSiDet_h

// Simple optical Si detector imitating a PIN photodiode

#include "Detector.h"
#include "G4VSensitiveDetector.hh"
#include "PhotoHitsV2.h"

class OpSiDet : public Detector, public G4VSensitiveDetector {

  public:

    OpSiDet(const G4String&, GeoParser *geo=0x0, G4LogicalVolume *top=0x0);

    G4LogicalVolume* CreateGeometry(G4double dx, G4double dy, G4double dz, G4VisAttributes *vis);
    void MakeBoundary(G4VPhysicalVolume *src_phys, G4VPhysicalVolume *opdet_phys);

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

    //called via Detector
    virtual const G4String& GetName() const { return fNam; }

  private:

    std::vector<G4double> LambdaNMtoEV(const std::vector<G4double>& lambda);

    G4String fNam; // name of sensitive logical volume



};

#endif

