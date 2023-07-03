
#ifndef CalPWO_h
#define CalPWO_h

// PbWO4 calorimeter

#include "Detector.h"
#include "G4VSensitiveDetector.hh"

#include "CalPWOHits.h"

class GeoParser;
class MCParticleAction;

class CalPWO : public Detector, public G4VSensitiveDetector {

  public:

    CalPWO(const G4String&, GeoParser*, G4LogicalVolume *top);
    virtual ~CalPWO() {}

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory*);

    //called via Detector
    virtual const G4String& GetName() const { return fNam; }
    virtual void Add(std::vector<Detector*> *vec);
    void CreateOutput(TTree*);
    void ClearEvent();
    void FinishEvent();

  private:

    void SetCrystalOptics(G4Material *mat);
    void SetCrystalSurface(G4LogicalVolume *vol);
    void SetCrystalBoundary(G4VPhysicalVolume *crystal, G4VPhysicalVolume *glass);

    std::vector<G4double> LambdaNMtoEV(const std::vector<G4double>& lambda);

    G4LogicalVolume* GetMotherVolume(G4String mother_nam, G4LogicalVolume *top);

    G4String fNam; // name of detector sensitive logical volume

    Detector *fPMT; // PMT photocathode

    const MCParticleAction *fStack; // stack for primary particles

    CalPWOHits::Coll fHits; // hit collection

};

#endif

