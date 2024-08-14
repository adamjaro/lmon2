
#ifndef FiberXYCal_h
#define FiberXYCal_h

#include "Detector.h"
#include "G4VSensitiveDetector.hh"

#include "CalPWOHits.h"

class FiberXYCal : public Detector, public G4VSensitiveDetector {

  public:

    FiberXYCal(const G4String&, GeoParser*, G4LogicalVolume*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void Add(std::vector<Detector*> *vec);
    void CreateOutput(TTree *tree) { fHits.CreateOutput(fNam, tree); }
    void ClearEvent() { fHits.ClearEvent(); }
    void FinishEvent() { fHits.FinishEvent(); }

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

  private:

    G4LogicalVolume* MakeCell(GeoParser *geo);

    void BuildCellMaterials();

    G4LogicalVolume* GetMotherVolume(G4String mother_nam, G4LogicalVolume *top);

    G4String fNam; // name of detector sensitive logical volume

    CalPWOHits::Coll fHits; // energy deposition in fibers

    Detector *fOpDet; // optical photon detector

};

#endif

