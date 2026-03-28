
#ifndef QCal2_h
#define QCal2_h

#include "Detector.h"
#include "G4VSensitiveDetector.hh"

#include "CalPWOHits.h"

class QCal2 : public Detector, public G4VSensitiveDetector {

  public:

    QCal2(const G4String&, GeoParser*, G4LogicalVolume*);

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

    G4String fNam; // name of detector sensitive logical volume

    Detector *fOpDet; // optical photon detector

    CalPWOHits::Coll fHits; // hits in fibers as in PWO

};

#endif

