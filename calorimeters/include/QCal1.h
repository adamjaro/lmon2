
#ifndef QCal1_h
#define QCal1_h

#include "Detector.h"
#include "G4VSensitiveDetector.hh"

#include "CalPWOHits.h"

class QCal1 : public Detector, public G4VSensitiveDetector {

  public:

    QCal1(const G4String&, GeoParser*, G4LogicalVolume*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    void CreateOutput(TTree*);
    void ClearEvent();
    void FinishEvent();
    virtual void Add(std::vector<Detector*> *vec);

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

  private:

    G4LogicalVolume* MakeCell(GeoParser *geo);

    G4String fNam; // name of detector sensitive logical volume

    CalPWOHits::Coll fHits; // hits in fibers as in PWO

    Detector *fOpDet; // optical photon detector

};

#endif

