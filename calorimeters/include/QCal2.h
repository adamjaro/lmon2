
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

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*) {return true;}

  private:

    G4String fNam; // name of detector sensitive logical volume

};

#endif

