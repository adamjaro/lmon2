
#ifndef DetectorConstruction_h
#define DetectorConstruction_h

// main detector construction

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4Step;
class G4Event;
class G4GenericMessenger;

#include "Detector.h"

class RootOut;
class vector;
class MCAttributes;

class DetectorConstruction : public G4VUserDetectorConstruction {

  public:

    DetectorConstruction();
    virtual ~DetectorConstruction();

    G4VPhysicalVolume* Construct();

    void BeginEvent(const G4Event *evt) const;
    void FinishEvent() const;

    void CreateOutput() const;

    void ConstructSDandField();

  private:

    std::vector<Detector*> *fDet; //all detectors

    RootOut *fOut; // output to ROOT TTree

    MCAttributes *fMCattr; // generated particles

    G4String fGeoName; // name of geometry input

    G4GenericMessenger *fMsg; // messenger for geometry input

};

#endif













