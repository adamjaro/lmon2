
#ifndef BeamAndrii20240604_h
#define BeamAndrii20240604_h

// Prototype beam pipe, , code example by Andrii Natochii, June 4, 2024

class GeoParser;

#include "Detector.h"

class BeamAndrii20240604: public Detector {

  public:

    BeamAndrii20240604(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; //segment name

};

#endif

