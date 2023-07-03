
#ifndef ConeBeam_h
#define ConeBeam_h

// Conical beampipe section

class GeoParser;

#include "Detector.h"

class ConeBeam: public Detector {

  public:

    ConeBeam(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; //segment name

};

#endif

