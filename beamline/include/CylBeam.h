
#ifndef CylBeam_h
#define CylBeam_h

// Cylindrical beampipe section

class GeoParser;

#include "Detector.h"

class CylBeam: public Detector {

  public:

    CylBeam(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; //segment name

};

#endif



