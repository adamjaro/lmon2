
#ifndef CylSegment_h
#define CylSegment_h

// Cylindrical construction segment

class GeoParser;

#include "Detector.h"

class CylSegment: public Detector {

  public:

    CylSegment(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; //segment name

};

#endif



