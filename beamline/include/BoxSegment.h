
#ifndef BoxSegment_h
#define BoxSegment_h

// Construction box

class GeoParser;

#include "Detector.h"

class BoxSegment: public Detector {

  public:

    BoxSegment(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4LogicalVolume* GetMotherVolume(G4String mother_nam, G4LogicalVolume *top);

    G4String fNam; //segment name

};

#endif

