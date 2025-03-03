#ifndef ReadGDML_h
#define ReadGDML_h

// reader for GDML geometry

#include "Detector.h"

class ReadGDML: public Detector {

  public:

    ReadGDML(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; //segment name

};

#endif


