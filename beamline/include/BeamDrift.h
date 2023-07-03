
#ifndef BeamDrift_h
#define BeamDrift_h

// Beam drift section between B2 and Q3 magnets

class GeoParser;
class G4GenericTrap;

#include "Detector.h"

class BeamDrift: public Detector {

  public:

    BeamDrift(G4String nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; //component name

    G4GenericTrap* MakeGT(
      G4double z0T, G4double x0T, G4double z0B, G4double x0B, G4double z1T, G4double x1T, G4double z1B, G4double x1B,
      G4double ysiz, G4String nam);

};

#endif

