
#ifndef VacDrift_h
#define VacDrift_h

// Vacuum drift section

class GeoParser;
class G4GenericTrap;

#include "Detector.h"

class VacDrift: public Detector {

  public:

    VacDrift(G4String nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; //component name

    G4GenericTrap* MakeGT(
      G4double z0T, G4double x0T, G4double z0B, G4double x0B, G4double z1T, G4double x1T, G4double z1B, G4double x1B,
      G4double ysiz, G4String nam);

    G4LogicalVolume* GetMotherVolume(G4String mother_nam, G4LogicalVolume *top);

};

#endif

