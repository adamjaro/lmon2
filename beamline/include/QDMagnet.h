
#ifndef QDMagnet_h
#define QDMagnet_h

// Universal quadrupole - dipole beam magnet

class G4MagneticField;

#include "Detector.h"

class QDMagnet: public Detector {

  public:

    QDMagnet(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4MagneticField* MakeDipoleField(G4double length, G4double angle);
    G4MagneticField* MakeQuadrupoleField(G4double length, G4double K1L);

    G4String fNam; //segment name

    G4double fGamma;
    G4double fMass;

};

#endif

