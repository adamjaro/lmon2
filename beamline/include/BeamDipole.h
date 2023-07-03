
#ifndef BeamDipole_h
#define BeamDipole_h

// beamline dipole magnet

#include "Detector.h"

class GeoParser;

class BeamDipole : public Detector {

  public:

    BeamDipole(G4String nam, GeoParser *geo, G4LogicalVolume*);
    virtual ~BeamDipole() {}

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; // magnet name

    void PrintField(G4FieldManager *fman);

    G4VisAttributes* ColorDecoder(GeoParser *geo);

};

#endif

