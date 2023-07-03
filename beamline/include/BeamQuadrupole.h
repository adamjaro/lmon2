
#ifndef BeamQuadrupole_h
#define BeamQuadrupole_h

//beam quadrupole magnet

#include "Detector.h"

class BeamQuadrupole : public Detector {

  public:

    BeamQuadrupole(G4String nam, GeoParser *geo, G4LogicalVolume*);
    virtual ~BeamQuadrupole() {}

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

    G4String fNam; // compoment name

    void PrintField(G4FieldManager *fman);

    G4VisAttributes* ColorDecoder(GeoParser *geo);

};

#endif

