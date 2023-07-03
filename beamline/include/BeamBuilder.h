
#ifndef BeamBuilder_h
#define BeamBuilder_h

//builder for a set of beamline components

#include <map>

#include "BuilderBase.h"

class BeamBuilder : public BuilderBase {

  public:

    BeamBuilder(G4LogicalVolume *top, GeoParser *geo);

    Detector* FindAndLoad(G4String type, G4String name);

  private:

    G4LogicalVolume *fTop; // top world volume

    GeoParser *fGeo; // geometry parser

    //factory function for individual components
    template<class det> Detector* MakeDet(G4String nam, GeoParser *geo, G4LogicalVolume *vol) {
      return new det(nam, geo, vol);
    }
    typedef Detector* (BeamBuilder::*MakeDetPtr)(G4String, GeoParser*, G4LogicalVolume*);

    std::map<G4String, MakeDetPtr> fComp; // component definitions
    std::map<G4String, MakeDetPtr> fDets; // local defined detectors

};

#endif
