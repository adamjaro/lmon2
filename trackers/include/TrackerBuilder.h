
#ifndef TrackerBuilder_h
#define TrackerBuilder_h

//builder for tracker detectors

#include <map>

#include "BuilderBase.h"

class TrackerBuilder : public BuilderBase {

  public:

    TrackerBuilder(G4LogicalVolume *top, GeoParser *geo);

    Detector* FindAndLoad(G4String type, G4String name);

  private:

    G4LogicalVolume *fTop; // top world volume

    GeoParser *fGeo; // geometry parser

    //factory function for individual detectors
    template<class det> Detector* MakeDet(G4String nam, GeoParser *geo, G4LogicalVolume *vol) {
      return new det(nam, geo, vol);
    }
    typedef Detector* (TrackerBuilder::*MakeDetPtr)(G4String, GeoParser*, G4LogicalVolume*);

    std::map<G4String, MakeDetPtr> fDets; // local defined detectors

};

#endif
