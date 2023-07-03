
#ifndef ComponentBuilder_h
#define ComponentBuilder_h

// constructs detectors and components

#include "BuilderBase.h"

class ComponentBuilder {

  public:

    ComponentBuilder(G4LogicalVolume *top, GeoParser *geo, std::vector<Detector*> *vdet);

  private:

    std::vector<BuilderBase*> fBuild; //all builders

};

#endif

