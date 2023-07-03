
#ifndef BuilderBase_h
#define BuilderBase_h

//abstract base class for component builders

#include "Detector.h"

class BuilderBase {

  public:

    virtual ~BuilderBase() = default;

    virtual Detector* FindAndLoad(G4String type, G4String name) = 0;

};//BuilderBase

#endif

