
#ifndef MCAttributes_h
#define MCAttributes_h

// IO for event attributes from MC generator

#include "Detector.h"

#include "MCAttribData.h"

class MCAttributes: public Detector {

  public:

    MCAttributes();

    void BeginEvent(const G4Event *evt);

    //Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void CreateOutput(TTree *tree) { fDat.CreateOutput(tree); }

  private:

    G4String fNam; // attributes name in detectors container

    MCAttribData fDat; //event data

};

#endif

