
#ifndef RunAction_h
#define RunAction_h

// generic run action

#include "G4UserRunAction.hh"

class G4Run;
class DetectorConstruction;

class RunAction : public G4UserRunAction {

  public:

    RunAction();
    virtual ~RunAction() {}

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

  private:

    const DetectorConstruction *fDet; // detector

    clock_t fStart; // start time

};

#endif

