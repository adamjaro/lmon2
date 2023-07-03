
#ifndef GeneratorAction_h
#define GeneratorAction_h

// generator reader

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"

#include <map>

class G4VPrimaryGenerator;
class G4GenericMessenger;

class GeneratorAction : public G4VUserPrimaryGeneratorAction {

  public:

    GeneratorAction();
    virtual ~GeneratorAction();

    virtual void GeneratePrimaries(G4Event*);

  private:

    G4String fGenType; // generator reader type

    G4GenericMessenger *fMsg; // messenger for generator reader type

    G4VPrimaryGenerator *fGen; // generator reader

    std::map<G4String, G4VPrimaryGenerator*> fGenAll; // all generators

};

#endif

