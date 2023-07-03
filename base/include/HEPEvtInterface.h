
#ifndef HEPEvtInterface_h
#define HEPEvtInterface_h

// reader for HEPEvt, based on G4HEPEvtInterfaceIons

#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4HEPEvtParticle.hh"

#include <fstream>
#include <vector>

class G4Event;
class G4GenericMessenger;

class HEPEvtInterface : public G4VPrimaryGenerator {

  public:

    HEPEvtInterface();
    virtual ~HEPEvtInterface() {}

    virtual void GeneratePrimaryVertex(G4Event *evt);

  private:

    void OpenInput(); // function to open the input

    std::ifstream fIn; // HEPEvt input
    G4String fInputName; // name of input file

    G4GenericMessenger *fMsg; // messenger for name of input file

    G4int vLevel; // verbose level

    std::vector<G4HEPEvtParticle*> fHPlist; // input particles

    int fIev; // event number for progress printout

};

#endif

