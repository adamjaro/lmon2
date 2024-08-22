
#ifndef HepMC3RootTreeReader_h
#define HepMC3RootTreeReader_h

// generator reader for HepMC3 ROOT tree

#include "G4VPrimaryGenerator.hh"

namespace HepMC3 {
  class ReaderRootTree;
}

class HepMC3RootTreeReader : public G4VPrimaryGenerator {

  public:

    HepMC3RootTreeReader();
    virtual ~HepMC3RootTreeReader();

    virtual void GeneratePrimaryVertex(G4Event*);

  private:

    void OpenInput();

    G4GenericMessenger *fMsg; // messenger for name of input file
    G4String fInputName; // name of input file

    std::shared_ptr<HepMC3::ReaderRootTree> fRead; // HepMC3 reader

    unsigned long fIev; // index of current event

    std::map<std::string, std::string> fHepmcAttrib; // event attributes

};

#endif


