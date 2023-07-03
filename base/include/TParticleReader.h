
#ifndef TParticleReader_h
#define TParticleReader_h

// generator reader for TParticle clones array

#include <set>

#include "Rtypes.h"

#include "G4VPrimaryGenerator.hh"
#include "G4VUserEventInformation.hh"

class G4Event;
class G4GenericMessenger;

class TFile;
class TTree;
class TClonesArray;
class MCAttribData;

class TParticleReader : public G4VPrimaryGenerator {

  public:

    TParticleReader();
    virtual ~TParticleReader() {}

    virtual void GeneratePrimaryVertex(G4Event*);

  private:

    void OpenInput(); // open the input file
    void AddSelectPdg(G4int); // add pdg selected to generate

    TFile *fIn; // input from ROOT
    G4String fInputName; // name of input file

    G4GenericMessenger *fMsg; // messenger for name of input file

    TTree *fTree; // input tree
    TClonesArray *fPart; // particles clones array
    Long64_t fIev; // index of current event

    MCAttribData *fDat; //event data

    std::set<Int_t> fSelPdg; // codes for particles selected to generate

};

#endif

