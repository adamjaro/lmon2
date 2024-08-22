
//_____________________________________________________________________________
//
// generator action
//
//_____________________________________________________________________________

//C++

//Geant
#include "G4GenericMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"

//local classes
#include "GeneratorAction.h"
//#include "TxReader.h"
//#include "Pythia6Reader.h"
#include "TParticleReader.h"
#include "HEPEvtInterface.h"
#include "use_hepmc3.h"
#ifdef USE_HEPMC3
  #include "HepMC3AsciiReader.h"
  #include "HepMC3RootTreeReader.h"
#endif

using namespace std;

//_____________________________________________________________________________
GeneratorAction::GeneratorAction() : G4VUserPrimaryGeneratorAction(), fGenType("gun"), fGen(nullptr) {

  //command for generator reader type
  fMsg = new G4GenericMessenger(this, "/lmon/input/");
  fMsg->DeclareProperty("type", fGenType);

  //prepare the generators
  //fGenAll.insert( make_pair("tx", new TxReader()) );
  //fGenAll.insert( make_pair("pythia6", new Pythia6Reader()) );
  fGenAll.insert( make_pair("tparticle", make_unique<TParticleReader>()) );
  fGenAll.insert( make_pair("gun", make_unique<G4ParticleGun>()) );
  fGenAll.insert( make_pair("gps", make_unique<G4GeneralParticleSource>()) );
  fGenAll.insert( make_pair("hepevt", make_unique<HEPEvtInterface>()) );
  #ifdef USE_HEPMC3
    fGenAll.insert( make_pair("hepmc_ascii", make_unique<HepMC3AsciiReader>()) );
    fGenAll.insert( make_pair("hepmc_root_tree", make_unique<HepMC3RootTreeReader>()) );
  #endif

}//GenReader

//_____________________________________________________________________________
GeneratorAction::~GeneratorAction() {

  delete fMsg;

}//~GenReader

//_____________________________________________________________________________
void GeneratorAction::GeneratePrimaries(G4Event *evt) {

  //G4cout << "GeneratorAction::GeneratePrimaries " << fGenType << G4endl;

  //select the reader at the first call
  if(!fGen) {
    map<G4String, unique_ptr<G4VPrimaryGenerator>>::iterator igen = fGenAll.find(fGenType);
    if (igen == fGenAll.end()) return;

    fGen = move( (*igen).second );
  }

  //generate the event
  fGen->GeneratePrimaryVertex(evt);

}//GeneratePrimaries


















