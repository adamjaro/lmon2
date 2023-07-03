
//_____________________________________________________________________________
//
// generator reader for TParticle clones array
//
//_____________________________________________________________________________

//C++
#include <string>

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"

//Geant
#include "G4GenericMessenger.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"

//local classes
#include "MCAttribData.h"
#include "TParticleReader.h"

using namespace std;

//_____________________________________________________________________________
TParticleReader::TParticleReader() : G4VPrimaryGenerator(), fIn(0), fMsg(0),
    fTree(0), fPart(0), fIev(0), fDat(0) {

  //default input
  fInputName = "/home/jaroslav/sim/eic-lgen/lgen.root";

  //command for input name
  fMsg = new G4GenericMessenger(this, "/lmon/input/tparticle/");
  fMsg->DeclareProperty("name", fInputName);

  //pdg selection
  fMsg->DeclareMethod("select", &TParticleReader::AddSelectPdg);

}//TParticleReader

//_____________________________________________________________________________
void TParticleReader::GeneratePrimaryVertex(G4Event *evt) {

  //open input on the first call
  if(!fIn) OpenInput();

  if( fIev%100000 == 0 ) {
    G4cout << "TParticleReader::GeneratePrimaryVertex, event number: " << fIev << G4endl;
  }

  //load the event
  Int_t nb = fTree->GetEntry(fIev++);
  if(nb <= 0) {
    G4cout << "TParticleReader::GeneratePrimaryVertex: no more events" << G4endl;
    return;
  }

  //particles loop
  for(Int_t i=0; i<fPart->GetEntriesFast(); i++) {
    TParticle *p = dynamic_cast<TParticle*>(fPart->At(i));

    //select according to pdg if required
    if( fSelPdg.empty() != true && fSelPdg.find(p->GetPdgCode()) == fSelPdg.end() ) continue;

    //primary vertex for the particle
    G4PrimaryVertex *vtx = new G4PrimaryVertex(p->Vx()*mm, p->Vy()*mm, p->Vz()*mm, 0);

    //G4cout << "  particle: " << p->GetPdgCode() << " " << p->Energy() << " " << p->Vz() << G4endl;

    //create the G4 particle
    G4PrimaryParticle *gp = new G4PrimaryParticle(p->GetPdgCode(), p->Px()*GeV, p->Py()*GeV, p->Pz()*GeV, p->Energy()*GeV);
    gp->SetPolarization(0, 0, 0);

    //put the G4 particle to the vertex
    vtx->SetPrimary(gp);

    //put the vertex to the event
    evt->AddPrimaryVertex(vtx);

  }//particles loop

  //fDat->Print("TParticleReader y:", "true_y");

  //add event data if available
  if(fDat) evt->SetUserInformation(new MCAttribData(*fDat));

}//GeneratePrimaryVertex

//_____________________________________________________________________________
void TParticleReader::OpenInput() {

  G4cout << "TParticleReader::OpenInput: " << fInputName << G4endl;

  TFile *cdd = TFile::CurrentFile(); // previously opened output file

  //open the input
  fIn = TFile::Open(fInputName.data(), "read");

  //test if file exists
  if(!fIn) {
    string description = "Can't open input: " + fInputName;
    G4Exception("TParticleReader::OpenInput", "InputNotOpen01", FatalException, description.c_str());
  }

  //connect the clones array
  fTree = dynamic_cast<TTree*>(fIn->Get("ltree"));
  fTree->SetBranchAddress("particles", &fPart);

  //attach the generator data
  fDat = new MCAttribData();
  fDat->ConnectInput(fTree);

  if(cdd) cdd->cd(); // move back to the previous file

}//OpenInput

//_____________________________________________________________________________
void TParticleReader::AddSelectPdg(G4int pdg) {

  //pdg selection

  fSelPdg.insert(pdg);

}//AddSelectPdg

