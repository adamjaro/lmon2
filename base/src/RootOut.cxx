
//_____________________________________________________________________________
//
// helper class for ROOT TTree output
//_____________________________________________________________________________

//C++
#include <string>

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TClass.h"
#include "TROOT.h"
#include "TSystem.h"

//Geant
#include "G4GenericMessenger.hh"

//local classes
#include "RootOut.h"

using namespace std;

//_____________________________________________________________________________
RootOut::RootOut(): fOut(0), fDetTree(0) {

  //default name for output file
  fOutputName = "./lmon.root";

  //command for name of output file
  fMsg = new G4GenericMessenger(this, "/lmon/output/");
  fMsg->DeclareProperty("name", fOutputName);

}//RootOut

//_____________________________________________________________________________
void RootOut::Open() {

  std::string nam(fOutputName.data());

  G4cout << "RootOut::Open, " << nam << G4endl;

  //create directory for the output
  if( nam.find_last_of("/") != string::npos ) {
    string dir = nam.substr(0, nam.find_last_of("/"));
    if(!dir.empty()) gSystem->MakeDirectory(dir.c_str());
  }

  //create the output file
  fOut = new TFile(nam.c_str(), "recreate");

  //test if file exists
  if(!fOut->IsOpen()) {
    G4String description = "Can't open output: " + fOutputName;
    G4Exception("DetectorConstruction::CreateOutput", "OutputNotOpen01", FatalException, description);
  }

  //output detector tree
  fDetTree = new TTree("DetectorTree", "DetectorTree");

  //G4cout << "RootOut::Open, DetectorTree created: " << fDetTree->GetCurrentFile()->GetName() << G4endl;

}//Open

//_____________________________________________________________________________
void RootOut::FillTree() {

  fDetTree->Fill();

}//FillTree

//_____________________________________________________________________________
void RootOut::Close() {

  if(fOut) fOut->cd();

  //write the tree
  if(fDetTree) fDetTree->Write(0, TObject::kOverwrite);

  //close the output file
  if(fOut) fOut->Close();

}//Close

















