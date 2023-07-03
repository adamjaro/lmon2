
//_____________________________________________________________________________
//
// Direct generator data to be put to output tree
//
//_____________________________________________________________________________

//C++
#include <map>
#include <vector>

//ROOT
#include "TTree.h"

//Geant
#include "G4ios.hh"

//local classes
#include "MCAttribData.h"

using namespace std;

//_____________________________________________________________________________
MCAttribData::MCAttribData(): G4VUserEventInformation() {

  //names of generator values
  vector<string> names({
    "true_x", "true_y", "true_Q2", "true_W2",
    "true_el_pT", "true_el_theta", "true_el_phi", "true_el_E", "true_el_Q2",
    "flux_photon_per_s", "power_W",
    "true_phot_theta", "true_phot_phi", "true_phot_E", "num_interactions"
  });

  //set the generator values
  vector<string>::const_iterator inam = names.cbegin();
  for(; inam != names.cend(); inam++) {

    fGenVal.insert( make_pair(*inam, new Double_t(0)) );
  }

}//MCAttribData

//_____________________________________________________________________________
MCAttribData::MCAttribData(const MCAttribData& d): G4VUserEventInformation(),
    fGenVal(d.fGenVal) {

}//MCEvtDat

//_____________________________________________________________________________
void MCAttribData::ConnectInput(TTree *t) {

  //connect the generator values from input tree

  map<string, Double_t*>::const_iterator ival = fGenVal.cbegin();
  for(; ival != fGenVal.cend(); ival++) {

    //branch name, proceed only if present
    string bnam = (*ival).first;
    if( !t->FindBranch(bnam.c_str()) ) continue;

    //set address inside the map
    t->SetBranchAddress(bnam.c_str(), (*ival).second);
  }

}//ConnectInput

//_____________________________________________________________________________
void MCAttribData::LoadGenVal(const MCAttribData& d) {

  //load generator values to local map

  map<string, Double_t*>::const_iterator ival = d.fGenVal.cbegin();
  for(; ival != d.fGenVal.cend(); ival++) {

    *fGenVal[ (*ival).first ] = *(*ival).second;
  }

}//LoadGenVal

//_____________________________________________________________________________
void MCAttribData::SetVal(const std::string& nam, Double_t val) {

  //set value for a given name
  *fGenVal[nam] = val;

}//SetVal

//_____________________________________________________________________________
void MCAttribData::CreateOutput(TTree *t) {

  //make output branches

  map<string, Double_t*>::const_iterator ival = fGenVal.cbegin();
  for(; ival != fGenVal.cend(); ival++) {

    string bnam = (*ival).first;
    string leaf = bnam + "/D";

    t->Branch(bnam.c_str(), (*ival).second, leaf.c_str());

    //G4cout << bnam << " " << leaf << G4endl;
  }

}//CreateOutput

//_____________________________________________________________________________
void MCAttribData::Print(string msg, string dat) {

  G4cout << msg << " " << *fGenVal[dat] << G4endl;

}//Print























