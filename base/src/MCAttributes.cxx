
//_____________________________________________________________________________
//
// IO for event attributes from MC generator
//
//_____________________________________________________________________________

//Geant
#include "G4Event.hh"

//local classes
#include "MCAttributes.h"

//_____________________________________________________________________________
MCAttributes::MCAttributes(): Detector(), fNam("MCEvent") {

}//MCAttributes

//_____________________________________________________________________________
void MCAttributes::BeginEvent(const G4Event *evt) {

  //generator data

  MCAttribData *dat = dynamic_cast<MCAttribData*>(evt->GetUserInformation());
  if(!dat) return;

  //load the input data
  fDat.LoadGenVal(*dat);

}//BeginEvent

