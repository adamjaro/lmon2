
//_____________________________________________________________________________
//
// standard action initialization,
// selects the event generator
//_____________________________________________________________________________

//local classes
#include "ActionInitialization.h"
#include "GeneratorAction.h"
#include "EventAction.h"
//#include "LgenReader.h"
#include "RunAction.h"
//#include "UniformGen.h"
//#include "GenReader.h"
#include "MCParticleAction.h"

//_____________________________________________________________________________
void ActionInitialization::Build() const {

  //select the generator
  SetUserAction(new GeneratorAction);
  //SetUserAction(new UniformGen);
  //SetUserAction(new GenReader);

  SetUserAction(new MCParticleAction);
  SetUserAction(new EventAction);
  SetUserAction(new RunAction);

}//build

