
#ifndef ActionInitialization_h
#define ActionInitialization_h

//standard action initialization

#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization {

  public:

    ActionInitialization() : G4VUserActionInitialization() {}
    virtual ~ActionInitialization() {}

    virtual void Build() const;

};

#endif

