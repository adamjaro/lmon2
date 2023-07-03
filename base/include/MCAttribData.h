
#ifndef MCAttribData_h
#define MCAttribData_h

//MC attribute data

class TTree;
#include "Rtypes.h"

#include "G4VUserEventInformation.hh"

//_____________________________________________________________________________
class MCAttribData: public G4VUserEventInformation {
  public:

    MCAttribData();
    MCAttribData(const MCAttribData& d);

    void ConnectInput(TTree *t);
    void LoadGenVal(const MCAttribData& d);
    void CreateOutput(TTree *t);

    void SetVal(const std::string& nam, Double_t val);

    void Print() const {} // reimplemented
    void Print(std::string msg, std::string dat);

  private:

    std::map<std::string, Double_t*> fGenVal;

};

#endif

