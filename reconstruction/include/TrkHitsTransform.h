
#ifndef TrkHitsTransform_h
#define TrkHitsTransform_h

// geometry transformation for TrkPlaneBasicHits

#include "TrkPlaneBasicHits.h"
class GeoParser;

class TrkHitsTransform {

  public:

    TrkHitsTransform(std::string nam, TTree *tree);

    void CreateOutput(TTree *tree);

    void LoadV6p3_rev3(GeoParser& geo, std::string section);

    void ProcessEvent();
    void PrintCounts();

    TrkPlaneBasicHits::Coll& GetHits() { return fHits; }

    const std::string& GetName() { return fNam; }

  private:

    TrkPlaneBasicHits::Coll fHits; // input hits

    std::string fNam; // detector name

    bool fWriteOut=false; // flag to write output tree

    //selection criteria
    Double_t fEmin; // keV, threshold in energy

    //position and rotation to transform the hits
    Double_t fX=0; // x position, mm
    Double_t fY=0; // y position, mm
    Double_t fZ=0; // z position, mm
    Double_t fTheta=0; // rotation in z-x plane, rad

    //counters for processed hits
    unsigned long fNall=0;
    unsigned long fNsel=0;

};

#endif

