
#ifndef HcalA262_h
#define HcalA262_h

// NIM A262 (1987) 229-242

#include "Detector.h"
#include "G4VSensitiveDetector.hh"

class GeoParser;

class HcalA262 : public Detector, public G4VSensitiveDetector {

  public:

    HcalA262(const G4String&, GeoParser*, G4LogicalVolume *top);
    virtual ~HcalA262() {}

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void ClearEvent();
    //virtual void FinishEvent();
    virtual void CreateOutput(TTree *tree);

  private:

    G4double BirksCorrectedEnergyDeposit(G4Step *step);

    G4String fNam; // name of detector sensitive logical volume

    G4int fNem; // number of EM layers to separate the sections

    Double_t fEdepEM; // deposited energy in EM section
    Double_t fEdepHAD; // deposited energy in HAD section
    std::vector<Float_t> *fELayer; // deposited energy in individual layers

    G4bool fUseBirksCorrection; // use Birks correction to deposited energy
    G4double fBirksCoefficient; // value of Birks coefficient in mm/MeV














};

#endif

