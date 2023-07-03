
#ifndef UcalA290_h
#define UcalA290_h

// NIM A290 (1990) 95-108

#include "Detector.h"
#include "G4VSensitiveDetector.hh"

class GeoParser;
class TH1F;

class UcalA290 : public Detector, public G4VSensitiveDetector {

  public:

    UcalA290(const G4String&, GeoParser*, G4LogicalVolume *top);
    virtual ~UcalA290() {}

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *, G4TouchableHistory*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void ClearEvent();
    virtual void FinishEvent();
    virtual void CreateOutput(TTree *tree);

  private:

    //Birks correction
    G4double BirksCorrectedEnergyDeposit(G4Step *step);

    //geometry utilities
    void MakeAlFront(G4LogicalVolume *modv, G4double al_z);
    G4LogicalVolume *MakeScinLayer(G4Material *mod_mat, G4double scin_z, G4double spacer_z);
    G4LogicalVolume *MakeAbsoLayer(G4double abso_z, G4double clad_z, G4String eh);
    G4LogicalVolume *MakeAbsoNoClad(G4double abso_z, G4String eh);
    void DefineMaterials();

    G4String fNam; // name of detector sensitive logical volume
    GeoParser *fGeo; // geometry parser

    G4double fModXY; // module transverse size, mm

    G4bool fUseBirksCorrection; // use Birks correction to deposited energy
    G4double fBirksCoefficient; // value of Birks coefficient in mm/MeV

    G4double fMaxTime; // maximal time for signal integration

    unsigned fStartHAC1; // first layer in HAC1 section
    unsigned fStartHAC2; // first layer in HAC2 section

    Double_t fEdepEMC; // deposited energy in EMC section
    Double_t fEdepHAC1; // deposited energy in HAC1 section
    Double_t fEdepHAC2; // deposited energy in HAC2 section
    std::vector<Float_t> *fELayer; // deposited energy in individual layers
    TH1F *fTimeStep; // time in step distribution in event, ns















};

#endif

