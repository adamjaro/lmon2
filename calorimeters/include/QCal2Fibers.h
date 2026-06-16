
#ifndef QCal2Fibers_h
#define QCal2Fibers_h

#include "Detector.h"
#include "G4VSensitiveDetector.hh"

#include "HitAtID.h"

class QCal2Fibers : public Detector, public G4VSensitiveDetector {

  public:

    QCal2Fibers(const G4String&, GeoParser*, G4LogicalVolume*);

    //called via Detector
    virtual const G4String& GetName() const {return fNam;}
    virtual void Add(std::vector<Detector*> *vec);
    void CreateOutput(TTree *tree) { fHits.CreateOutput(fNam, tree); }
    void ClearEvent() { fHits.ClearEvent(); }
    void FinishEvent() { fHits.FinishEvent(); }

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

  private:

    G4LogicalVolume* MakeCell(GeoParser *geo);
    G4LogicalVolume* MakeFiberYZ(Double_t L, Double_t yL, Double_t r, Double_t ds, const G4String& nam, G4Material *mat, Double_t theta);
    G4LogicalVolume* MakeStraightFib(G4double cladD, G4double coreD, G4double Lz,
      G4Material *clad_mat, G4Material *core_mat, const G4String& clad_nam, const G4String& core_nam, GeoParser *geo);
    //void MakeFiberLguide(std::vector<std::array<G4double, 4>>& pos);

    G4String fNam; // name of detector sensitive logical volume

    Detector *fOpDet; // optical photon detector

    HitAtID::Coll fHits; // hits in fibers as in PWO

    G4double fMaxOptEn; // maximal energy for optical photon, eV
    G4double fMinOptEn; // minimal energy for optical photon, eV

    G4int fUpHist; // hierarchy for hits in fibers

};

#endif


















