
#ifndef PhotoCathPMT_h
#define PhotoCathPMT_h

// PMT photocathode

#include "Detector.h"
#include "G4VSensitiveDetector.hh"
#include "PhotoHitsV2.h"

class GeoParser;
class MCParticleAction;

class PhotoCathPMT : public Detector, public G4VSensitiveDetector {

  public:

    PhotoCathPMT(const G4String&, GeoParser*, G4LogicalVolume *top=0x0);
    virtual ~PhotoCathPMT() {}

    G4LogicalVolume* CreateGeometry(G4double radius, G4double dz, G4VisAttributes *vmain, G4VisAttributes *vcath);
    G4VPhysicalVolume* GetGlassPhysVol() { return fGlassPhysVol; }

    //called via G4VSensitiveDetector
    virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory*);

    //called via Detector
    virtual const G4String& GetName() const { return fNam; }
    virtual void CreateOutput(TTree*);
    virtual void ClearEvent();
    virtual void FinishEvent();

  private:

    void SetOptics(G4Material *mat, G4LogicalVolume *glass_vol, G4Material *cath_mat, G4LogicalVolume *cath_vol);
    void SetBoundary(G4VPhysicalVolume *cath_phys);

    std::vector<G4double> LambdaNMtoEV(const std::vector<G4double>& lambda);

    G4String fNam; // name of sensitive photocathode logical volume

    G4VPhysicalVolume *fGlassPhysVol; // physical volume for the glass layer

    //PhotoHits fHits; // prototype hit collection
    PhotoHitsV2::Coll fHits; // hit collection

    const MCParticleAction *fStack; // stack for primary particles

};











#endif

