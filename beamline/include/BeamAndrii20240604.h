
#ifndef BeamAndrii20240604_h
#define BeamAndrii20240604_h

// Prototype beam pipe, , code example by Andrii Natochii, June 4, 2024

class GeoParser;
class BeamAndrii_SimParameters;
class G4MagneticField;

#include "Detector.h"

class BeamAndrii20240604: public Detector {

  public:

    BeamAndrii20240604(const G4String& nam, GeoParser *geo, G4LogicalVolume *top);

    //Detector
    virtual const G4String& GetName() const {return fNam;}

  private:

  	void BuildParticles();
  	G4MagneticField* GetMagField(G4int,G4LogicalVolume*);
    void CreateBeamPipeSolid();
    void BuildMagnets();

    G4String fNam; //segment name

    std::unique_ptr<BeamAndrii_SimParameters> simPar;

	G4Material* mat_vac; // vacuum material
	G4VSolid* beampipeVac_solid; // beam pipe SV = beampipeUni_solid - (all magnets + IP beam pipe)
	G4VSolid* beampipeUni_solid; // beam pipe SV = original beam pipe + extension
	G4LogicalVolume* world_log; // world LV
	G4VSolid* bwd_beampipeVac_solid; // backward SV
	G4VSolid* farbwd_beampipeVac_solid; // far-backward SV
	G4VSolid* ip_box; // IP beam pipe virtual volume
	G4LogicalVolume* ip_log; // IP beam pipe LV
	G4bool visibility; // components visibility
	G4bool addAbsorber; // build SR photon absorber

	G4double mass;
	G4double gamma;

};

#endif

