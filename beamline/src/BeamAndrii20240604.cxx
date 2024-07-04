
//_____________________________________________________________________________
//
// Prototype beam pipe, code example by Andrii Natochii, June 4, 2024
//
//_____________________________________________________________________________

//ROOT

//Geant
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Polycone.hh"
#include "G4SystemOfUnits.hh"
#include "G4Torus.hh"
#include "G4UnionSolid.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "G4NistManager.hh"


//local classes
#include "BeamAndrii20240604.h"
#include "GeoParser.h"

//_____________________________________________________________________________
BeamAndrii20240604::BeamAndrii20240604(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "BeamAndrii20240604: " << fNam << G4endl;

  G4LogicalVolume *world_log = top;

  G4Material *mat_vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");


	G4bool addAbsorber = true; // build SR photon absorber
	G4bool visibility = true; // components visibility


	G4Box* beampipe_cut = new G4Box("beampipe_cut",2*m, 2*m, 10*m);


	//---------------------------------------------------------------------------------------------------------------//
	//- Build an extended bwd beam pipe
	const G4double B2eR_bend = 0.020 * rad;
	const G4double B2eR_len = 5.5 * m;
	const G4double B2eR_cut_offset = 1 * cm;
	const G4double B2eR_beampipe_radius = 5 * cm;
	G4ThreeVector B2eR_start(0*cm,0*cm,1675*cm);
	G4double torus_radius = (0.5*B2eR_len/sin(B2eR_bend/2.)); 
	G4ThreeVector B2eR_end = B2eR_start - G4ThreeVector(torus_radius,0,B2eR_start.z());	
	B2eR_end.rotateY(B2eR_bend);
	B2eR_end += G4ThreeVector(torus_radius,0,B2eR_start.z());  

	const G4int numZPlanes = 8;
	const G4double zPlane[numZPlanes] = {
		9280.0 * mm, 9288.0 * mm, 9364.0 * mm, 11364.0 * mm, 11414.0 * mm, 
		15000.0 * mm, 16000.0 * mm, 25000.0 * mm};
	const G4double rInner[numZPlanes] = {0,0,0,0,0,0,0,0};
	const G4double rOuter[numZPlanes] = {
		128.6/2 * mm, 128.6/2 * mm, 154.0/2 * mm, 154.0/2 * mm, 196.0/2 * mm, 
		196.0/2 * mm, B2eR_beampipe_radius, B2eR_beampipe_radius};
	G4Polycone* polycone = new G4Polycone("polycone",
		0, 2.0 * M_PI, numZPlanes, zPlane, rInner, rOuter);

	// Torus - bent beam pipe
	G4RotationMatrix* torus_rot = new G4RotationMatrix();
	torus_rot->rotateX(-M_PI * 0.5 * rad);
	torus_rot->rotateY(M_PI * rad);
	
	G4VSolid* torus_solid = new G4Torus("torus_solid",0,B2eR_beampipe_radius,
		torus_radius,0,B2eR_bend);
	G4VSolid* bwd_beampipeVac_solid = new G4UnionSolid("bwd_beampipeVac_solid",
		polycone,torus_solid,
		torus_rot,G4ThreeVector(torus_radius, 0, B2eR_start.z()));

	//---------------------------------------------------------------------------------------------------------------//
	// B2eR beam pipe extension
	G4RotationMatrix* torus_ext_rot = new G4RotationMatrix();
	torus_ext_rot->rotateY(-B2eR_bend);

	G4double torus_ext_len = 10 * cm;
	G4VSolid* torus_ext_solid = new G4Tubs("torus_ext_solid",0,B2eR_beampipe_radius, torus_ext_len/2.0, 0, 2.0 * M_PI * rad);
	bwd_beampipeVac_solid = new G4UnionSolid("bwd_beampipeVac_solid",
		bwd_beampipeVac_solid,torus_ext_solid,
		torus_ext_rot,G4ThreeVector(	B2eR_end.x() + 0.5*torus_ext_len*sin(B2eR_bend),
						B2eR_end.y(),
						B2eR_end.z() + 0.5*torus_ext_len*cos(B2eR_bend)));

	// Cut the BWD beam pipe at the end of B2eR
	bwd_beampipeVac_solid = new G4SubtractionSolid("bwd_beampipeVac_solid",
		bwd_beampipeVac_solid,beampipe_cut,0,G4ThreeVector(0,0,B2eR_end.z() + B2eR_cut_offset + 10*m));
	//---------------------------------------------------------------------------------------------------------------//
	// FarBWD beam pipe
	G4double farbwd_len = 20 * m;
	G4double farbwd_ext = 5 * m;
	G4double lumi_halflen = 2.25 * m;
	G4RotationMatrix* farbwd_rot = new G4RotationMatrix();
	farbwd_rot->rotateY(-B2eR_bend);
	G4RotationMatrix* lumi_rot = new G4RotationMatrix();
	lumi_rot->rotateY(B2eR_bend);
	G4RotationMatrix* cut_rot = new G4RotationMatrix();
	cut_rot->rotateY(B2eR_bend);

	// main beam pipe
	G4VSolid* farbwd_tube = new G4Tubs("farbwd_tube", 0,B2eR_beampipe_radius, (farbwd_len+farbwd_ext)/2., 0, 2.0 * M_PI * rad);
	// luminosity beam pipe
	G4VSolid* lumi_tube = new G4Tubs("farbwd_tube", 0,B2eR_beampipe_radius, lumi_halflen, 0, 2.0 * M_PI * rad);
	// ante-chamber
	G4Box* ante_chmbr_box = new G4Box("ante_chmbr_box",22*cm/2.,5*cm/2.,farbwd_len/2.);
	// unite main and ante-chamber
	farbwd_tube = new G4UnionSolid("farbwd_tube",
		farbwd_tube,ante_chmbr_box,0,G4ThreeVector(22*cm/2.,0,0));
	// distance between the main beam pipe center (placement coordinates) and the cut
	G4double L = (B2eR_end.z() + 0.5*farbwd_len*cos(B2eR_bend)-(B2eR_end.z() + B2eR_cut_offset))/cos(B2eR_bend);
	// half-width of the cut in the main beam pipe
	G4double M = B2eR_end.x() + 0.5*farbwd_len*sin(B2eR_bend) - L*sin(B2eR_bend);
	// displacement of the lumi pipe w.r.t. the main beam pipe
	G4double dL = M * sin(B2eR_bend); // along main beam pipe Z-axis
	G4double dK = M * cos(B2eR_bend); // along main beam pipe X-axis
	// unite lumi and main+ante-chamber in the main beam pipe coordinate system
	farbwd_tube = new G4UnionSolid("farbwd_tube",farbwd_tube,lumi_tube,lumi_rot,G4ThreeVector(-dK,0,-L-dL));

	// Cut the FarBWD beam pipe at the end of B2eR
		G4VSolid* farbwd_beampipeVac_solid = new G4SubtractionSolid("farbwd_beampipeVac_solid",
		farbwd_tube,beampipe_cut,cut_rot,G4ThreeVector(0,0,-(L+(10*m)/cos(B2eR_bend))));

	// Place the assembly
	G4LogicalVolume* farbwd_beampipeVac_log = new G4LogicalVolume(farbwd_beampipeVac_solid,mat_vac,"farbwd_beampipeVac_log");
	//farbwd_beampipeVac_log->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,simPar->GetTrackTcut())); 
	new G4PVPlacement(farbwd_rot, 
		G4ThreeVector(	B2eR_end.x() + 0.5*farbwd_len*sin(B2eR_bend),
				B2eR_end.y(),
				B2eR_end.z() + 0.5*farbwd_len*cos(B2eR_bend)),
		farbwd_beampipeVac_log,"vac_farbwd_beampipeVac",world_log,false,0,true);

	// Visualization
	G4VisAttributes* farbwd_beampipeVac_vis = new G4VisAttributes(G4Color::Gray());
	farbwd_beampipeVac_vis->SetVisibility(visibility);
	farbwd_beampipeVac_log->SetVisAttributes(farbwd_beampipeVac_vis);
	//---------------------------------------------------------------------------------------------------------------//
	// Place an absorber inside the FarBWD vacuum
	if(addAbsorber && farbwd_beampipeVac_log)
	{
		G4Tubs* abs_tube = new G4Tubs("abs_tube_farbwd",0,B2eR_beampipe_radius,1*cm,0, 2.0 * M_PI * rad);
		G4LogicalVolume* abs_log = new G4LogicalVolume(abs_tube,mat_vac, "abs_log_farbwd");
		new G4PVPlacement(0,G4ThreeVector(0,0,(farbwd_len+farbwd_ext)/2. - 1 * cm),
			abs_log,"abs_farbwd",farbwd_beampipeVac_log,false,0,true);

		//- Visualization
		G4VisAttributes* abs_vis = new G4VisAttributes(G4Color::Magenta());
		abs_vis->SetVisibility(visibility);
		abs_log->SetVisAttributes(abs_vis);
	}

	//---------------------------------------------------------------------------------------------------------------//
	// Place a lumi window
	G4double lumiWindStart_posz = B2eR_end.z() + B2eR_cut_offset + lumi_halflen;
	G4cout<<"[INFO] lumiWindStart_posz = "<<lumiWindStart_posz/cm<<" [cm]"<<G4endl;
	G4Box* lumiWind_box = new G4Box("lumiWind_box",B2eR_beampipe_radius,B2eR_beampipe_radius,1 * cm);
	G4LogicalVolume* lumiWind_log = new G4LogicalVolume(lumiWind_box,mat_vac,"lumiWind_log");
	new G4PVPlacement(0,G4ThreeVector(0,0,lumiWindStart_posz+ 1*cm),lumiWind_log,"lumiWind",world_log,false,0,true);



}//BeamAndrii20240604

















