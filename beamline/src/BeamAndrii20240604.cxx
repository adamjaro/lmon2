
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
#include "G4IntersectionSolid.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4QuadrupoleMagField.hh"
//#include "SolenoidMagField.hh"
#include <G4ParticleTable.hh>

//local classes
#include "BeamAndrii20240604.h"
#include "BeamAndrii_SimParameters.h"
#include "GeoParser.h"

//_____________________________________________________________________________
BeamAndrii20240604::BeamAndrii20240604(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "BeamAndrii20240604: " << fNam << G4endl;

	simPar = make_unique<BeamAndrii_SimParameters>();
	simPar->InitDefault();
	simPar->ReadXML( geo->GetS(fNam, "xmlFileName") );
	simPar->PrintBeampipeFileName();
	simPar->PrintIpBeampipeParameters();
	simPar->PrintGeoParameters();
	simPar->PrintMagParameters();
	simPar->PrintAbsParameters();

  //cout << geo->GetS(fNam, "xmlFileName") << endl;

  mat_vac = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  world_log = top;

	addAbsorber = true;
	visibility = true;

  CreateBeamPipeSolid();
  BuildMagnets();

}//BeamAndrii20240604

//_____________________________________________________________________________
void BeamAndrii20240604::CreateBeamPipeSolid() {

  G4cout << "BeamAndrii20240604::CreateBeamPipeSolid" << G4endl;

	if(!mat_vac)
	{
		throw runtime_error(
		"[ERROR] DetectorConstruction::CreateBeamPipeSolid ==> Materials are not created\n");
	}
/*
	//- Original (beam pipe)
	auto mesh = CADMesh::TessellatedMesh::FromOBJ(simPar->GetBeampipeFileName().c_str());
	G4VSolid* beampipeOrig_solid = mesh->GetSolid();

	//- Extention
	G4VSolid* beampipeExtFwd_solid = new G4Tubs("beampipeExtFwd_solid", 0, 
		simPar->GetExtR(), simPar->GetExtL()/2, 0, 2.0 * M_PI);

	//- Union (beam pipe with an extention)
	beampipeUni_solid = new G4UnionSolid("beampipeUni_solid", 
		beampipeOrig_solid, beampipeExtFwd_solid, 0, 
		G4ThreeVector(0, 0, simPar->GetExtE()-simPar->GetExtL()/2));
*/

	//- Extention only (JA)
  beampipeUni_solid = new G4Tubs("beampipeExtFwd_solid", 0, 
		simPar->GetExtR(), simPar->GetExtL()/2, 0, 2.0 * M_PI);

	G4Box* beampipe_cut = new G4Box("beampipe_cut",2*m, 2*m, 10*m);

	//- Cut the vacuum on the BWD side 
	beampipeUni_solid = new G4SubtractionSolid("beampipeUni_solid",
		beampipeUni_solid,beampipe_cut,0,G4ThreeVector(0,0,9280.0 * mm + 10*m));

	//- Copy
	beampipeVac_solid = beampipeUni_solid;

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
	bwd_beampipeVac_solid = new G4UnionSolid("bwd_beampipeVac_solid",
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
	farbwd_beampipeVac_solid = new G4SubtractionSolid("farbwd_beampipeVac_solid",
		farbwd_tube,beampipe_cut,cut_rot,G4ThreeVector(0,0,-(L+(10*m)/cos(B2eR_bend))));

	// Place the assembly
	G4LogicalVolume* farbwd_beampipeVac_log = new G4LogicalVolume(farbwd_beampipeVac_solid,mat_vac,"farbwd_beampipeVac_log");
	farbwd_beampipeVac_log->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,simPar->GetTrackTcut())); 

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

}//CreateBeamPipeSolid

//_____________________________________________________________________________
void BeamAndrii20240604::BuildMagnets() {

  G4cout << "BeamAndrii20240604::BuildMagnets ";
  G4cout << mat_vac <<" "<< beampipeVac_solid <<" "<< beampipeUni_solid <<" "<< world_log <<" "<< bwd_beampipeVac_solid << G4endl;

  //return;

	if(!mat_vac || !beampipeVac_solid || !beampipeUni_solid || !world_log || !bwd_beampipeVac_solid)
	{
		throw runtime_error(
		"[ERROR] DetectorConstruction::BuildMagnets ==> Materials, beam pipe solids, or world LV are not created\n");
	}
	// Particles
	BuildParticles();

	//- Magnets (fill the beam pipe volume with the magnetic field)
	for(int iMag = 0; iMag < simPar->GetMagNum(); iMag++)
	{
		// magnet virtual volume
		G4VSolid* mag_vol;
		// magnet center coordinates
		G4ThreeVector mag_center;
		// rotation matrix
		G4RotationMatrix* mag_rot = new G4RotationMatrix();
		// magnet solid object
		G4VSolid* mag_solid;

    //if( simPar->GetMagType(iMag) == "SOL" ) continue;

    if( simPar->GetMagName(iMag) != "Q1ER" && simPar->GetMagName(iMag) != "Q2ER" &&
      simPar->GetMagName(iMag) != "D2ER" && simPar->GetMagName(iMag) != "Q3ER" ) continue;

    //if( simPar->GetMagName(iMag) != "D2ER" ) continue;
    //if( simPar->GetMagName(iMag) != "Q1ER" ) continue;
    //if( simPar->GetMagName(iMag) != "Q2ER" ) continue;

    G4cout << "Magnet: " << simPar->GetMagType(iMag) << " " << simPar->GetMagName(iMag) << G4endl;
    //continue;

		if(simPar->GetMagType(iMag) == "D-pole" || simPar->GetMagType(iMag) == "Q-pole")
		{
			mag_vol = new G4Box(simPar->GetMagName(iMag)+"_vol",
				1*m/2,1*m/2,simPar->GetMagLength(iMag)/2); 

			// vector of start coordinates in the global ref. system shifted w.r.t. the ip
			G4ThreeVector magVector1(	simPar->GetMagStartPosX(iMag) - simPar->GetIpPosX(),
							simPar->GetMagStartPosY(iMag) - simPar->GetIpPosY(),
							simPar->GetMagStartPosZ(iMag) - simPar->GetIpPosZ());
			// vector of end coordinates in the global ref. system shifted w.r.t. the ip
			G4ThreeVector magVector2(	simPar->GetMagEndPosX(iMag) - simPar->GetIpPosX(),
							simPar->GetMagEndPosY(iMag) - simPar->GetIpPosY(),
							simPar->GetMagEndPosZ(iMag) - simPar->GetIpPosZ());
			// rotate coordinate vectors from global to local/ip system
			magVector1.rotateY(-simPar->GetIpTheta());
			magVector2.rotateY(-simPar->GetIpTheta());
			// check if the magnet length is the same as in the XML
			G4double L = sqrt(pow(magVector2.x()-magVector1.x(),2)+pow(magVector2.z()-magVector1.z(),2));
			if((L - simPar->GetMagLength(iMag))/simPar->GetMagLength(iMag) > 1e-4)
			{
				G4cout<<"\n\n=================================================================="<<G4endl;
				G4cout<<"[WARNING] DetectorConstruction::BuildMagnets"<<G4endl;
				G4cout<<"Not accurate magnet length calculation"<<G4endl;
				G4cout<<"L (calc) = "<<L/m<<" [m]"<<G4endl;
				G4cout<<"L (xml)  = "<<simPar->GetMagLength(iMag)/m<<" [m]"<<G4endl;
				G4cout<<"==================================================================\n"<<G4endl;
			}
				

      G4cout << "length: " << simPar->GetMagLength(iMag) << " " << L << G4endl;

			// coordinates of the magnet center
			mag_center = G4ThreeVector(	magVector1.x() + (magVector2.x()-magVector1.x())/2,
							magVector1.y() + (magVector2.y()-magVector1.y())/2,
							magVector1.z() + (magVector2.z()-magVector1.z())/2);

			// orientation angle of the magnet
			G4double rotY_angle = atan((magVector2.x()-magVector1.x())/(magVector2.z()-magVector1.z()))*rad;
			mag_rot->rotateY(-rotY_angle);

      G4cout << "center, angle: " << mag_center.x() << " " << mag_center.y() << " " << mag_center.z() << " " << rotY_angle << G4endl;

			if(simPar->GetMagName(iMag) == "D2ER")
			{
				mag_solid = new G4IntersectionSolid(simPar->GetMagName(iMag)+"_solid",
					bwd_beampipeVac_solid, mag_vol,mag_rot,mag_center);
			}
			else
			{
				mag_solid = new G4IntersectionSolid(simPar->GetMagName(iMag)+"_solid",
					beampipeUni_solid, mag_vol,mag_rot,mag_center);
			}
		}
		else if(simPar->GetMagType(iMag) == "SOL")
		{
			if(!simPar->IsBuildSolenoid()){continue;}

			mag_vol = new G4Tubs(simPar->GetMagName(iMag)+"_vol",
				0.0,simPar->GetSolRmax(),(simPar->GetSolZmax()-simPar->GetSolZmin())/2.0,0.,2*M_PI*rad); 
			// coordinates of the magnet center
			mag_center = simPar->GetSolOrigVec(); 

			if(!ip_box)
			{
				throw runtime_error(
				"[ERROR] DetectorConstruction::BuildMagnets ==> IP beam pipe or magnet are not created\n");
			}

			G4VSolid* vac_wo_ip_solid = new G4SubtractionSolid("vac_wo_ip_solid", 
				beampipeUni_solid, ip_box,0,
			G4ThreeVector(	0,
					0,
					simPar->GetIpBeampipeZmin()+
						(simPar->GetIpBeampipeZmax()-simPar->GetIpBeampipeZmin())/2));
			mag_solid = new G4IntersectionSolid(simPar->GetMagName(iMag)+"_solid",
				vac_wo_ip_solid, mag_vol,mag_rot,mag_center);
		}
		else
		{
			throw runtime_error(
			"[ERROR] DetectorConstruction::BuildMagnets ==> Wrong magnet type\n");
		}

		// magnet logic volume
		G4LogicalVolume* mag_log = new G4LogicalVolume(mag_solid,mat_vac, simPar->GetMagName(iMag)+"_log");

    G4cout << "Logical volume name: " << mag_log->GetName() << endl;

		// set tracking time limit in the volume
		mag_log->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,simPar->GetTrackTcut()));
		// magnet placement
		if(simPar->GetMagName(iMag) == "Q3ER")
		{
			new G4PVPlacement(mag_rot,mag_center,
				mag_log,"mag_"+simPar->GetMagName(iMag),world_log,false,0,true);
		}
		else
		{
			new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
				mag_log,"mag_"+simPar->GetMagName(iMag),world_log,false,0,true);
		}

		//- Create magnetic field
		G4MagneticField* mag_field = GetMagField(iMag,mag_log); 

		// apply the field to logic volume
		G4FieldManager* localFieldMgr = new G4FieldManager(mag_field);
		localFieldMgr->SetDetectorField(mag_field);
		localFieldMgr->CreateChordFinder(mag_field);
		mag_log->SetFieldManager(localFieldMgr, true);

		// apply the solenoid field to the IP beam pipe
		if(simPar->GetMagType(iMag) == "SOL" && simPar->IsBuildSolenoid())
		{
			if(!ip_log)
			{
				throw runtime_error(
				"[ERROR] DetectorConstruction::BuildMagnets ==> IP beam pipe is not built\n");
			}
			ip_log->SetFieldManager(localFieldMgr, true); 
		}

		//- Vacuum (subtract the magnet volume from the beam pipe volume)
		if(mag_center.z() > 15 * m) // Far-BWD side 
		{
			bwd_beampipeVac_solid = new G4SubtractionSolid("bwd_beampipeVac_solid",
				bwd_beampipeVac_solid,mag_vol,mag_rot,mag_center);

			if(simPar->GetMagName(iMag) == "D2ER")
			{
				G4VSolid* yoke_solid = new G4Box("yoke_solid",
					80 * cm / 2., 40 * cm / 2., simPar->GetMagLength(iMag)/2);
				G4Box* yoke_cut = new G4Box("yoke_cut",
					70 * cm / 2., 10 * cm / 2. + 3 * mm, 100 * m);
				yoke_solid = new G4SubtractionSolid(simPar->GetMagName(iMag)+"_yoke",
					yoke_solid, yoke_cut, 0, G4ThreeVector(0,0,0)); 

				yoke_cut = new G4Box("yoke_cut",
					20 * cm / 2., 30 * cm / 2., 100 * m);
				yoke_solid = new G4SubtractionSolid(simPar->GetMagName(iMag)+"_yoke",
					yoke_solid, yoke_cut, 0, 
					G4ThreeVector(30 * cm / 2. + 20 * cm / 2.,0,0)); 
				yoke_solid = new G4SubtractionSolid(simPar->GetMagName(iMag)+"_yoke",
					yoke_solid, yoke_cut, 0, 
					G4ThreeVector(-30 * cm / 2. - 20 * cm / 2.,0,0)); 

				

				G4LogicalVolume* yoke_log = new G4LogicalVolume(yoke_solid,mat_vac,
					simPar->GetMagName(iMag)+"_yoke_log");
				new G4PVPlacement(mag_rot,mag_center,yoke_log,"yoke_"+simPar->GetMagName(iMag),
					world_log,false,0,true);
				G4VisAttributes* yoke_vis = new G4VisAttributes(G4Color::Magenta());
				yoke_vis->SetVisibility(visibility);
				yoke_log->SetVisAttributes(yoke_vis);
			}
		}
		else
		{
			beampipeVac_solid = new G4SubtractionSolid("beampipeVac_solid",
				beampipeVac_solid,mag_vol,mag_rot,mag_center);
		}
	}

}//BuildMagnets

//_____________________________________________________________________________
void BeamAndrii20240604::BuildParticles() {

        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particleDef;
        particleDef = particleTable->FindParticle(simPar->GetBeamName().c_str());
        mass = particleDef->GetPDGMass();
        gamma = simPar->GetBeamGamma();

}//BuildParticles

//_____________________________________________________________________________
G4MagneticField* BeamAndrii20240604::GetMagField(G4int iMag, G4LogicalVolume* mag_log) {

	if(mass < 0 || gamma < 0)
	{
		throw runtime_error(
		"[ERROR] DetectorConstruction::GetMagField ==> Particles are not built\n");
	}

	G4MagneticField* mag_field;
	if(simPar->GetMagType(iMag) == "D-pole")
	{
		// By[T]*rho[m] = 3.33564*pc[GeV]
		// rho[m] = l[m]/(2*sin(theta[rad]/2)) (RBEND)
		// pc = mc^2*sqrt(gamma^2 - 1)
		G4double By = ((3.33564*(mass/GeV)*sqrt(pow(gamma,2)-1.0))/
			((simPar->GetMagLength(iMag)/m)/(2.*sin(simPar->GetMagAngle(iMag)/2.))))*tesla;
		mag_field = new G4UniformMagField(G4ThreeVector(0,By,0));

    G4cout << "dipole field: " << mass << " " << simPar->GetBeamName() << " " << gamma << " ";
    G4cout << simPar->GetMagLength(iMag) << " " << simPar->GetMagAngle(iMag) << " " << By << G4endl;

		//- Visualization
		G4VisAttributes* mag_vis = new G4VisAttributes(G4Color::Blue());
		mag_vis->SetVisibility(visibility);
		mag_log->SetVisAttributes(mag_vis);
	}
	else if(simPar->GetMagType(iMag) == "Q-pole")
	{
		// k1[1/m^2] = 0.2998*grad[T/m]/pc[GeV], k1[1/m^2] = k1l[1/m]/l[m]
		G4double grad = (3.33564*(simPar->GetMagK1L(iMag)/(1/m))/(simPar->GetMagLength(iMag)/m)*
			(mass/GeV)*sqrt(pow(gamma,2)-1.0))*(tesla/m);
		mag_field = new G4QuadrupoleMagField(grad);

    G4cout << "quadrupole field: " << mass << " " << simPar->GetBeamName() << " " << gamma << " ";
    G4cout << simPar->GetMagLength(iMag) << " " << simPar->GetMagK1L(iMag) << " " << grad << G4endl;
		//- Visualization
		G4VisAttributes* mag_vis = new G4VisAttributes(G4Color::Yellow());
		mag_vis->SetVisibility(visibility);
		mag_log->SetVisAttributes(mag_vis);
	}
	else if(simPar->GetMagType(iMag) == "SOL" && simPar->IsBuildSolenoid())
	{
/*
		// rotation matrix
		G4RotationMatrix* sol_rot = new G4RotationMatrix();
		// flip XZ coordinates from EIC to Geant4
		if(simPar->IsFlipXZSolenoid()){sol_rot->rotateY(M_PI*rad);}

		mag_field = new SolenoidMagField(	simPar->GetSolFileName(),
							simPar->GetSolRmax(),
							simPar->GetSolZmin(),
							simPar->GetSolZmax(),
							simPar->IsPrintSolenoid(),
							simPar->GetSolOrigVec(),
							sol_rot); 
		//- Visualization
		G4VisAttributes* mag_vis = new G4VisAttributes(G4Colour::White());
		mag_vis->SetVisibility(visibility);
		mag_log->SetVisAttributes(mag_vis);
*/

	}
	else
	{
		throw runtime_error(
		"[ERROR] DetectorConstruction::GetMagField ==> Wrong magnet type\n");
	}

	return mag_field;

}//GetMagField



/*
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

*/



















