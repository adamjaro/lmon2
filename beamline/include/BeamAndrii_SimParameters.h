#ifndef BeamAndrii_SimParameters_h
#define BeamAndrii_SimParameters_h

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

#include <vector>

using namespace std;

class BeamAndrii_SimParameters
{
public: 
	BeamAndrii_SimParameters();
	~BeamAndrii_SimParameters();
private:
    	G4double _world_dim_x;
    	G4double _world_dim_y;
    	G4double _world_dim_z;

	G4String _beampipeFileName;

	G4double _ext_length;
	G4double _ext_radius;
	G4double _ext_endpos;	

	G4double _ipPosX;
	G4double _ipPosY;
	G4double _ipPosZ;
	G4double _ipTheta;

	G4String _solFileName;
	G4bool _solBuild;
	G4bool _solPrint;
	G4bool _solFlipXZ;
	G4double _solRmax;
	G4double _solZmin;
	G4double _solZmax;
	G4double _solOrigPosX;
	G4double _solOrigPosY;
	G4double _solOrigPosZ;

	vector<G4String> _magName;
	vector<G4String> _magType;
	vector<G4double> _magLength;
	vector<G4double> _magStartPosX;
	vector<G4double> _magStartPosY;
	vector<G4double> _magStartPosZ;
	vector<G4double> _magEndPosX;
	vector<G4double> _magEndPosY;
	vector<G4double> _magEndPosZ;
	vector<G4double> _magAngle;
	vector<G4double> _magK1L;

	G4String _beamMonName;

	vector<G4double> _absSizeX;
	vector<G4double> _absSizeY;
	vector<G4double> _absSizeZ;
	vector<G4double> _absPosX;
	vector<G4double> _absPosY;
	vector<G4double> _absPosZ;

    	G4double _beamStartPosX;
    	G4double _beamStartPosY;
    	G4double _beamStartPosZ;
    	G4double _beamEndPosX;
    	G4double _beamEndPosY;
    	G4double _beamEndPosZ;
	G4double _beamEmitX;
	G4double _beamEmitY;
	G4double _beamAlphaX;
	G4double _beamAlphaY;
	G4double _beamBetaX;
	G4double _beamBetaY;
	G4double _beamEtaX;
	G4double _beamEtaY;
	G4double _beamEtaPrimeX;
	G4double _beamEtaPrimeY;

	G4double _beamMomSpread;
    	G4double _beamGamma;
    	G4String _beamName;
	G4String _beamType;
	G4double _beamTailKx;
	G4double _beamTailKy;
	G4double _beamTailMin;
	G4double _beamTailMax;

	G4String _generator;
	G4String _hepmcFileName;

	G4bool _storeTrj;
	G4double _storeTrjZmin;
	G4double _storeTrjZmax;
	G4double _storeTrjEmin;

	G4int _SRtype;
	G4int _reflType;
	G4bool _reflPrint;
	G4double _sigmaRatio;
        G4double _surfaceRoughness;
	G4double _autoCorrLength;

	G4double _gammaEcut;
	G4double _trackTcut;
	G4bool _straightTracksKillStatus;
	G4double _gammaDmin;
	G4double _gammaZend;

	G4double _ipBeampipeZmin;
	G4double _ipBeampipeZmax;

	G4double _meanFreePathFactor;
	G4double _bwdTrkCutZ;
public:
	// print functions
	void InitDefault();	
	void PrintGeoParameters();
	void PrintBeampipeFileName();
	void PrintBeamParameters();
	void PrintMagParameters();
	void PrintSolParameters();
	void CleanMagParameters();
	void PrintAbsParameters();
	void PrintIpBeampipeParameters();
	void CleanAbsParameters();
	void PrintSRprocParameters();
	void PrintGammaCut();
	void PrintMeanFreePathFactor();
	void ReadXML(G4String fName);
	void PrintStoreTrajectories();
	void PrintTrkGeoCut();
	void PrintGenerator();

	// get functions
	G4double GetWorldSizeX(){return _world_dim_x;}
	G4double GetWorldSizeY(){return _world_dim_y;}
	G4double GetWorldSizeZ(){return _world_dim_z;}

	G4String GetBeampipeFileName(){return _beampipeFileName;}

	G4double GetExtL(){return _ext_length;}
	G4double GetExtR(){return _ext_radius;}
	G4double GetExtE(){return _ext_endpos;}

	G4double GetIpBeampipeZmin(){return _ipBeampipeZmin;}
	G4double GetIpBeampipeZmax(){return _ipBeampipeZmax;}

	G4double GetIpPosX(){return _ipPosX;}
	G4double GetIpPosY(){return _ipPosY;}
	G4double GetIpPosZ(){return _ipPosZ;}
	G4double GetIpTheta(){return _ipTheta;}

	G4bool IsBuildSolenoid(){return _solBuild;}
	G4bool IsPrintSolenoid(){return _solPrint;}
	G4bool IsFlipXZSolenoid(){return _solFlipXZ;}
	G4String GetSolFileName(){return _solFileName;}
	G4double GetSolRmax(){return _solRmax;}
	G4double GetSolZmin(){return _solZmin;}
	G4double GetSolZmax(){return _solZmax;}
	G4ThreeVector GetSolOrigVec(){return G4ThreeVector(_solOrigPosX,_solOrigPosY,_solOrigPosZ);}

	G4int GetMagNum(){return _magName.size();}
	G4String GetMagName(G4int i){return _magName.at(i);}
	G4String GetMagType(G4int i){return _magType.at(i);}
	G4double GetMagLength(G4int i){return _magLength.at(i);}
	G4double GetMagStartPosX(G4int i){return _magStartPosX.at(i);}
	G4double GetMagStartPosY(G4int i){return _magStartPosY.at(i);}
	G4double GetMagStartPosZ(G4int i){return _magStartPosZ.at(i);}
	G4double GetMagEndPosX(G4int i){return _magEndPosX.at(i);}
	G4double GetMagEndPosY(G4int i){return _magEndPosY.at(i);}
	G4double GetMagEndPosZ(G4int i){return _magEndPosZ.at(i);}
	G4double GetMagAngle(G4int i){return _magAngle.at(i);}
	G4double GetMagK1L(G4int i){return _magK1L.at(i);}

	G4String GetBeamMonName(){return _beamMonName;}
	
	G4int GetAbsNum(){return _absSizeX.size();}
	G4double GetAbsSizeX(G4int i){return _absSizeX.at(i);}
	G4double GetAbsSizeY(G4int i){return _absSizeY.at(i);}
	G4double GetAbsSizeZ(G4int i){return _absSizeZ.at(i);}
	G4double GetAbsPosX(G4int i){return _absPosX.at(i);}
	G4double GetAbsPosY(G4int i){return _absPosY.at(i);}
	G4double GetAbsPosZ(G4int i){return _absPosZ.at(i);}

    	G4String GetBeamName(){return _beamName;}
    	G4String GetBeamType(){return _beamType;}
	G4double GetBeamTailKx(){return _beamTailKx;}
	G4double GetBeamTailKy(){return _beamTailKy;}
	G4double GetBeamTailMin(){return _beamTailMin;}
	G4double GetBeamTailMax(){return _beamTailMax;}
	G4double GetBeamStartPosX(){return _beamStartPosX;}
	G4double GetBeamStartPosY(){return _beamStartPosY;}
	G4double GetBeamStartPosZ(){return _beamStartPosZ;}
	G4double GetBeamEndPosX(){return _beamEndPosX;}
	G4double GetBeamEndPosY(){return _beamEndPosY;}
	G4double GetBeamEndPosZ(){return _beamEndPosZ;}
    	G4double GetBeamGamma(){return _beamGamma;}
	G4double GetBeamMomSpread(){return _beamMomSpread;}
	G4double GetBeamAlphaX(){return _beamAlphaX;}
	G4double GetBeamAlphaY(){return _beamAlphaY;}
	G4double GetBeamBetaX(){return _beamBetaX;}
	G4double GetBeamBetaY(){return _beamBetaY;}
	G4double GetBeamEtaX(){return _beamEtaX;}
	G4double GetBeamEtaY(){return _beamEtaY;}
	G4double GetBeamEtaPrimeX(){return _beamEtaPrimeX;}
	G4double GetBeamEtaPrimeY(){return _beamEtaPrimeY;}
	G4double GetBeamEmitX(){return _beamEmitX;}
	G4double GetBeamEmitY(){return _beamEmitY;}

	G4String GetGenerator(){return _generator;}
	G4String GetHepMcFileName(){return _hepmcFileName;}

	G4int GetSRtype(){return _SRtype;}
	G4int GetReflType(){return _reflType;}
	G4bool GetReflPrint(){return _reflPrint;}
	G4double GetSigmaRatio(){return _sigmaRatio;}
	G4double GetSurfRoughness(){return _surfaceRoughness;}
	G4double GetAutoCorrLength(){return _autoCorrLength;}

	G4double GetGammaEcut(){return _gammaEcut;}
	G4double GetTrackTcut(){return _trackTcut;}
	G4bool GetStraightTracksKillStatus(){return _straightTracksKillStatus;}
	G4double GetStraightTracksDmin(){return _gammaDmin;}
	G4double GetStraightTracksZend(){return _gammaZend;}

	G4bool GetStoreTrajectories(){return _storeTrj;}
	G4double GetStoreTrajectoriesZmin(){return _storeTrjZmin;}
	G4double GetStoreTrajectoriesZmax(){return _storeTrjZmax;}
	G4double GetStoreTrajectoriesEmin(){return _storeTrjEmin;}

	G4double GetMeanFreePathFactor(){return _meanFreePathFactor;}
	G4double GetTrackBwdGeoCut(){return _bwdTrkCutZ ;}

private:
	// set functions
	void SetWorldSizeX(G4double val){_world_dim_x = val; return;}
	void SetWorldSizeY(G4double val){_world_dim_y = val; return;}
	void SetWorldSizeZ(G4double val){_world_dim_z = val; return;}

	void SetBeampipeFileName(G4String val){_beampipeFileName = val; return;}

	void SetExtL(G4double val){_ext_length = val; return;}
	void SetExtR(G4double val){_ext_radius = val; return;}
	void SetExtE(G4double val){_ext_endpos = val; return;}

	void SetIpBeampipeZmin(G4double val){_ipBeampipeZmin = val; return;}
	void SetIpBeampipeZmax(G4double val){_ipBeampipeZmax = val; return;}

	void SetIpPosX(G4double val){_ipPosX = val; return;}
	void SetIpPosY(G4double val){_ipPosY = val; return;}
	void SetIpPosZ(G4double val){_ipPosZ = val; return;}
	void SetIpTheta(G4double val){_ipTheta = val; return;}

	void SetSolFileName(G4String val){_solFileName = val; return;}
	void SetSolRmax(G4double val){_solRmax = val; return;}
	void SetSolZmin(G4double val){_solZmin = val; return;}
	void SetSolZmax(G4double val){_solZmax = val; return;}
	void SetSolOrigPosX(G4double val){_solOrigPosX = val; return;}
	void SetSolOrigPosY(G4double val){_solOrigPosY = val; return;}
	void SetSolOrigPosZ(G4double val){_solOrigPosZ = val; return;}

	void SetMagParameters(G4String magName, G4String magType, G4double magLength, 
		G4double magStartPosX, G4double magStartPosY, G4double magStartPosZ,
		G4double magEndPosX, G4double magEndPosY, G4double magEndPosZ,
		G4double magAngle, G4double magK1L);

	void SetBeamMonName(G4String name){_beamMonName = name; return;}

	void SetAbsParameters(G4double absSizeX, G4double absSizeY, G4double absSizeZ,
	G4double absPosX, G4double absPosY, G4double absPosZ);

    	void SetBeamName(G4String val){_beamName = val; return;}
    	void SetBeamType(G4String val){_beamType = val; return;}
	void SetBeamTailKx(G4double val){_beamTailKx = val; return;}
	void SetBeamTailKy(G4double val){_beamTailKy = val; return;}
    	void SetBeamTailMin(G4double val){_beamTailMin = val; return;}
    	void SetBeamTailMax(G4double val){_beamTailMax = val; return;}
	void SetBeamStartPosX(G4double val){_beamStartPosX = val; return;}
	void SetBeamStartPosY(G4double val){_beamStartPosY = val; return;}
	void SetBeamStartPosZ(G4double val){_beamStartPosZ = val; return;}
	void SetBeamEndPosX(G4double val){_beamEndPosX = val; return;}
	void SetBeamEndPosY(G4double val){_beamEndPosY = val; return;}
	void SetBeamEndPosZ(G4double val){_beamEndPosZ = val; return;}
    	void SetBeamGamma(G4double val){_beamGamma = val; return;}
	void SetBeamMomSpread(G4double val){_beamMomSpread = val; return;}
	void SetBeamAlphaX(G4double val){_beamAlphaX = val; return;}
	void SetBeamAlphaY(G4double val){_beamAlphaY = val; return;}
	void SetBeamBetaX(G4double val){_beamBetaX = val; return;}
	void SetBeamBetaY(G4double val){_beamBetaY = val; return;}
	void SetBeamEtaX(G4double val){_beamEtaX = val; return;}
	void SetBeamEtaY(G4double val){_beamEtaY = val; return;}
	void SetBeamEtaPrimeX(G4double val){_beamEtaPrimeX = val; return;}
	void SetBeamEtaPrimeY(G4double val){_beamEtaPrimeY = val; return;}
	void SetBeamEmitX(G4double val){_beamEmitX = val; return;}
	void SetBeamEmitY(G4double val){_beamEmitY = val; return;}

	void SetGenerator(G4String val){_generator = val; return;}
	void SetHepMcFileName(G4String val){_hepmcFileName = val; return;}

	void SetSrType(G4int val){_SRtype = val; return;}
	void SetReflType(G4int val){_reflType = val; return;}
	void SetReflPrint(G4int val){if(val == 1){_reflPrint=true;}else{_reflPrint=false;}return;}
	void SetSigmaRatio(G4double val){_sigmaRatio = val; return;}
	void SetSurfRoughness(G4double val){_surfaceRoughness = val; return;}
	void SetAutoCorrLength(G4double val){_autoCorrLength = val; return;}

	void SetGammaEcut(G4double val){_gammaEcut = val; return;}
	void SetTrackTcut(G4double val){_trackTcut = val; return;}
	void SetStraightTracksKillStatus(G4int val)
	{
		if(val == 1){_straightTracksKillStatus=true;}
		else{_straightTracksKillStatus=false;} 
		return;
	}
	void SetStraightTracksDmin(G4double val){_gammaDmin = val; return;}
	void SetStraightTracksZend(G4double val){_gammaZend = val; return;}

	G4bool FileExists(const G4String& fName);

	void SetStoreTrajectories(G4int val){if(val == 1){_storeTrj=true;}else{_storeTrj=false;}return;}
	void SetStoreTrajectoriesZmin(G4double val){_storeTrjZmin = val;return;}
	void SetStoreTrajectoriesZmax(G4double val){_storeTrjZmax = val;return;}
	void SetStoreTrajectoriesEmin(G4double val){_storeTrjEmin = val;return;}

	void SetMeanFreePathFactor(G4double val){_meanFreePathFactor = val; return;}
	void SetTrackBwdGeoCut(G4double val){_bwdTrkCutZ = val; return;}
};

#endif
