#include "BeamAndrii_SimParameters.h"  
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <sys/stat.h>

using namespace boost::property_tree;

BeamAndrii_SimParameters::BeamAndrii_SimParameters()
{}    	

BeamAndrii_SimParameters::~BeamAndrii_SimParameters()
{}

void BeamAndrii_SimParameters::ReadXML(G4String fName)
{
    	// Create empty property tree object
    	ptree tree;

    	// Parse the XML into the property tree.
	if(!FileExists(fName)){return;}
    	read_xml(fName, tree);

    	// Use get_child to find the node containing the a specific data, and iterate over
    	// its children. If the path cannot be resolved, get_child throws.
    	// A C++11 for-range loop would also work.

	// world
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.world"))
	{
        	// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "x" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetWorldSizeX(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "y" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetWorldSizeY(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "z" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetWorldSizeZ(stod(child_node.second.data())*m);
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n"); 
    	}

	// beam pipe
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.beampipe"))
	{
		// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "fname")
			SetBeampipeFileName(child_node.second.data());
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// ip beam pipe
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.ipbeampipe"))
	{
		// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "zmin" &&
			child_node.second.get<string>("<xmlattr>.unit") == "cm")
			SetIpBeampipeZmin(stod(child_node.second.data())*cm);
		else if(child_node.second.get<string>("<xmlattr>.name") == "zmax" &&
			child_node.second.get<string>("<xmlattr>.unit") == "cm")
			SetIpBeampipeZmax(stod(child_node.second.data())*cm);
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// extension
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.extension"))
	{
        	// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "l" &&
			child_node.second.get<string>("<xmlattr>.unit") == "cm")
			SetExtL(stod(child_node.second.data())*cm);
		else if(child_node.second.get<string>("<xmlattr>.name") == "r" &&
			child_node.second.get<string>("<xmlattr>.unit") == "cm")
			SetExtR(stod(child_node.second.data())*cm);
		else if(child_node.second.get<string>("<xmlattr>.name") == "e" &&
			child_node.second.get<string>("<xmlattr>.unit") == "cm")
			SetExtE(stod(child_node.second.data())*cm);
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n"); 
    	}

	// SR model
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.srmodel"))
	{
		// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "model")
			SetSrType(stoi(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "type")
			SetReflType(stoi(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "print")
			SetReflPrint(stoi(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "sigma" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetSigmaRatio(stod(child_node.second.data()));	
		else if(child_node.second.get<string>("<xmlattr>.name") == "roughness" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetSurfRoughness(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "corrlength" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetAutoCorrLength(stod(child_node.second.data())*m);
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// Store SR trajectories
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.storetrj"))
	{
		// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "store")
			SetStoreTrajectories(stoi(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "Zmin" &&
			child_node.second.get<string>("<xmlattr>.unit") == "cm")
			SetStoreTrajectoriesZmin(stod(child_node.second.data())*cm);
		else if(child_node.second.get<string>("<xmlattr>.name") == "Zmax" &&
			child_node.second.get<string>("<xmlattr>.unit") == "cm")
			SetStoreTrajectoriesZmax(stod(child_node.second.data())*cm);
		else if(child_node.second.get<string>("<xmlattr>.name") == "Emin" &&
			child_node.second.get<string>("<xmlattr>.unit") == "GeV")
			SetStoreTrajectoriesEmin(stod(child_node.second.data())*GeV);
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// gamma cuts
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.gammacut"))
	{
		// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "ecut" &&
			child_node.second.get<string>("<xmlattr>.unit") == "eV")
			SetGammaEcut(stod(child_node.second.data())*eV);
		else if(child_node.second.get<string>("<xmlattr>.name") == "trackkill")
			SetStraightTracksKillStatus(stoi(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "Dmin" &&
			child_node.second.get<string>("<xmlattr>.unit") == "mm")
			SetStraightTracksDmin(stod(child_node.second.data())*mm);
		else if(child_node.second.get<string>("<xmlattr>.name") == "Zend" &&
			child_node.second.get<string>("<xmlattr>.unit") == "cm")
			SetStraightTracksZend(stod(child_node.second.data())*cm);
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// mean free path factor 
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.meanfreepath"))
	{
		// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "MeanFreePathFactor" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetMeanFreePathFactor(stod(child_node.second.data()));
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// tracking time cut
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.timecut"))
	{
		// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "tcut" &&
			child_node.second.get<string>("<xmlattr>.unit") == "s")
			SetTrackTcut(stod(child_node.second.data())*s);
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// tracking geometry cut on the BWD side
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.bwdgeocut"))
	{
		// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "bwdTrkCutZ" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetTrackBwdGeoCut(stod(child_node.second.data())*m);
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// beam
	BOOST_FOREACH(ptree::value_type &child_node, tree.get_child("xml.generator"))
	{
        	// The data function is used to access the data stored in a node.
		if(child_node.second.get<string>("<xmlattr>.name") == "genname")
			SetGenerator(child_node.second.data());
		else if(child_node.second.get<string>("<xmlattr>.name") == "part")
			SetBeamName(child_node.second.data());
		else if(child_node.second.get<string>("<xmlattr>.name") == "posx1" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamStartPosX(stod(child_node.second.data())*m); 
		else if(child_node.second.get<string>("<xmlattr>.name") == "posy1" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamStartPosY(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "posz1" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamStartPosZ(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "posx2" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamEndPosX(stod(child_node.second.data())*m); 
		else if(child_node.second.get<string>("<xmlattr>.name") == "posy2" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamEndPosY(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "posz2" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamEndPosZ(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "gamma" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamGamma(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "spread" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamMomSpread(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "alphax" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamAlphaX(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "alphay" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamAlphaY(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "betax" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamBetaX(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "betay" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamBetaY(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "etax" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamEtaX(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "etay" &&
			child_node.second.get<string>("<xmlattr>.unit") == "m")
			SetBeamEtaY(stod(child_node.second.data())*m);
		else if(child_node.second.get<string>("<xmlattr>.name") == "etapx" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamEtaPrimeX(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "etapy" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamEtaPrimeY(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "emitx" &&
			child_node.second.get<string>("<xmlattr>.unit") == "nm")
			SetBeamEmitX(stod(child_node.second.data())*nm);
		else if(child_node.second.get<string>("<xmlattr>.name") == "emity" &&
			child_node.second.get<string>("<xmlattr>.unit") == "nm")
			SetBeamEmitY(stod(child_node.second.data())*nm);
		else if(child_node.second.get<string>("<xmlattr>.name") == "type")
			SetBeamType(child_node.second.data());
		else if(child_node.second.get<string>("<xmlattr>.name") == "tailkx" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamTailKx(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "tailky" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamTailKy(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "tailmin" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamTailMin(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "tailmax" &&
			child_node.second.get<string>("<xmlattr>.unit") == "")
			SetBeamTailMax(stod(child_node.second.data()));
		else if(child_node.second.get<string>("<xmlattr>.name") == "hepmcfilename")
			SetHepMcFileName(child_node.second.data());
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n"); 
    	}

	// absorbers
	// clean vectors
	CleanAbsParameters();
	G4double absSizeX, absSizeY, absSizeZ, absPosX, absPosY, absPosZ;
	
	BOOST_FOREACH(ptree::value_type &child_tree, tree.get_child("xml.absorbers"))
	{
		if (child_tree.first == "absorber")
		{
			BOOST_FOREACH(ptree::value_type &child_node, child_tree.second)
			{
				if(child_node.second.get<string>("<xmlattr>.name") == "x" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					absSizeX = stod(child_node.second.data())*cm; 
				else if(child_node.second.get<string>("<xmlattr>.name") == "y" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					absSizeY = stod(child_node.second.data())*cm; 
				else if(child_node.second.get<string>("<xmlattr>.name") == "z" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					absSizeZ = stod(child_node.second.data())*cm; 
				else if(child_node.second.get<string>("<xmlattr>.name") == "posx" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					absPosX = stod(child_node.second.data())*cm; 
				else if(child_node.second.get<string>("<xmlattr>.name") == "posy" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					absPosY = stod(child_node.second.data())*cm; 
				else if(child_node.second.get<string>("<xmlattr>.name") == "posz" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					absPosZ = stod(child_node.second.data())*cm; 
				else
					throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
			}
			// fill vectors with the new numbers
			SetAbsParameters(absSizeX,absSizeY,absSizeZ,absPosX,absPosY,absPosZ);
		}
		else
			throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
	}

	// magnets
	// clean vectors
	CleanMagParameters();
	G4String magName, magType;
	G4double magLength, magStartPosX, magStartPosY, magStartPosZ, magEndPosX, magEndPosY, magEndPosZ;
	G4double magAngle, magK1L;

	BOOST_FOREACH(ptree::value_type &child_tree, tree.get_child("xml.magnets"))
	{
		if (child_tree.first == "magnet")
		{
			BOOST_FOREACH(ptree::value_type &child_node, child_tree.second)
			{
				if(child_node.second.get<string>("<xmlattr>.name") == "name")
					magName = child_node.second.data();
				else if(child_node.second.get<string>("<xmlattr>.name") == "type")
					magType = child_node.second.data();
				else if(child_node.second.get<string>("<xmlattr>.name") == "length" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					magLength = stod(child_node.second.data())*m; 
				else if(child_node.second.get<string>("<xmlattr>.name") == "posx1" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					magStartPosX = stod(child_node.second.data())*m; 
				else if(child_node.second.get<string>("<xmlattr>.name") == "posy1" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					magStartPosY = stod(child_node.second.data())*m;
				else if(child_node.second.get<string>("<xmlattr>.name") == "posz1" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					magStartPosZ = stod(child_node.second.data())*m;
				else if(child_node.second.get<string>("<xmlattr>.name") == "posx2" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					magEndPosX = stod(child_node.second.data())*m; 
				else if(child_node.second.get<string>("<xmlattr>.name") == "posy2" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					magEndPosY = stod(child_node.second.data())*m;
				else if(child_node.second.get<string>("<xmlattr>.name") == "posz2" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					magEndPosZ = stod(child_node.second.data())*m;
				else if(child_node.second.get<string>("<xmlattr>.name") == "angle" &&
					child_node.second.get<string>("<xmlattr>.unit") == "rad")
					magAngle = stod(child_node.second.data())*rad;
				else if(child_node.second.get<string>("<xmlattr>.name") == "k1l" &&
					child_node.second.get<string>("<xmlattr>.unit") == "1/m")
					magK1L = stod(child_node.second.data())*1/m;
				else
					throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
			}
			// fill vectors with the new numbers
			SetMagParameters(magName,magType,magLength,magStartPosX,magStartPosY,magStartPosZ,
				magEndPosX,magEndPosY,magEndPosZ,magAngle,magK1L);
		}
		else if (child_tree.first == "ip")
		{
			BOOST_FOREACH(ptree::value_type &child_node, child_tree.second)
			{
				if(child_node.second.get<string>("<xmlattr>.name") == "posx" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					SetIpPosX(stod(child_node.second.data())*m); 
				else if(child_node.second.get<string>("<xmlattr>.name") == "posy" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					SetIpPosY(stod(child_node.second.data())*m); 
				else if(child_node.second.get<string>("<xmlattr>.name") == "posz" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					SetIpPosZ(stod(child_node.second.data())*m); 
				else if(child_node.second.get<string>("<xmlattr>.name") == "theta" &&
					child_node.second.get<string>("<xmlattr>.unit") == "rad")
					SetIpTheta(stod(child_node.second.data())*rad); 
				else
					throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
			}
		}
		else if (child_tree.first == "monitor")
		{
			BOOST_FOREACH(ptree::value_type &child_node, child_tree.second)
			{
				if(child_node.second.get<string>("<xmlattr>.name") == "magname")
					SetBeamMonName(child_node.second.data()); 
				else
					throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
			}
		}
		else if (child_tree.first == "solenoid")
		{
			BOOST_FOREACH(ptree::value_type &child_node, child_tree.second)
			{
				if(child_node.second.get<string>("<xmlattr>.name") == "name")
					magName = child_node.second.data();
				else if(child_node.second.get<string>("<xmlattr>.name") == "type")
					magType = child_node.second.data();
				else if(child_node.second.get<string>("<xmlattr>.name") == "build")
					_solBuild = stoi(child_node.second.data());
				else if(child_node.second.get<string>("<xmlattr>.name") == "fname")
					SetSolFileName(child_node.second.data());
				else if(child_node.second.get<string>("<xmlattr>.name") == "rmax" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					SetSolRmax(stod(child_node.second.data())*m);
				else if(child_node.second.get<string>("<xmlattr>.name") == "zmin" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					SetSolZmin(stod(child_node.second.data())*m);
				else if(child_node.second.get<string>("<xmlattr>.name") == "zmax" &&
					child_node.second.get<string>("<xmlattr>.unit") == "m")
					SetSolZmax(stod(child_node.second.data())*m);
				else if(child_node.second.get<string>("<xmlattr>.name") == "print")
					_solPrint = stoi(child_node.second.data());
				else if(child_node.second.get<string>("<xmlattr>.name") == "origposx" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					SetSolOrigPosX(stod(child_node.second.data())*cm);
				else if(child_node.second.get<string>("<xmlattr>.name") == "origposy" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					SetSolOrigPosY(stod(child_node.second.data())*cm);
				else if(child_node.second.get<string>("<xmlattr>.name") == "origposz" &&
					child_node.second.get<string>("<xmlattr>.unit") == "cm")
					SetSolOrigPosZ(stod(child_node.second.data())*cm);
				else if(child_node.second.get<string>("<xmlattr>.name") == "flipxz")
					_solFlipXZ = stoi(child_node.second.data());
				else
					throw runtime_error("[ERROR] SimParameters::ReadXML ==> Bad XML data\n");
			}
			// fill vectors with the new numbers
			SetMagParameters(magName,magType,0,0,0,0,0,0,0,0,0);
		}
	}
	
	return;
}

void BeamAndrii_SimParameters::SetMagParameters(G4String magName, G4String magType, G4double magLength, 
	G4double magStartPosX, G4double magStartPosY, G4double magStartPosZ,
	G4double magEndPosX, G4double magEndPosY, G4double magEndPosZ,
	G4double magAngle, G4double magK1L)
{
	_magName.push_back(magName);
	_magType.push_back(magType);
	_magLength.push_back(magLength);
	_magStartPosX.push_back(magStartPosX);
	_magStartPosY.push_back(magStartPosY);
	_magStartPosZ.push_back(magStartPosZ);
	_magEndPosX.push_back(magEndPosX);
	_magEndPosY.push_back(magEndPosY);
	_magEndPosZ.push_back(magEndPosZ);
	_magAngle.push_back(magAngle);
	_magK1L.push_back(magK1L);

	return;
}

void BeamAndrii_SimParameters::CleanMagParameters()
{
	_magName.clear();
	_magType.clear();
	_magLength.clear();
	_magStartPosX.clear();
	_magStartPosY.clear();
	_magStartPosZ.clear();
	_magEndPosX.clear();
	_magEndPosY.clear();
	_magEndPosZ.clear();
	_magAngle.clear();
	_magK1L.clear();

	return;
}

void BeamAndrii_SimParameters::SetAbsParameters(G4double absSizeX, G4double absSizeY, G4double absSizeZ,
	G4double absPosX, G4double absPosY, G4double absPosZ)
{
	_absSizeX.push_back(absSizeX);
	_absSizeY.push_back(absSizeY);
	_absSizeZ.push_back(absSizeZ);
	_absPosX.push_back(absPosX);
	_absPosY.push_back(absPosY);
	_absPosZ.push_back(absPosZ);

	return;
}

void BeamAndrii_SimParameters::CleanAbsParameters()  
{
	_absSizeX.clear();
	_absSizeY.clear();
	_absSizeZ.clear();
	_absPosX.clear();
	_absPosY.clear();
	_absPosZ.clear();

	return;
}

void BeamAndrii_SimParameters::InitDefault()
{
	_world_dim_x =  10.0*m;
    	_world_dim_y =  10.0*m;
    	_world_dim_z =  80.0*m;

	_beampipeFileName = "../geometry/beam_pipe_20231011_Zup.obj";

	_ext_length = 370*cm;
	_ext_radius = 20.0*cm;
	_ext_endpos = -3350.0*cm;	

	_ipBeampipeZmin = -0.67*m; 
	_ipBeampipeZmax = 0.686923*m;

	_ipPosX = -589.77164158437*m;
	_ipPosY = 0*m;
	_ipPosZ = 5.6843418860808e-14*m;
	_ipTheta = -3.1335926535898*rad;

	_solFileName = "../geometry/MARCO_v.6.4.1.1.3_1.7T_Magnetic_Field_Map_2022_11_14_rad_coords_cm_T.txt";
	_solBuild = true;
	_solPrint = true;
	_solFlipXZ = true;
	_solRmax =  5.0*m;
	_solZmin = -5.0*m;
	_solZmax =  5.0*m;
	_solOrigPosX = 0.0*cm;
	_solOrigPosY = 0.0*cm;
	_solOrigPosZ = 8.0*cm; // in the custom-build Geant4 coord. system (in the eic-shell = -8 cm)

	_magName = {"D2EF2","D2EF1","D1EF","Q1EF","Q0EF"};
	_magType = {"D-pole","D-pole","D-pole","Q-pole","Q-pole"};
	_magLength = {	2.7260154313337*m,
			2.7260154313337*m,
			0.89140008356876*m,
			1.61*m,
			1.2*m};
	_magStartPosX = {	-589.55061312277*m,
				-589.51933973323*m,
				-589.52285141176*m,
				-589.67268263994*m,
				-589.7156421817*m};
	_magStartPosY = {0*m,0*m,0*m,0*m,0*m};
	_magStartPosZ = {	37.15030382285*m,
				34.166071692045*m,
				31.181677766382*m,
				12.369604162111*m,
				6.9997760011947*m};
	_magEndPosX = {	-589.5206719963*m,
			-589.52117182359*m,
			-589.52931400515*m,
			-589.68556250255*m,
			-589.7252420793*m};
	_magEndPosY = {0*m,0*m,0*m,0*m,0*m};
	_magEndPosZ = {	34.42446825757*m,
			31.440072307701*m,
			30.290301193385*m,
			10.759655681836*m,
			5.7998144009899*m};
	_magAngle = {	0.011655839693698*rad,
			0.011655839693698*rad,
			0.0015*rad,
			0*rad,
			0*rad};
	_magK1L = {	0*1/m,
			0*1/m,
			0*1/m,
			0.1611384495*1/m,
			-0.261830778*1/m};

	_beamMonName = "Q0EF";
	
	_absSizeX = {30*cm,10*cm,10*cm,50*cm};
	_absSizeY = {30*cm,10*cm,10*cm,50*cm};
	_absSizeZ = {5*cm,5*cm,5*cm,5*cm};
	_absPosX = {0*cm,-12*cm,14*cm,0*cm};
	_absPosY = {0*cm,0*cm,0*cm,0*cm};
	_absPosZ = {1527*cm,475*cm,-550*cm,-3720*cm};
	
    	_beamStartPosX = -589.55313475589*m; // "SQ2EF_6 - X"
    	_beamStartPosY = 0*m; // "SQ2EF_6 - Y"
    	_beamStartPosZ = 37.300282625907*m; // "SQ2EF_6 - Z"
    	_beamEndPosX = -589.55061312277*m; // "OWW_SH - X"
    	_beamEndPosY = 0*m; // "OWW_SH - Y"
    	_beamEndPosZ = 37.15030382285*m; // "OWW_SH - Z"
	_beamEmitX = 24.0*nm;
	_beamEmitY =  2.0*nm;
	_beamAlphaX = -7.031713181; // "OWW_SH - ALX"
	_beamAlphaY = -1.226925635; // "OWW_SH - ALY"
	_beamBetaX = 161.7921567*m; // "OWW_SH - BETX"
	_beamBetaY = 121.6383193*m; // "OWW_SH - BETY"
	_beamEtaX =  -0.07618056746*m; // "OWW_SH - DX"
	_beamEtaY = 0.0*m; // "OWW_SH - DY"
	_beamEtaPrimeX = 0.02481186422; // "OWW_SH - DPX"
	_beamEtaPrimeY = 0.0*m; // "OWW_SH - DPY"
	_beamMomSpread = 0.109e-2;
	_beamGamma = 34924.26476;
	_beamName = "e-";
	_beamType = "core";
	_beamTailKx = 0;
	_beamTailKy = 0;
	_beamTailMin = 0;
	_beamTailMax = 0;

	_generator = "Beam";
	_hepmcFileName = "";

	_SRtype = 1;
	_reflType = 3;
	_reflPrint = true;
	_sigmaRatio = 5e-3;
	_surfaceRoughness = 50e-9*m;
	_autoCorrLength = 10000e-9*m;

	_gammaEcut = -1*keV;
	_trackTcut =  1*s;
	_straightTracksKillStatus = false;
	_gammaDmin = 0;
	_gammaZend = 0;

	_storeTrj = false;
	_storeTrjZmin = -67.0*cm;
	_storeTrjZmax =  68.7*cm;
	_storeTrjEmin = 10e-6*GeV;

	_meanFreePathFactor = 1.0;
	_bwdTrkCutZ = 1e9*m;

	return;
}

void BeamAndrii_SimParameters::PrintTrkGeoCut()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintTrkGeoCut ==> Tracking geometry cut"<<G4endl;
	G4cout<<" BWD cut Z = "<<GetTrackBwdGeoCut()/m<<" [m]"<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintGeoParameters()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintGeoParameters ==> Geometry parameters"<<G4endl;
	G4cout<<" World size X = "<<GetWorldSizeX()/m<<" [m]"<<G4endl;
	G4cout<<" World size Y = "<<GetWorldSizeY()/m<<" [m]"<<G4endl;
	G4cout<<" World size Z = "<<GetWorldSizeZ()/m<<" [m]"<<G4endl;
	G4cout<<" Extension length = "<<GetExtL()/cm<<" [cm]"<<G4endl;
	G4cout<<" Extension radius = "<<GetExtR()/cm<<" [cm]"<<G4endl;
	G4cout<<" Extension end position = "<<GetExtE()/cm<<" [cm]"<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintStoreTrajectories()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintStoreTrajectories ==> Store SR photon trajectories"<<G4endl;
	G4cout<<"Store trajectory flag = "<<GetStoreTrajectories()<<G4endl;
	G4cout<<"	0 - false"<<G4endl;
	G4cout<<"	1 - true"<<G4endl;
	G4cout<<"Zmin = "<<GetStoreTrajectoriesZmin()/cm<<" [cm]"<<G4endl;
	G4cout<<"Zmax = "<<GetStoreTrajectoriesZmax()/cm<<" [cm]"<<G4endl;
	G4cout<<"Emin = "<<GetStoreTrajectoriesEmin()/GeV<<" [GeV]"<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintGenerator()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintGenerator ==> Particle Geneator: "<<GetGenerator()<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	if(GetGenerator() == "Beam")
		PrintBeamParameters();
	else
		G4cout<<"[INFO] HepMC file name: "<<GetHepMcFileName()<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;		
	return;
}

void BeamAndrii_SimParameters::PrintBeamParameters()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintBeamParameters ==> Beam parameters"<<G4endl;
	G4cout<<"Beam particles = "<<GetBeamName()<<G4endl;
	G4cout<<"Start position X = "<<GetBeamStartPosX()/m<<" [m]"<<G4endl;
	G4cout<<"Start position Y = "<<GetBeamStartPosY()/m<<" [m]"<<G4endl;
	G4cout<<"Start position Z = "<<GetBeamStartPosZ()/m<<" [m]"<<G4endl;
	G4cout<<"End position X = "<<GetBeamEndPosX()/m<<" [m]"<<G4endl;
	G4cout<<"End position Y = "<<GetBeamEndPosY()/m<<" [m]"<<G4endl;
	G4cout<<"End position Z = "<<GetBeamEndPosZ()/m<<" [m]"<<G4endl;
	G4cout<<"Gamma = "<<GetBeamGamma()<<" []"<<G4endl;
	G4cout<<"Momentum spread = "<<GetBeamMomSpread()<<" []"<<G4endl;
	G4cout<<"Alpha X = "<<GetBeamAlphaX()/m<<" []"<<G4endl;
	G4cout<<"Alpha Y = "<<GetBeamAlphaY()/m<<" []"<<G4endl;
	G4cout<<"Beta X = "<<GetBeamBetaX()/m<<" [m]"<<G4endl;
	G4cout<<"Beta Y = "<<GetBeamBetaY()/m<<" [m]"<<G4endl;
	G4cout<<"Dispersion X = "<<GetBeamEtaX()/m<<" [m]"<<G4endl;
	G4cout<<"Dispersion Y = "<<GetBeamEtaY()/m<<" [m]"<<G4endl;
	G4cout<<"Dispersion PX = "<<GetBeamEtaPrimeX()/m<<" []"<<G4endl;
	G4cout<<"Dispersion PY = "<<GetBeamEtaPrimeY()/m<<" []"<<G4endl;
	G4cout<<"Emittance X = "<<GetBeamEmitX()/nm<<" [nm]"<<G4endl;
	G4cout<<"Emittance Y = "<<GetBeamEmitY()/nm<<" [nm]"<<G4endl;
	G4cout<<"Beam type = "<<GetBeamType()<<G4endl;
	G4cout<<"Tail Kx = "<<GetBeamTailKx()<<G4endl;
	G4cout<<"Tail Ky = "<<GetBeamTailKy()<<G4endl;
	G4cout<<"Tail Min = "<<GetBeamTailMin()<<" [nSigma]"<<G4endl;
	G4cout<<"Tail Max = "<<GetBeamTailMax()<<" [nSigma]"<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintSolParameters()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintSolParameters ==> Solenoid parameters"<<G4endl;
	G4cout<<"Build field = "<<IsBuildSolenoid()<<G4endl;
	G4cout<<"File name = "<<GetSolFileName()<<G4endl;
	G4cout<<"Max R = "<<GetSolRmax()/m<<" [m]"<<G4endl;
	G4cout<<"Min Z = "<<GetSolZmin()/m<<" [m]"<<G4endl;
	G4cout<<"Max Z = "<<GetSolZmax()/m<<" [m]"<<G4endl;
	G4cout<<"Print field = "<<IsPrintSolenoid()<<G4endl;
	G4ThreeVector solOrigVec = GetSolOrigVec();
	G4cout<<"Origin position = ("
		<<solOrigVec.x()/cm<<";"<<solOrigVec.y()/cm<<";"<<solOrigVec.z()/cm<<") [cm]"<<G4endl;
	G4cout<<"Flip XZ = "<<IsFlipXZSolenoid()<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintMagParameters()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintMagParameters ==> Magnet parameters"<<G4endl;
	G4cout<<"IP position X = "<<GetIpPosX()/m<<" [m]"<<G4endl;
	G4cout<<"IP position Y = "<<GetIpPosY()/m<<" [m]"<<G4endl;
	G4cout<<"IP position Z = "<<GetIpPosZ()/m<<" [m]"<<G4endl;
	G4cout<<"IP orientation angle = "<<GetIpTheta()/rad<<" [rad]"<<G4endl;
	G4cout<<G4endl;
	for(G4int iMag = 0; iMag < GetMagNum(); iMag++)
	{
		G4cout<<"ID = "<<iMag<<G4endl;
		G4cout<<"  Name = "<<GetMagName(iMag)<<G4endl;
		G4cout<<"  Type = "<<GetMagType(iMag)<<G4endl;

		if(GetMagType(iMag) == "SOL")
		{
			PrintSolParameters();
			continue;
		}

		G4cout<<"  Length = "<<GetMagLength(iMag)/m<<" [m]"<<G4endl;
		G4cout<<"  Start position X = "<<GetMagStartPosX(iMag)/m<<" [m]"<<G4endl;
		G4cout<<"  Start position Y = "<<GetMagStartPosY(iMag)/m<<" [m]"<<G4endl;
		G4cout<<"  Start position Z = "<<GetMagStartPosZ(iMag)/m<<" [m]"<<G4endl;
		G4cout<<"  End position X = "<<GetMagEndPosX(iMag)/m<<" [m]"<<G4endl;
		G4cout<<"  End position Y = "<<GetMagEndPosY(iMag)/m<<" [m]"<<G4endl;
		G4cout<<"  End position Z = "<<GetMagEndPosZ(iMag)/m<<" [m]"<<G4endl;
		G4cout<<"  Dipole angle = "<<GetMagAngle(iMag)/rad<<" [rad]"<<G4endl;
		G4cout<<"  Quadrupole strength K1L = "<<GetMagK1L(iMag)/(1/m)<<" [1/m)]"<<G4endl;
	}
	G4cout<<"Beam monitor name = "<<GetBeamMonName()<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintAbsParameters()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintAbsParameters ==> Absorber parameters"<<G4endl;
	for(G4int iAbs = 0; iAbs < GetAbsNum(); iAbs++)
	{
		G4cout<<"ID = "<<iAbs<<G4endl;
		G4cout<<"  Size X = "<<GetAbsSizeX(iAbs)/cm<<" [cm]"<<G4endl; 
		G4cout<<"  Size Y = "<<GetAbsSizeY(iAbs)/cm<<" [cm]"<<G4endl; 
		G4cout<<"  Size Z = "<<GetAbsSizeZ(iAbs)/cm<<" [cm]"<<G4endl; 
		G4cout<<"  Position X = "<<GetAbsPosX(iAbs)/cm<<" [cm]"<<G4endl; 
		G4cout<<"  Position Y = "<<GetAbsPosY(iAbs)/cm<<" [cm]"<<G4endl; 
		G4cout<<"  Position Z = "<<GetAbsPosZ(iAbs)/cm<<" [cm]"<<G4endl; 
	}
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintSRprocParameters()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintSRprocParameters"<<G4endl;
	G4cout<<" SRtype = "<<GetSRtype()<<G4endl;
	G4cout<<" 	1 - G4SynchrotronRadiation"<<G4endl;
	G4cout<<" 	2 - G4SynchrotronRadiationInMat"<<G4endl;
	G4cout<<" SR photon reflection type = "<<GetReflType()<<G4endl;
	G4cout<<" 	0 - Geant4 (Xray,Névot-Croce)"<<G4endl;
	G4cout<<" 	1 - Synrad+ (Gamma,Debye-Waller)"<<G4endl;
	G4cout<<" 	2 - Synrad+ (Gamma,perturb norm - old model)"<<G4endl;
	G4cout<<" 	3 - Synrad+ (Gamma,Debye-Waller,perturb refl - new model)"<<G4endl;
	G4cout<<" Print surface material data = "<<GetReflPrint()<<G4endl;
	G4cout<<" 	0 - no"<<G4endl;
	G4cout<<" 	1 - yes"<<G4endl;
	G4cout<<" Sigma ratio = "<<GetSigmaRatio()<<G4endl;
	G4cout<<" Surface roughness = "<<GetSurfRoughness()/m<<" [m]"<<G4endl;
	G4cout<<" Surface autocorrelation length = "<<GetAutoCorrLength()/m<<" [m]"<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintGammaCut()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintGammaCut:"<<G4endl;
	G4cout<<" Ecut = "<<GetGammaEcut()/eV<<" [eV]"<<G4endl;
	G4cout<<" Tcut = "<<GetTrackTcut()/s<<" [s]"<<G4endl;
	G4cout<<" StraightTracksKillStatus = "<<GetStraightTracksKillStatus()<<G4endl;
	G4cout<<"  - Dmin = "<<GetStraightTracksDmin()/mm<<" [mm]"<<G4endl;
	G4cout<<"  - Zend = "<<GetStraightTracksZend()/cm<<" [cm]"<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintMeanFreePathFactor()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintMeanFreePathFactor ==> MeanFreePathFactor = "<<GetMeanFreePathFactor()<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

void BeamAndrii_SimParameters::PrintBeampipeFileName()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintBeampipeFileName ==> "<<GetBeampipeFileName()<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}

G4bool BeamAndrii_SimParameters::FileExists(const G4String& fName) 
{
  	struct stat buffer;   
  	return (stat (fName.c_str(), &buffer) == 0); 
}

void BeamAndrii_SimParameters::PrintIpBeampipeParameters()
{
	G4cout<<"\n\n=================================================================="<<G4endl;
	G4cout<<"[INFO] SimParameters::PrintIpBeampipeParameters ==> IP beam pipe parameters"<<G4endl;
	G4cout<<"FWD Z (zmin) = "<<GetIpBeampipeZmin()/m<<" [m]"<<G4endl;
	G4cout<<"BWD Z (zmax) = "<<GetIpBeampipeZmax()/m<<" [m]"<<G4endl;
	G4cout<<"==================================================================\n\n"<<G4endl;

	return;
}
