
//_____________________________________________________________________________
//
// construction of beamline components
//
//_____________________________________________________________________________

//C++
#include <vector>

//ROOT
//#include "Rtypes.h"

//Geant
#include "G4LogicalVolume.hh"

//local classes
#include "GeoParser.h"
#include "BeamBuilder.h"

//local detectors and components
#include "BeamQuadrupole.h"
#include "BeamDipole.h"
#include "CylBeam.h"
#include "ConeBeam.h"
#include "CylSegment.h"
#include "BeamDrift.h"
#include "BoxSegment.h"
#include "ParticleCounter.h"
#include "VacDrift.h"
#include "BeamAndrii20240604.h"
#include "QDMagnet.h"
#include "ReadGDML.h"

//macros
#define ADD_COMPONENT(comp) (fComp.insert( make_pair(#comp, &BeamBuilder::MakeDet<comp>) ))
#define ADD_DETECTOR(det) (fDets.insert( make_pair(#det, &BeamBuilder::MakeDet<det>) ))

//_____________________________________________________________________________
BeamBuilder::BeamBuilder(G4LogicalVolume *top, GeoParser *geo): BuilderBase(),
  fTop(top), fGeo(geo) {

  //G4cout << "BeamBuilder::BeamBuilder" << G4endl;

  //passive components
  ADD_COMPONENT( BeamQuadrupole );
  ADD_COMPONENT( BeamDipole );
  ADD_COMPONENT( CylBeam );
  ADD_COMPONENT( ConeBeam );
  ADD_COMPONENT( CylSegment );
  ADD_COMPONENT( BeamDrift );
  ADD_COMPONENT( BoxSegment );
  ADD_COMPONENT( VacDrift );
  ADD_COMPONENT( BeamAndrii20240604 );
  ADD_COMPONENT( QDMagnet );
  ADD_COMPONENT( ReadGDML );

  //sensitive detectors
  ADD_DETECTOR( ParticleCounter );

}//BeamBuilder

//_____________________________________________________________________________
Detector* BeamBuilder::FindAndLoad(G4String type, G4String name) {

  //G4cout << "BeamBuilder::FindAndLoad, " << type << " " << name << G4endl;

  //factory component and detector construction
  Detector *det = 0x0; 
  std::map<G4String, MakeDetPtr>::iterator idet;

  //component
  idet = fComp.find(type);
  if( idet != fComp.end() ) {
    (this->*(*idet).second)(name, fGeo, fTop);
  }

  //detector
  idet = fDets.find(type);
  if( idet != fDets.end() ) {
    det = (this->*(*idet).second)(name, fGeo, fTop);
  }

  return det;

}//FindAndLoad


















