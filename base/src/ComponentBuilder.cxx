
//_____________________________________________________________________________
//
// construction of individual detectors and components
// is here
//_____________________________________________________________________________


//ROOT
#include "Rtypes.h"

//Geant
#include "G4LogicalVolume.hh"

//local classes
#include "Detector.h"
#include "GeoParser.h"
#include "ComponentBuilder.h"
#include "../../calorimeters/include/CaloBuilder.h"
#include "../../beamline/include/BeamBuilder.h"
#include "../../trackers/include/TrackerBuilder.h"

using namespace std;

//_____________________________________________________________________________
ComponentBuilder::ComponentBuilder(G4LogicalVolume *top, GeoParser *geo, std::vector<Detector*> *vdet) {

  //register the individual builders
  fBuild.push_back( new CaloBuilder(top, geo) );
  fBuild.push_back( new BeamBuilder(top, geo) );
  fBuild.push_back( new TrackerBuilder(top, geo) );

  //geometry loop
  for(unsigned int i=0; i<geo->GetN(); i++) {

    //detector type and name
    G4String type = geo->GetType(i);
    G4String name = geo->GetName(i);

    //construct detector or component of type 'type'
    Detector *det = 0x0;

    //builder loop
    for(vector<BuilderBase*>::iterator ib = fBuild.begin(); ib != fBuild.end(); ib++) {

      //load the detector from the builder if present
      det = (*ib)->FindAndLoad(type, name);

      //add detector to all detectors
      if( det ) det->Add(vdet);      

    }//builder loop
  }//geometry loop

}//ComponentBuilder



















