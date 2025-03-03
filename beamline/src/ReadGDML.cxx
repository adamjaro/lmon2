//_____________________________________________________________________________
//
// Reader for GDML geometry
//
//_____________________________________________________________________________

//Geant
#include "G4LogicalVolume.hh"
#include "G4GDMLParser.hh"
#include "G4PVPlacement.hh"

//local classes
#include "GeoParser.h"
#include "ReadGDML.h"

using namespace std;

//_____________________________________________________________________________
ReadGDML::ReadGDML(const G4String& nam, GeoParser *geo, G4LogicalVolume *top):
    Detector(), fNam(nam) {

  G4cout << "ReadGDML: " << fNam << G4endl;

  G4GDMLParser parser;

  //parser.Read("/home/jaroslav/sim/Andrii/20250103/IR6_v6.3.1_Nov2024_2x50mm.gdml", false);
  //parser.Read("/home/jaroslav/sim/geant/geant4-v11.1.2/examples/extended/persistency/gdml/G01/axes.gdml");
  //parser.Read("/home/jaroslav/sim/geant/geant4-v11.1.2/examples/extended/persistency/gdml/G01/solids.gdml");
  //parser.Read("/home/jaroslav/sim/geant/geant4-v11.1.2/examples/extended/persistency/gdml/G01/scale.gdml");
  //parser.Read("/home/jaroslav/sim/geant/geant4-v11.1.2/examples/extended/persistency/gdml/G01/divisionvol.gdml");
  //parser.Read("/home/jaroslav/sim/geant/geant4-v11.1.2/examples/extended/persistency/gdml/G01/parameterized.gdml");
  //parser.Read("/home/jaroslav/sim/geant/geant4-v11.1.2/examples/extended/persistency/gdml/G01/pTube.gdml");
  //parser.Read("/home/jaroslav/sim/geant/geant4-v11.1.2/examples/extended/persistency/gdml/G01/auxiliary.gdml");

  parser.Read( geo->GetS(fNam, "gdml") );

  //G4cout << parser.GetWorldVolume()->GetName() << G4endl;

  new G4PVPlacement(0, G4ThreeVector(0, 0, 0), parser.GetWorldVolume()->GetLogicalVolume(), fNam, top, false, 0);

}//ReadGDML






















