
#ifndef LogicVolFinder_h
#define LogicVolFinder_h

class LogicVolFinder {

  public:

    static G4LogicalVolume* GetMotherVolume(G4String place_into, G4LogicalVolume *top, GeoParser *geo, G4String det_name);
    static G4LogicalVolume* GetMotherVolume2(G4String p1, G4String p2, G4LogicalVolume *top, GeoParser *geo, G4String det_name);

    static G4LogicalVolume* GetMotherVolume(G4String mother_nam, G4LogicalVolume *top);
    static G4LogicalVolume* GetMotherVolume2(G4String m1, G4String m2, G4LogicalVolume *top);

};

#endif

