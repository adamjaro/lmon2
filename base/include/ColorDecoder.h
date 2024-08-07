
#ifndef ColorDecoder_h
#define ColorDecoder_h

class GeoParser;

class ColorDecoder {

  public:

    ColorDecoder(G4String col): fCol(col) {}

    void SetCol(G4String col) { fCol = col; }
    void SetAuxEdgeVisible(G4bool v) { fAuxEdge = v; }

    G4VisAttributes* MakeVis(GeoParser *geo, G4String nam, G4String par);

  private:

    G4String fCol; // default visibility, red:green:blue:alpha
    G4bool fAuxEdge=true; // flag for aux edge visible

};

#endif

