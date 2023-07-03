
#ifndef LoadIN_h
#define LoadIN_h

// geometry in .in format

#include "GeoParser.h"

class LoadIN {

  public:

    LoadIN(GeoParser& geo): fGeo(geo) {}

    void ReadInput(G4String input);

  private:

    typedef boost::tokenizer<boost::char_separator<char>>::iterator token_it;

    void Include(token_it &it); // another geometry input
    void AddNew(token_it &it); // new detector or element
    void AddConst(token_it &it, token_it end); // new constant
    void AddPar(token_it &it, token_it end); // new geometry parameter

    GeoParser& fGeo;

    G4String fIncludeDir; // directory for include files, if set

};

#endif

