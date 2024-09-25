
#ifndef GeoParser_h
#define GeoParser_h

//geometry data for detector and component construction

#include <map>
#include <boost/tokenizer.hpp>
#include "G4String.hh"

class GeoParser {

  public:

    void AddNew(G4String type, G4String name);
    void AddPar(G4String pnam, G4String val);
    void AddConst(G4String name, G4String value);

    //detectors and elements
    unsigned int GetN() const {return fDet.size();}
    const G4String& GetType(unsigned int i) const {return fDet[i].first;}
    const G4String& GetName(unsigned int i) const {return fDet[i].second;}

    //top volume name
    G4String GetTopName();

    //geometry parameters
    G4String GetS(G4String name, G4String par);
    G4double GetD(G4String name, G4String par);
    G4int GetI(G4String name, G4String par);
    G4bool GetB(G4String name, G4String par);

    //units for optional parameters
    class Unit {
      public:
        Unit(): apply(false), u(0) {}
        Unit(G4double uu): apply(true), u(uu) {}
      private:
        bool apply; // flag to indicate that the unit was provided
        G4double u; // unit for the parameter
        friend class GeoParser;
    };

    //optional parameters
    G4bool GetOptD(G4String name, G4String par, G4double& val, const Unit& un = Unit());
    G4bool GetOptI(G4String name, G4String par, G4int& val);
    G4bool GetOptB(G4String name, G4String par, G4bool& val);
    G4bool GetOptS(G4String name, G4String par, G4String& val);
    G4bool HasParameter(G4String name, G4String par);

    //constants for development
    G4String GetConst(std::string nam);

    void PrintAll();

  private:

    typedef boost::tokenizer< boost::char_separator<char> >::iterator token_it;

    template<typename T> T GetNum(G4String name, G4String par);
    template<typename T> bool GetOptNum(G4String name, G4String par, T& val);

    std::vector< std::pair<G4String, G4String> > fDet; // detectors and components

    //constant for geometry parameters
    class Constant {
    public:

      int type;
      std::string val_s;
      double val_d;

      enum con_type {kString=0, kDouble}; // string or numerical type

      Constant(): type(kString), val_s(""), val_d(0) {}

      //constant value as numeric type
      template<typename T> T GetNum() {

        T x;

        if( type == kString ) {
          std::istringstream st(val_s);
          st >> x;
        } else {
          x = val_d;
        }

        return x;

      }//GetNum

      //constant value as string
      std::string GetS() {

        std::string x;

        if( type == kString ) {
          x = val_s;
        } else {
          x = std::to_string(val_d);
        }

        return x;

      }//GetS

    };//Constant

    Constant Evaluate(G4String val);

    std::map<G4String, Constant> fConst; // constants for geometry

    std::map<G4String, Constant> fPar; // geometry parameters

};

#endif



















