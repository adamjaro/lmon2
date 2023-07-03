
//_____________________________________________________________________________
//
// Structure holding geometry data for detector and component construction
//
//_____________________________________________________________________________

//C++
#include <fstream>
#include <string>
#include <boost/tokenizer.hpp>

//ROOT
#include "TFormula.h"

//Geant
#include "G4String.hh"
#include "G4ios.hh"
#include "globals.hh"

//local classes
#include "GeoParser.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
void GeoParser::AddNew(G4String type, G4String name) {

  //new detector

  fDet.push_back( make_pair(type, name) );

}//AddNew

//_____________________________________________________________________________
void GeoParser::AddPar(G4String pnam, G4String val) {

  //geometry parameter 'par' for detector 'det' as pnam = det.par

  fPar.insert( make_pair(pnam, Evaluate(val)) );

}//AddPar

//_____________________________________________________________________________
void GeoParser::AddConst(G4String name, G4String value) {

  //add new constant

  //name and the constant from Evaluate function
  fConst.insert( make_pair(name, Evaluate(value)) );

}//AddConst

//_____________________________________________________________________________
GeoParser::Constant GeoParser::Evaluate(G4String val) {

  //initialize the constant to be returned
  Constant c;
  c.type = Constant::kString;
  c.val_s = val;

  //arithmetic parameters for possible expression
  map<string, double> arithmetic_pars;

  //parse for possible constants in the constant value
  tokenizer< char_separator<char> > val_sep(val, char_separator<char>("", "+*-/()"));
  int ntok = 0;
  stringstream ss;

  //token loop
  for(token_it i = val_sep.begin(); i != val_sep.end(); i++) {

    //current token
    string tok = *i;

    //remove spaces
    tok.erase( remove_if(tok.begin(), tok.end(), ::isspace), tok.end() );

    //possible constant for the token
    map<G4String, Constant>::iterator iconst = fConst.find(tok);

    //no constant found, use the token as it is
    if( iconst == fConst.end() ) {
      ss << tok;
      ntok++;
      continue;
    }

    //constant is present for the token, substitute depending on type
    Constant& c_subst = iconst->second;

    //substitute for the case of just one token present
    c = c_subst;

    //set the arithmetic expression for the case of more tokens
    if( c_subst.type == Constant::kString ) {

      //string constant, set its string value
      ss << c_subst.val_s;

    } else {

      //double constant, set as the parameter
      ss << "[" << iconst->first << "]"; // constant name for parameter syntax in TFormula
      arithmetic_pars.insert( make_pair(iconst->first, c_subst.val_d) );
    }

    ntok++;

  }//token loop

  //arithmetic expression for the value of the constant, calculate the value and set type to double
  if(ntok > 1) {
    TFormula form("form", ss.str().c_str(), false);

    //parameter loop
    for( pair<const string, double> i: arithmetic_pars ) {

      form.AddParameter(i.first, i.second);

    }//parameter loop

    //constant value from the formula
    c.type = Constant::kDouble; // set type to double
    c.val_d = form.Eval(0); // set the double value

    return c;
  }

  //string constant or single substitution in the token loop
  return c;

}//Evaluate

//_____________________________________________________________________________
G4String GeoParser::GetTopName() {

  //name of the top volume

  vector< pair<G4String, G4String> >::reverse_iterator i = fDet.rbegin();
  while(i != fDet.rend()) {

    if( (*i).first == "top" ) return (*i).second;

    i++;
  }

  return "";

}//GetTopName

//_____________________________________________________________________________
G4String GeoParser::GetS(G4String name, G4String par) {

  //load geometry parameter 'par' as a string from the map for detector 'name'

  map<G4String, Constant>::iterator i = fPar.find(name+"."+par);

  if( i == fPar.end() ) {
    //parameter not found
    string description = "Parameter '" + par + "' not found for '" + name + "'";
    G4Exception("GeoParser::GetS", "ParameterNotFound01", FatalException, description.c_str());
  }

  return i->second.GetS();

}//GetS

//_____________________________________________________________________________
template<typename T> T GeoParser::GetNum(G4String name, G4String par) {

  //get parameter 'par' value as numeric type T for detector named 'name'

  map<G4String, Constant>::iterator i = fPar.find(name+"."+par);

  if( i == fPar.end() ) {
    string description = "Parameter '" + par + "' not found for '" + name + "'";
    G4Exception("GeoParser::GetNum", "ParameterNotFound01", FatalException, description.c_str());
  }

  return i->second.GetNum<T>();

}//GetNum

//_____________________________________________________________________________
G4double GeoParser::GetD(G4String name, G4String par) {

  //parameter 'par' for detector 'name' as G4double

  return GetNum<G4double>(name, par);

}//GetD

//_____________________________________________________________________________
G4int GeoParser::GetI(G4String name, G4String par) {

  //parameter 'par' for detector 'name' as G4int

  return GetNum<G4int>(name, par);

}//GetI

//_____________________________________________________________________________
G4bool GeoParser::GetB(G4String name, G4String par) {

  //parameter 'par' for detector 'name' as G4bool

  return GetNum<G4bool>(name, par);

}//GetB

//_____________________________________________________________________________
template<typename T> bool GeoParser::GetOptNum(G4String name, G4String par, T& val) {

  //get value 'val' for optional parameter 'par' as T for detector named 'name'

  map<G4String, Constant>::iterator i = fPar.find(name+"."+par);
  if( i == fPar.end() ) return false;

  val = i->second.GetNum<T>();

  return true;

}//GetOptPar

//_____________________________________________________________________________
G4bool GeoParser::GetOptD(G4String name, G4String par, G4double& val, const Unit& un) {

  //optional G4double parameter

  if( !GetOptNum<G4double>(name, par, val) ) return false;

  //apply the units if provided
  if(un.apply) val *= un.u;

  return true;

}//GetOptD

//_____________________________________________________________________________
G4bool GeoParser::GetOptI(G4String name, G4String par, G4int& val) {

  //optional G4int parameter

  return GetOptNum<G4int>(name, par, val);

}//GetOptI

//_____________________________________________________________________________
G4bool GeoParser::GetOptB(G4String name, G4String par, G4bool& val) {

  //optional G4bool parameter

  return GetOptNum<G4bool>(name, par, val);

}//GetOptI

//_____________________________________________________________________________
G4bool GeoParser::GetOptS(G4String name, G4String par, G4String& val) {

  //optional G4String parameter

  map<G4String, Constant>::iterator i = fPar.find(name+"."+par);

  if( i == fPar.end() ) return false;

  val = i->second.GetS();

  return true;

}//GetOptS

//_____________________________________________________________________________
G4String GeoParser::GetConst(std::string nam) {

  //retrieve the value of the constant for geometry development

  map<G4String, Constant>::iterator iconst = fConst.find(nam);
  if( iconst == fConst.end() ) {
    return "";
  }

  return iconst->second.GetS();

}//GetConst

//_____________________________________________________________________________
void GeoParser::PrintAll() {

  //print all geometry data

  G4cout << "Detectors:" << G4endl;

  for(pair<G4String, G4String>& i: fDet) {

    G4cout << i.first << " " << i.second << G4endl;
  }

  G4cout << "Parameters:" << G4endl;

  for(pair<const G4String, Constant>& i: fPar) {

    G4cout << i.first << " " << i.second.type << " " << i.second.GetS() << G4endl;
  }


  G4cout << "Constants, name-type-num-str:" << G4endl;

  for(pair<const G4String, Constant>& i: fConst) {

    G4cout << i.first << " " << i.second.type << " " << i.second.GetNum<double>() << " " << i.second.GetS() << G4endl;
  }

}//PrintAll











