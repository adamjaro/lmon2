
//_____________________________________________________________________________
//
// Geometry in .in format, used with the first lmon version
//
//_____________________________________________________________________________

//C++
#include <fstream>
#include <string>
#include <boost/tokenizer.hpp>

//Geant
#include "G4String.hh"
#include "G4ios.hh"
#include "globals.hh"

//local classes
#include "LoadIN.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
void LoadIN::ReadInput(G4String input) {

  //load geometry input

  G4cout << "LoadIN::ReadInput, " << input << G4endl;

  ifstream in(input);

  if(in.fail()) {
    string description = "Can't open input: '" + input + "'";
    G4Exception("GeoParser::LoadInput", "InputNotOpen01", FatalException, description.c_str());
  }

  char_separator<char> sep(" ");
  string line;

  //input loop
  while( getline(in, line) ) {

    //remove leading white spaces
    if( line.find_first_not_of(" \t") != string::npos ) {
      line = line.substr( line.find_first_not_of(" \t") );
    }

    //skip comments and empty lines
    if(line.empty() || line.find("#") == 0) continue;

    //remove trailing comments
    line = line.substr(0, line.find_first_of("#"));

    tokenizer<char_separator<char>> cline(line, sep);
    tokenizer<char_separator<char>>::iterator it = cline.begin();

    //command at the given line
    string cmd = *it;

    if( cmd == "new" ) {
      //new element
      AddNew(it);

    } else if( cmd == "include" ) {
      //constant in geometry
      Include(it);

    } else if( cmd == "include_directory" ) {
      //constant in geometry
      fIncludeDir = *(++it);

    } else if( cmd == "const" ) {
      //constant in geometry
      AddConst(it, cline.end());

    } else if( cmd.find(".") != string::npos ) {
      //constant in geometry
      AddPar(it, cline.end());

    }

  }//input loop

  in.close();


}//ReadInput

//_____________________________________________________________________________
void LoadIN::Include(token_it &it) {

  //include another geometry input

  //file to include
  G4String infile = *(++it);

  //directory for include files, if set
  if( !fIncludeDir.empty() ) {
    infile = fIncludeDir + "/" + infile;
  }

  //read the included file
  ReadInput(infile);

}//Include

//_____________________________________________________________________________
void LoadIN::AddNew(token_it &it) {

  //add new element

  G4String type = *(++it);
  G4String name = *(++it);

  fGeo.AddNew(type, name);

}//AddNew

//_____________________________________________________________________________
void LoadIN::AddConst(token_it &it, token_it end) {

  //add new constant

  it++; // skip the const statement
  G4String name = *(it++); // constant name

  //constant value, spaces are allowed
  G4String value;
  while( it != end ) {
    value += *(it++);
  }

  fGeo.AddConst(name, value);

}//AddConst

//_____________________________________________________________________________
void LoadIN::AddPar(token_it &it, token_it end) {

  //add new geometry parameter

  G4String name = *(it++); // detector.parameter
  it++; // skip the equal sign

  //parameter value
  G4String value;
  while( it != end ) {
    value += *(it++);
  }

  fGeo.AddPar(name, value);

}//AddPar



















