
//_____________________________________________________________________________
//
// Geometry from XML
//
//_____________________________________________________________________________

//C++
#include <vector>
#include <map>
#include <fstream>

//Boost
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/detail/xml_parser_utils.hpp>

//Geant
#include "G4ios.hh"

//local classes
#include "LoadXML.h"

using namespace std;
using namespace boost::property_tree;

//_____________________________________________________________________________
void LoadXML::ReadInput(string inp) {

  //read the XML in inp

  G4cout << "LoadXML::ReadInput, " << inp << G4endl; 

  //open the input
  ifstream in(inp);

  //load the xml to property_tree
  ptree tree;
  read_xml(in, tree);

  //top branch named main
  ptree& main = tree.get_child("main");

  //recursive read through branches
  ReadPTree("main", main);

}//ReadInput

//_____________________________________________________________________________
void LoadXML::ReadPTree(const ptree::key_type& key, const ptree& tree) {

  //recursive read through branches

  if( key == "comment" ) return;

  //next subbranch
  const ptree& tree_child = tree.get_child_optional( xml_parser::xmlattr<ptree::key_type>() ).get();

  //test individual load functions for attributes for the subbranch
  LoadDetector(key, tree_child); //'detector' token, next detector or component
  LoadInclude(key, tree_child); //'include' token, include command
  LoadIncludeDir(key, tree_child); //'include_directory' token, path to include files
  LoadConstant(key, tree_child); //'const' token, constant in geometry

  //call for next subbranch
  for(ptree::const_iterator it = tree.begin(); it != tree.end(); it++) {

    //skip attributes (handled by load functions), text and comments
    if( it->first == xml_parser::xmlattr<ptree::key_type>() ) continue;
    if( it->first == xml_parser::xmltext<ptree::key_type>() ) continue;
    if( it->first == xml_parser::xmlcomment<ptree::key_type>() ) continue;

    ReadPTree(it->first, it->second); //next subbranch
  }

  //finish the detector after all recursive reads for the detector
  if( key == "detector" ) {
    FinishDetector();
  }

}//ReadPTree

//_____________________________________________________________________________
void LoadXML::LoadDetector(const ptree::key_type& key, const ptree& tree) {

  if( key == "detector" ) {

    //call directly for detector
    fDetectorLoad = true;

  } else if( fDetectorLoad == true ) {

    //call for subbranch in detector

    //counter for nested attributes in the detector
    if( fDetNestedAttrCount.find(key) == fDetNestedAttrCount.end() ) {
      fDetNestedAttrCount[key] = 0;
    } else {
      fDetNestedAttrCount[key]++;
    }
  }

  //test for loading the detector
  if( !fDetectorLoad ) return;

  //detector attributes loop
  for(ptree::const_iterator it = tree.begin(); it != tree.end(); it++) {

    //parameter name from attributes
    string par_name = it->first;

    //direct attribute or nested branch for the detector
    if( key == "detector" ) {

      //direct attribute for the detector
      fDetParam.push_back( make_pair(par_name, it->second.get_value<string>()) );
    } else {

      //nested branch in detector, order in the vector by enum:
      //enum kdet_nes {knes_key=0, knes_cnt, knes_nam, knes_val}
      fDetNestedAttr.push_back({key, to_string(fDetNestedAttrCount[key]), par_name, it->second.get_value<string>()});
    }
  }//detector attributes loop

}//LoadDetector

//_____________________________________________________________________________
void LoadXML::FinishDetector() {

  //get type and name for the detector to be finished
  string type, name;

  for(pair<string, string>& i: fDetParam) {
    if( i.first == "type" ) type = i.second;
    if( i.first == "name" ) name = i.second;
  }

  //add new detector to the geometry
  fGeo.AddNew(type, name);

  //direct detector attributes
  for(pair<string, string>& i: fDetParam) {
    if( i.first == "type" or i.first == "name" ) continue;

    //add the attribute to the geometry
    fGeo.AddPar(name+"."+i.first, i.second);
  }

  //nested detector attributes
  for(vector<string>& i: fDetNestedAttr) {

    //one or more nested branches of the same key (branch name)
    if( fDetNestedAttrCount[i[knes_key]] <= 0 ) {

      //just one nested branch of the given key, add directly to the geometry
      fGeo.AddPar(name+"."+i[knes_key]+"_"+i[knes_nam], i[knes_val]);
    } else {

      //more nested branches of the same key, append the counter to add to geometry
      fGeo.AddPar(name+"."+i[knes_key]+i[knes_cnt]+":"+i[knes_nam], i[knes_val]);
    }
  }

  //clear the structures for detector parameters, ready for the next detector
  fDetParam.clear();
  fDetNestedAttr.clear();
  fDetNestedAttrCount.clear();

  //clear the flag for detector reading
  fDetectorLoad = false;

}//FinishDetector

//_____________________________________________________________________________
void LoadXML::LoadInclude(const ptree::key_type& key, const ptree& tree) {

  //include another xml

  if( fDetectorLoad or key != "include" ) return;

  //ptree loop
  for(ptree::const_iterator it = tree.begin(); it != tree.end(); it++) {

    //look for ref key
    if( it->first != "ref" ) continue;

    //file name to be included
    string include_name = it->second.get_value<string>();

    //directory for include files, if set
    if( !fIncludeDir.empty() ) {
      include_name = fIncludeDir + "/" + include_name;
    }

    //load the included xml
    ReadInput( include_name );

  }//ptree loop

}//LoadInclude

//_____________________________________________________________________________
void LoadXML::LoadIncludeDir(const ptree::key_type& key, const ptree& tree) {

  //directory for include files

  if( fDetectorLoad or key != "include_directory" ) return;

  for(ptree::const_iterator it = tree.begin(); it != tree.end(); it++) {

    if( it->first == "name" ) fIncludeDir = it->second.get_value<string>();
  }

}//LoadIncludeDir

//_____________________________________________________________________________
void LoadXML::LoadConstant(const ptree::key_type& key, const ptree& tree) {

  if( fDetectorLoad or key != "const" ) return;

  //constant name and value
  string name, value;

  for(ptree::const_iterator it = tree.begin(); it != tree.end(); it++) {

    if( it->first == "name" ) name = it->second.get_value<string>();
    if( it->first == "value" ) value = it->second.get_value<string>();
  }

  //add the constant to the geometry
  fGeo.AddConst(name, value);

}//LoadConstant


















