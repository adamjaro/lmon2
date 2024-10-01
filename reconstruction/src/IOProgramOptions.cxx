
//_____________________________________________________________________________
//
// Helper class for input/output from boost program_options
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include "glob.h"

//Boost
#include <boost/program_options.hpp>

//ROOT
#include "TChain.h"

#include "TFile.h"

//local classes
#include "IOProgramOptions.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
unique_ptr<TChain> IOProgramOptions::MakeChain(const string& inp, string tnam) {

  //inputs, remove quote marks if present
  string input = GetStr(inp);
  glob_t glob_inputs;
  glob(input.c_str(), GLOB_TILDE, NULL, &glob_inputs);

  //input tree
  unique_ptr<TChain> tree = make_unique<TChain>( tnam.c_str() );
  for(size_t i=0; i<glob_inputs.gl_pathc; i++) {
    cout << "IOProgramOptions::MakeChain, adding input: " << glob_inputs.gl_pathv[i] << endl;
    tree->Add( glob_inputs.gl_pathv[i] );
  }

  return tree;

}//MakeChain

//_____________________________________________________________________________
unique_ptr<TFile> IOProgramOptions::MakeFile(const std::string& out) {

  //output name from program_options, remove quotes in present
  string nam = GetStr(out);

  cout << "IOProgramOptions::MakeFile, " << nam << endl;

  //create the TFile
  return make_unique<TFile>(nam.c_str(), "recreate");

}//MakeFile

//_____________________________________________________________________________
string IOProgramOptions::GetStr(std::string par) {

  //string for a given program option

  string res = fMap[par].as<string>();
  res.erase(remove(res.begin(), res.end(), '\"'), res.end());

  return res;

}//GetStr






























