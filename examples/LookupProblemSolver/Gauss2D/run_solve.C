
//_____________________________________________________________________________
void run_solve() {

  //load shared library for the LPS
  gSystem->Load("../LPS_library/build/libLPS.so");

  //set ROOT environment to work with LPS

  //boost
  gInterpreter->Declare("#include <boost/program_options.hpp>");

  //header for LPS
  gInterpreter->AddIncludePath("../../../reconstruction/include/");
  gInterpreter->Declare("#include \"LookupProblemSolver.h\"");

  gROOT->ProcessLine(".x solve.C");

}//run_solve

