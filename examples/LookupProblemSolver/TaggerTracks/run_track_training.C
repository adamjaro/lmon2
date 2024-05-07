
//_____________________________________________________________________________
void run_track_training() {

  //load shared library for the LPS
  gSystem->Load("../LPS_library/build/libLPS.so");

  //set ROOT environment to work with LPS

  //boost
  gInterpreter->Declare("#include <boost/program_options.hpp>");

  //header for LPS
  gInterpreter->AddIncludePath("../../../reconstruction/include/");
  gInterpreter->Declare("#include \"LookupProblemSolver.h\"");

  //header for tracks representation
  gInterpreter->Declare("#include \"TaggerTracks.h\"");

  //run the training macro
  gROOT->ProcessLine(".x track_training.C");

}//run_track_training




