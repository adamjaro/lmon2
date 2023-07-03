
//C++ headers
#include <iostream>
#include <string>
#include <boost/program_options.hpp>

//explicit instantiation to reduce warnings from old Geant 4.10.00 on EIC nodes
//template std::vector<std::basic_string<char> > boost::program_options::to_internal(const std::vector<std::basic_string<char> >&);

//Geant headers
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"

//local headers
#include "DetectorConstruction.h"
#include "ActionInitialization.h"

namespace po = boost::program_options;
using namespace std;

//_____________________________________________________________________________
int main(int argc, char* argv[]) {

  //command line arguments
  po::options_description opt("Program arguments");
  opt.add_options()("vis", po::value<string>(), "visualization macro");
  opt.add_options()("mac,m", po::value<string>(), "batch macro");
  opt.add_options()("seed,s", po::value<long>(), "random seed");
  opt.add_options()("phys", po::value<string>(), "physics list");
  opt.add_options()("optics", po::value<bool>(), "turn optics on");

  //parse the arguments
  po::variables_map opt_map;
  po::store(po::parse_command_line(argc, argv, opt), opt_map);

  //visualization macro
  G4String vis_mac("visualize.mac"); // default name
  if(opt_map.count("vis")) {
    //name from the argument
    vis_mac = G4String( opt_map["vis"].as<string>() );
  }
  //batch macro
  G4String macro;
  if(opt_map.count("mac")) {
    macro = G4String( opt_map["mac"].as<string>() );
  }
  //random seed
  long random_seed = 0;
  if(opt_map.count("seed")) {
    random_seed = opt_map["seed"].as<long>();
  }
  //physics_list
  G4String physics_name("FTFP_BERT");
  if(opt_map.count("phys")) {
    physics_name = G4String( opt_map["phys"].as<string>() );
  }

  //UI session
  G4String session;
  G4UIExecutive *ui = 0x0;
  if( !macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }

  //random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  //random seed
  G4cout << "Input seed: " << random_seed << G4endl;;
  G4Random::setTheSeed(random_seed);
  G4cout << "Start seed: " << G4Random::getTheSeed() << G4endl;

  //default run manager
  G4RunManager *runManager = new G4RunManager;

  //detector construction
  runManager->SetUserInitialization(new DetectorConstruction);

  //physics
  G4PhysListFactory phys_factory;
  G4VModularPhysicsList *physicsList = phys_factory.GetReferencePhysList(physics_name);

  //turn optics on if requested
  if(opt_map.count("optics") and opt_map["optics"].as<bool>() == true) {
    G4OpticalPhysics *optics = new G4OpticalPhysics();
    physicsList->RegisterPhysics(optics);
  }

  //set the physics list
  runManager->SetUserInitialization(physicsList);

  //action
  runManager->SetUserInitialization(new ActionInitialization);

  //user interface manager
  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  //UI session for interactive or batch mode
  G4String command = "/control/execute ";
  if ( macro.size() > 0 ) {
    // batch mode
    UImanager->ApplyCommand(command+macro);
  } else {
    //interactive
    G4VisExecutive *visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand(command+vis_mac);
    ui->SessionStart();
    delete ui;
  }

  //job termination
  delete runManager;

  return 0;

}//main





























