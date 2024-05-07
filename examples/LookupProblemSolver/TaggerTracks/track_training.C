
//_____________________________________________________________________________
void track_training() {

  //tracks input for LPS training
  TFile input("trk_training.root", "read");

  //input tree with tracks and true kinematics
  TTree *event_tree = dynamic_cast<TTree*>( input.Get("event") );

  //attach the true kinematics (electron energy and theta and phi angles)
  Double_t true_en, true_theta, true_phi;
  event_tree->SetBranchAddress("true_el_E", &true_en);
  event_tree->SetBranchAddress("true_el_theta", &true_theta);
  event_tree->SetBranchAddress("true_el_phi", &true_phi);

  //input tracks
  TaggerTracks::Coll tracks;
  tracks.ConnectInput("s1_tracks", event_tree);

  //program options from boost for LPS configuration
  using namespace boost;
  program_options::options_description opt("opt");

  //LPS instance, solution dimensionality is 3 (energy, theta, phi)
  LookupProblemSolver solver(3, "Tagger_1", &opt);

  //independent quantities as track geometry parameters
  solver.MakeQuantity("x"); // x position, mm
  solver.MakeQuantity("y"); // y position, mm
  solver.MakeQuantity("theta_x"); // theta_x angle, rad
  solver.MakeQuantity("theta_y"); // theta_y angle, rad

  //load 'conf.ini' with optional configuration for the layers
  program_options::variables_map opt_map;
  program_options::store(parse_config_file("conf.ini", opt), opt_map);

  //optional non-uniform segmentation in quantities
  //solver.SetUseMedCellFinder();

  //initialize the LPS for training
  solver.Initialize(opt_map);

  //input event loop
  Long64_t nev = event_tree->GetEntries();
  for(Long64_t iev=0; iev<nev; iev++) {

    //load the current event
    event_tree->GetEntry(iev);
    tracks.LoadInput();

    //select events with one track
    if( tracks.GetN() != 1 ) continue;

    //get the track (it is known there is one)
    const TaggerTracks::Track& trk = tracks.GetUnit(0);

    //set independent quantitis for LPS from track geometry
    vector<double> quant = {trk.x, trk.y, trk.theta_x, trk.theta_y};

    //set the known solution from true kinematics, length is 3
    vector<double> sol = {true_en, true_theta, true_phi};

    //add training input to LPS
    solver.AddInput(quant, sol);
  }

  //finalize the LPS for export
  solver.Finalize();

  //export to file
  TFile out("lps_tracks.root", "recreate");

  solver.Export();

  out.Close();

  cout << "All done" << endl;





}//track_training


