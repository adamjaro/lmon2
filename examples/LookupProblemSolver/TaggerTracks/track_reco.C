
//_____________________________________________________________________________
void track_reco() {

  //tracks physics (photoproduction) input for kinematics reconstruction
  TFile input("trk_physics.root", "read");

  //input tree
  TTree *event_tree = dynamic_cast<TTree*>( input.Get("event") );

  //attach true kinematics to compare with solutions by LPS
  Double_t true_en, true_theta, true_phi;
  event_tree->SetBranchAddress("true_el_E", &true_en);
  event_tree->SetBranchAddress("true_el_theta", &true_theta);
  event_tree->SetBranchAddress("true_el_phi", &true_phi);

  //input tracks
  TaggerTracks::Coll tracks;
  tracks.ConnectInput("s1_tracks", event_tree);


  //LPS instance for obtaining solutions
  LookupProblemSolver solver(3, "Tagger_1");

  //import the trained LPS
  TFile inp_lps("lps_tracks.root", "read");
  solver.Import(inp_lps);


  //output file and tree on track kinematics by LPS solutions

  //create the file and tree
  TFile out("rec_tracks.root", "recreate");
  TTree otree("rec_trk", "rec_trk");

  //branches with reconstructed kinematics by LPS
  Double_t rec_en, rec_theta, rec_phi;
  otree.Branch("rec_en", &rec_en, "rec_en/D");
  otree.Branch("rec_theta", &rec_theta, "rec_theta/D");
  otree.Branch("rec_phi", &rec_phi, "rec_phi/D");

  //branches with true kinematics for comparison
  otree.Branch("true_en", &true_en, "true_en/D");
  otree.Branch("true_theta", &true_theta, "true_theta/D");
  otree.Branch("true_phi", &true_phi, "true_phi/D");

  //event loop
  Long64_t nev = event_tree->GetEntries();
  for(Long64_t iev=0; iev<nev; iev++) {

    //get the event and load the tracks
    event_tree->GetEntry(iev);
    tracks.LoadInput();

    //tracks loop
    for(const TaggerTracks::Track& trk: tracks.GetReadData()) {

      //independent quantities for LPS from track geometry
      vector<double> quant = {trk.x, trk.y, trk.theta_x, trk.theta_y};

      //prepare the vector container for the solutions
      vector<double> sol;
      sol.resize(3); // length is 3

      //obtain the solution and check its status (passed to sol)
      bool stat = solver.Solve(quant, sol);
      if( !stat ) continue;

      //put the obtained solution to output tree
      rec_en = sol[0];
      rec_theta = sol[1];
      rec_phi = sol[2];

      //fill the output tree (also true kinematics will be there)
      otree.Fill();
    }
  }

  //write the output tree and close the output file
  otree.Write(0, TObject::kOverwrite);

  out.Close();

  cout << "All done" << endl;

















}//track_reco



