
//_____________________________________________________________________________
void training() {

  //two-dimensional Gaussian as a training function
  TF2 gauss2("gauss2","exp(-(x^2)-(y^2))", -3, 3, -3, 3);

  //program options from boost for LPS configuration
  using namespace boost;
  program_options::options_description opt("opt");

  //LPS instance, solution dimensionality is 1 (single value)
  LookupProblemSolver solver(1, "Gauss2D", &opt);

  //independent quantities, x and y of the Gaussian
  solver.MakeQuantity("x");
  solver.MakeQuantity("y");

  //load 'conf.ini' with optional configuration for the layers
  program_options::variables_map opt_map;
  program_options::store(parse_config_file("conf.ini", opt), opt_map);

  //optional non-uniform segmentation in quantities
  //solver.SetUseMedCellFinder();

  //initialize the LPS for training
  solver.Initialize(opt_map);

  //linear function to generate values for x and y,
  //non-uniform from 0.3 to 1 over x from -3 to 3
  TF1 linear("linear", "(1.3/2)+(0.7/6)*x", -3, 3);

  //number of training values
  int nxy = 10000;

  //training loop
  for(int i=0; i<nxy; i++) {

    //x and y as independent quantities
    double x = linear.GetRandom();
    double y = linear.GetRandom();

    //Gaussian value as known solution
    double gv = gauss2.Eval(x, y);

    //set the quantities and solution to vector containers for LPS
    vector<double> quant = {x, y}; // dimensionality (length) is 2
    vector<double> sol = {gv}; // length is 1

    //add training input to LPS
    solver.AddInput(quant, sol);
  }

  //finalize the LPS for export
  solver.Finalize();

  //export to file
  TFile out("lps_gaussian.root", "recreate");

  solver.Export();

  out.Close();



  //draw the Gaussian which was used for the training
  TCanvas can("can", "can", 768, 768);

  gauss2.Draw("colz");

  can.SaveAs("gauss.pdf");












}//training


