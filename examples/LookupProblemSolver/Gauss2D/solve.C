
//_____________________________________________________________________________
void solve() {

  //LPS instance without program_options for obtaining solutions
  LookupProblemSolver solver(1, "Gauss2D"); // dimensionality and name match training

  //import trained LPS
  TFile inp("lps_gaussian.root", "read");
  solver.Import(inp);

  //number of points to obtain solutions
  int nxy = 1000;

  //graph to draw the solutions
  TGraph2D g2(nxy);

  //point loop
  for(int ix=0; ix<nxy; ix++) {

    //random x and y uniformly from -3 to 3
    double x = -3 + 6*gRandom->Rndm();
    double y = -3 + 6*gRandom->Rndm();

    //put x and y (independent quantities) to vector container
    vector<double> quant = {x, y}; // length is 2

    //prepare vector container for solution
    vector<double> sol;
    sol.resize(1); // set length as 1

    //obtain the solution for x and y (passed to 'sol')
    bool stat = solver.Solve(quant, sol);

    //check status if solution was found
    if( !stat ) continue;

    //mark the obtained solution in output graph
    g2.SetPoint(ix, x, y, sol[0]);

  }

  //draw the graph
  TCanvas can("can", "can", 768, 768);

  g2.Draw("surf2");

  can.SaveAs("gauss2.pdf");


}//solve






















