
void run_geo_s() {

  gROOT->ProcessLine(".include ../../../../calorimeters/include");

  gROOT->ProcessLine(".L ../../../../calorimeters/src/FiberYZ.cxx");
  gROOT->ProcessLine(".x run_geo.C");

}

