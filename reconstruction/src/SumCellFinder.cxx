
//_____________________________________________________________________________
//
// Provides variable binning for input TH1D histogram
// using even sum per bin from interpolated input histogram
// 
// Usage: the finder inherits from TH1D, creates a new one
// from the input histogram provided. Number of divisions over the full range
// of the input can be set by ndiv (1e6 by default)
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include <vector>

//ROOT
#include "TH1D.h"

//local classes
#include "SumCellFinder.h"

using namespace std;

//_____________________________________________________________________________
SumCellFinder::SumCellFinder(TH1D *hx, const char *name_title, int ndiv): TH1D() {
//SumCellFinder::SumCellFinder(shared_ptr<TH1D> hx, const char *name_title, int ndiv): TH1D() {

  cout << "SumCellFinder: " << hx->GetName() << endl;

  //number of bins from the original flat histogram
  UInt_t nbins = hx->GetNbinsX();

  //content integral per one bin
  Double_t nev_bin = hx->Integral("width")/nbins;
  //cout << "Events integral per bin: " << nev_bin << endl;

  //range from the original histogram
  Double_t xlow = hx->GetBinLowEdge(1);
  Double_t xup = hx->GetBinLowEdge(hx->GetNbinsX())+hx->GetBinWidth(hx->GetNbinsX());

  //division for increment integrals
  Double_t dx = (xup-xlow)/ndiv;

  //bins for the final histogram
  vector<Double_t> bins;
  bins.push_back( xlow );

  //finder loop
  while( bins.size() < nbins ) {

    Double_t xpos = bins.back();
    Double_t sum = 0;

    //increment integral in steps to reach original integral per bin
    while( sum < nev_bin and xpos < xup ) {
      xpos += dx;
      sum += dx*hx->Interpolate(xpos);
    }

    //new bin edge found
    bins.push_back( xpos );

    if( xpos > xup ) break;

  }//finder loop

  //last upper edge
  if( bins.back() < xup ) {
    bins.push_back( xup );
  }

  cout << "Bins: " << bins.size()-1 << ", range: " << xlow << " -> " << xup << endl;
  //for(size_t i=0; i<bins.size(); i++) {
    //cout << i << " " << bins[i] << endl;
    //cout << i+1 << ":" << bins[i] << " ";
  //}
  //cout << endl;

  //bins for TH1D
  SetBins(bins.size()-1, bins.data());

  //name and title for TH1D
  SetNameTitle(name_title, name_title);

}//SumCellFinder

















