
//_____________________________________________________________________________
//
// Provides variable binning for input TH1D histogram based on median cell
// split
// 
// Usage: construct with a TH1D input histogram (fixed binning) and number
// of cells and obtain an empty histogram with variable binning
// by call to MakeH1D.
//
// Other public function are present for development purposes.
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include <vector>
#include <map>

//ROOT
#include "TH1D.h"
#include "TF1.h"

//local classes
#include "MedCellFinder.h"

using namespace std;

//_____________________________________________________________________________
MedCellFinder::MedCellFinder(TH1D *hx, std::size_t ncells): fInpH1(hx) {

  cout << "MedCellFinder: " << hx->GetName() << endl;

  //input function from hx
  TF1 fx("f_inp", this, &MedCellFinder::EvalInpH1, 0, 1, 0);
  fFuncInpH1 = fx;

  fFuncInpH1.SetNpx(2000);

  //initialize the cells
  fCells.resize(ncells);

  fCells[0].xmin = hx->GetBinLowEdge(1);
  fCells[0].xmax = hx->GetBinLowEdge(hx->GetNbinsX())+hx->GetBinWidth(hx->GetNbinsX());
  fCells[0].is_active = true;

  //cell split loop
  while( fLastCell+2 < fCells.size() ) {

    //active cell with largest integral
    Double_t cell_i = -1;
    size_t idx = 0;

    for(size_t i=0; i<fLastCell+1; i++) {
      Cell& c = fCells[i];
      if( !c.is_active ) continue;

      Double_t in = fFuncInpH1.Integral(c.xmin, c.xmax);
      if( in < cell_i ) continue;

      cell_i = in;
      idx = i;
    }

    //split the cell with largest integral
    Cell& c_split = fCells[idx];
    c_split.is_active = false;
    fFuncInpH1.SetRange(c_split.xmin, c_split.xmax);

    Double_t med[1], quantile[1]={0.5}; 
    fFuncInpH1.GetQuantiles(1, med, quantile);

    //cout << idx << " " << cell_i << " " << med[0] << endl;

    //two new cells below and above the medial
    Cell& c_lo = fCells[fLastCell+1];
    Cell& c_hi = fCells[fLastCell+2];

    c_lo.xmin = c_split.xmin;
    c_lo.xmax = med[0];
    c_lo.is_active = true;

    c_hi.xmin = med[0];
    c_hi.xmax = c_split.xmax;
    c_hi.is_active = true;

    fLastCell += 2;

  }//cell split loop

  //indices for active cells, sorted by xmin
  map<Double_t, size_t> pos_idx;

  //active cells loop
  for(size_t i=0; i<fLastCell+1; i++) {
    const Cell& c = fCells[i];
    if( !c.is_active ) continue;

    pos_idx.insert( make_pair(c.xmin, i) );
  }//active cells loop

  //map loop
  for(const pair<Double_t, size_t>& i: pos_idx) {

    fActiveCellIdx.push_back(i.second);

  }//map loop

}//MedCellFinder

//_____________________________________________________________________________
Double_t MedCellFinder::EvalInpH1(Double_t *x, Double_t*) {

  Double_t interp = fInpH1->Interpolate(x[0]);

  if( interp < 0 ) { return 0; }

  return interp;

}//EvalInpH1

//_____________________________________________________________________________
TH1D MedCellFinder::MakeH1D(const char *name_title) {

  vector<Double_t> bins;

  //active cells loop
  for(const size_t& i: fActiveCellIdx) {
    const Cell& c = fCells[i];

    bins.push_back( c.xmin );
  }//active cells loop

  //end of last bin
  bins.push_back( fCells[fActiveCellIdx.back()].xmax );

  TH1D hx(name_title, name_title, bins.size()-1, bins.data());

  return hx;

}//MakeH1D

//_____________________________________________________________________________
void MedCellFinder::PrintCells() {

  cout << "MedCellFinder::PrintCells" << endl;

  for(size_t i=0; i<fLastCell+1; i++) {
    const Cell& c = fCells[i];

    cout << i << " " << c.is_active << " " << c.xmin << " " << c.xmax << endl;
  }

  cout << "Active only:" << endl;
  for(const size_t& i: fActiveCellIdx) {
    const Cell& c = fCells[i];

    cout << i << " " << c.xmin << " " << c.xmax << endl;
  }

}//PrintCells

















