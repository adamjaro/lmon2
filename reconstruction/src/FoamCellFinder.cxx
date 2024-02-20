
//_____________________________________________________________________________
//
// Provides variable binning for input TH1D histogram based on TFoam
// 
// Usage: construct with a TH1D input histogram (fixed binning) and number
// of cells for FOAM and obtain an empty histogram with variable binning
// by call to MakeH1D.
//
// Other public function are present for development purposes.
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include <map>

//ROOT
#include "TFoamVect.h"
#include "TFoamCell.h"
#include "TH1D.h"
#include "TRandom3.h"

//local classes
#include "FoamCellFinder.h"

using namespace std;

//_____________________________________________________________________________
FoamCellFinder::FoamCellFinder(const TH1D& hx, Int_t ncells): TFoam(hx.GetName()) {

  cout << "FoamCellFinder: " << hx.GetName() << endl;

  fInpMin = hx.GetBinLowEdge(1);
  fInpMax = hx.GetBinLowEdge(hx.GetNbinsX())+hx.GetBinWidth(hx.GetNbinsX());

  //cout << fInpMin << " " << fInpMax << endl;

  fInpHist.SetNameTitle( (string(hx.GetName())+"_inp").c_str(), (string(hx.GetName())+"_inp").c_str() );
  fInpHist.SetBins(hx.GetNbinsX(), 0, 1);

  for(Int_t i=1; i<hx.GetNbinsX()+1; i++) {

    fInpHist.SetBinContent(i, hx.GetBinContent(i));
    fInpHist.SetBinError(i, hx.GetBinError(i));

  }

  fFuncInpH1.SetH1(&fInpHist);

  SetkDim(1);
  SetnCells(ncells);
  SetPseRan(&fRand);
  SetRho(&fFuncInpH1);
  Initialize();

}//FoamCellFinder

//_____________________________________________________________________________
void FoamCellFinder::GetCellPos1D(size_t i, Double_t& pos, Double_t& siz) {

  TFoamVect cell_pos(1), cell_siz(1);

  GetCell(i)->GetHcub(cell_pos, cell_siz);

  pos = cell_pos.GetCoord(0);
  siz = cell_siz.GetCoord(0);

  //cout << "GetCellPos1D: " << i << " " << cell_pos.GetCoord(0) << endl;
  TFoamCell *cell = GetCell(i);
  //cout << "GetCellPos1D: " << i << " " << cell->GetDau0() << " " << cell->GetDau1() << endl;

}//GetCellPos1D

//_____________________________________________________________________________
vector<Long_t> FoamCellFinder::GetActiveIdxSort() {

  map<Double_t, Long_t> pos_idx;

  //active cells loop
  for(const Long_t& i: fCellsAct) {

    Double_t pos=0, siz=0;
    GetCellPos1D(i, pos, siz);

    pos_idx.insert( make_pair(pos, i) );

  }//active cells loop

  vector<Long_t> idx;

  //map loop
  for(const pair<Double_t, Long_t>& i: pos_idx) {

    idx.push_back(i.second);

  }//map loop


  return idx;

}//GetActiveIdxSort

//_____________________________________________________________________________
TH1D FoamCellFinder::MakeH1D(const char *name_title) {

  //cout << "FoamCellFinder::MakeH1D" << endl;

  vector<Long_t> cell_idx = GetActiveIdxSort();

  vector<Double_t> bins;
  Double_t pos=0, siz=0;

  //cells loop
  for(const Long_t& i: cell_idx) {

    GetCellPos1D(i, pos, siz);

    if( fInpHist.GetNbinsX() > 1 ) {
      pos = pos*(fInpMax-fInpMin) + fInpMin;
    }

    bins.push_back(pos);

    //cout << bins.back() << endl;

  }//cells loop

  //end of last bin
  GetCellPos1D(cell_idx.back(), pos, siz);

  if( fInpHist.GetNbinsX() > 1 ) {
    pos = pos*(fInpMax-fInpMin) + fInpMin;
    siz *= fInpMax-fInpMin;
  }

  bins.push_back(pos+siz);

  //cout << bins.back() << endl;

  TH1D hx(name_title, name_title, bins.size()-1, bins.data());

  //cout << "FoamCellFinder::MakeH1D" << endl;

  return hx;

}//MakeH1D































