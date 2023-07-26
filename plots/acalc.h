
// acceptance in a given kinematics variable

#if !defined(__CLING__) || defined(__ROOTCLING__)

//C++
#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <functional>

//ROOT
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TMath.h"

using namespace std;

#endif

//_____________________________________________________________________________
class acalc {
  public:

  ULong64_t nev; // number of events to process
  Double_t prec; // precision
  Double_t delt; // size of steps
  Double_t bmin; // minimal bin width
  Int_t modif; // modifier to input value

  //_____________________________________________________________________________
  acalc(TTree *t, string nam, string sel, string sel2="", int logic_sel=0): tree(t),
    prec(1e-2), delt(1e-6), nev(0), bmin(-1), modif(-1), val_sel2(0) {

    //input tree
    tree->SetBranchAddress(nam.c_str(), &val);

    //val_sel = new Bool_t;
    val_sel = new Int_t;
    tree->SetBranchAddress(sel.c_str(), val_sel);

    if( !sel2.empty() ) {
      //val_sel2 = new Bool_t;
      val_sel2 = new Int_t;
      tree->SetBranchAddress(sel2.c_str(), val_sel2);
    }

    //modifiers to input value
    mods[0] = &acalc::theta_to_eta;
    mods[1] = &acalc::log10_Q2;
    mods[2] = &acalc::pitheta;
    mods[3] = &acalc::mlt_theta;


    //logical function for two selection variables
    switch( logic_sel ) {

      case 0:
        select_func = logical_or<bool>();
        break;

      case 1:
        select_func = logical_and<bool>();
        break;

    }


  }//acalc

  //_____________________________________________________________________________
  TGraphAsymmErrors get() {

    //all and selected events
    list<Double_t> valAll, valSel;

    //tree loop
    if( nev==0 ) nev = tree->GetEntries();
    for(ULong64_t iev=0; iev<nev; iev++) {
      tree->GetEntry(iev);

      //modifiers to input value
      if( modif > -1 ) {
        val = (this->*mods[modif])(val);
      }

      //all events
      valAll.push_back(val);

      //make the selection
      bool select = false;

      //one or two selection variables
      if( val_sel2 ) {

        //two selection variables
        select = select_func(*val_sel, *val_sel2);
      } else {

        //single variable
        select = *val_sel;
      }

      if( !select ) continue;

      //event is selected
      valSel.push_back(val);

    }//tree loop

    //bins
    vector<Double_t> bins;
    get_bins(valAll, valSel, bins);

    //event distributions by get_acc.C
    TH1D hAll = h1_from_list("hAll", valAll, bins);
    TH1D hSel = h1_from_list("hSel", valSel, bins);

    cout << "nsel: " << hSel.GetEntries() << endl;
    cout << "nall: " << hAll.GetEntries() << endl;
    cout << "nsel/nall: " << hSel.GetEntries()/hAll.GetEntries() << endl;

    //acceptance
    TGraphAsymmErrors acc(&hSel, &hAll);

    cout << "Integral: " << integrate(acc) << endl;

    return acc;

  }//get

private:

  //_____________________________________________________________________________
  void get_bins(list<Double_t> &valAll, list<Double_t> &valSel, vector<Double_t> &bins) {

    //bin edges for the input lists

    //sort the lists
    valAll.sort();
    valSel.sort();

    //range in selected events
    Double_t minval = *(valSel.cbegin());
    Double_t maxval = *(valSel.crbegin());
    cout << "Range for bin edges: " << minval << " " << maxval << endl;

    //iterators for all and selected lists
    list<Double_t>::const_iterator itAll = valAll.cbegin();
    list<Double_t>::const_iterator itSel = valSel.cbegin();

    //move both lists to start of selected events
    while( *itAll < minval ) itAll++;
    while( *itSel < minval ) itSel++;

    //starting bin edge
    bins.push_back(minval);

    //counters for selected and all events
    UInt_t nsel = 0;
    UInt_t nall = 0;

    //current maximum
    Double_t cmax = minval + delt;

    //lists loop
    while(1) {
      //move to current maximum in step
      while( itSel != valSel.cend() and *itSel < cmax ) {
        nsel++;
        itSel++;
      }
      while( itAll != valAll.cend() and *itAll < cmax ) {
        nall++;
        itAll++;
      }

      //relative Binomial error in current step
      Double_t rel = 1.;
      if( nsel != 0 and nall != 0 and nsel != nall) {

        Double_t xsel = Double_t(nsel);
        Double_t xall = Double_t(nall);
        rel = TMath::Sqrt( (xall-xsel)/(xall*xsel) );
      }

      //test current step as bin edge candidate
      bool accept = false;

      //maximal allowed error
      if( rel < prec ) accept = true;

      //minimal bin width
      if( bmin > 0 and (cmax-bins[bins.size()-1]) < bmin ) {
        accept = false;
      }

      //result for current step as bin edge candidate
      if( accept ) {
        //mark edge from current step
        bins.push_back(cmax);

        //reset the counters for next bin
        nall = 0;
        nsel = 0;
        //cout << "bin edge: " << cmax << endl;
      }

      //increment current maximum for next step
      cmax += delt;

      //test for end of range in selected events
      if( cmax > maxval ) {
        //set last found edge to the end of values interval and finish
        bins.push_back(maxval);
        break;
      }

    }//lists loop

  }//get_bins

  //_____________________________________________________________________________
  void h1_from_list(list<Double_t> &inp, TH1D &hx) {

    //fill a TH1D histogram from the input list

    hx.Sumw2();

    //list loop
    list<Double_t>::const_iterator it = inp.cbegin();
    while(it != inp.cend()) {
      hx.Fill(*it);
      it++;
    }//list loop

  }//h1_from_list

  //_____________________________________________________________________________
  TH1D h1_from_list(string nam, list<Double_t>& inp, vector<Double_t>& bins) {

    //create TH1D from input list

    TH1D hx(nam.c_str(), nam.c_str(), bins.size()-1, bins.data());

    h1_from_list(inp, hx);

    return hx;

  }//h1_from_list

  //_____________________________________________________________________________
  Double_t integrate(const TGraphAsymmErrors& g) {

    //graph Riemann integral

    Double_t ig = 0;

    for(Int_t ip=0; ip<g.GetN(); ip++) {

      Double_t xp, yp;
      g.GetPoint(ip, xp, yp);

      ig += yp*(g.GetErrorXlow(ip) + g.GetErrorXhigh(ip));
    }

    return ig;

  }//integrate

  //_____________________________________________________________________________
  Double_t theta_to_eta(Double_t theta) {

    //pseudorapidity from polar angle

    return -TMath::Log(TMath::Tan(theta/2.));

  }//theta_to_eta

  //_____________________________________________________________________________
  Double_t log10_Q2(Double_t Q2) {

    //log_10(Q^2)

    return TMath::Log10(Q2);

  }//log10_Q2

  //_____________________________________________________________________________
  Double_t pitheta(Double_t theta) {

    //pi - theta  in mrad

    return (TMath::Pi()-theta)*1e3;

  }//pitheta

  //_____________________________________________________________________________
  Double_t mlt_theta(Double_t theta) {

    // -log_10(pi - theta)

    return -TMath::Log10(TMath::Pi()-theta);

  }//mlt_theta

  Double_t (acalc::*mods[4])(Double_t); // modifiers to input value

  TTree *tree; // input tree
  Double_t val; // kinematics variable
  //Bool_t *val_sel; // selection variable
  //Bool_t *val_sel2; // selection variable 2
  Int_t *val_sel; // selection variable
  Int_t *val_sel2; // selection variable 2

  function<bool(bool, bool)> select_func; // selection function

};//acalc
















