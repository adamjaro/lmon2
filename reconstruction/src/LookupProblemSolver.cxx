
//_____________________________________________________________________________
//
// LookupProblemSolver, generalized from EThetaPhiRecoV2
//
//_____________________________________________________________________________

//C++
#include <vector>
#include <string>
#include <iostream>

//Boost
#include <boost/program_options.hpp>

//ROOT
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"

//local classes
#include "LookupProblemSolver.h"
#include "FoamCellFinder.h"
#include "MedCellFinder.h"
#include "SumCellFinder.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
LookupProblemSolver::LookupProblemSolver(size_t nsol, std::string nam, program_options::options_description& opt):
  fNsol(nsol), fNam(nam) {

  //program options for layers configuration with default values
  opt.add_options()
    ((fNam+".nlayers").c_str(), program_options::value<int>()->default_value(3), "Number of layers")
    ((fNam+".layer0_bins").c_str(), program_options::value<int>()->default_value(100), "First layer #0 segmentation")
    ((fNam+".layer_dec").c_str(), program_options::value<int>()->default_value(40), "Decrement in bins to each next layer")
    ((fNam+".layer0_bins_pow2").c_str(), program_options::value<int>()->default_value(0), "Segmentation in powers of 2")
  ;

  //initialize the solutions input
  fSolCache.assign(fNsol, 0);

}//LookupProblemSolver

//_____________________________________________________________________________
void LookupProblemSolver::MakeQuantity(std::string qnam) {

  fQuantConf.push_back( Quantity(qnam) );

}//MakeQuantity

//_____________________________________________________________________________
void LookupProblemSolver::Initialize(program_options::variables_map& opt_map) {

  cout << "LPS::Initialize, " << fNam << endl;

  //input cache
  fCacheFile = make_unique<TFile>( (fNam+"_tmp.root").c_str(), "recreate" );
  fCacheTree = make_unique<TTree>( (fNam+"_cache").c_str(), (fNam+"_cache").c_str() );

  //quantity loop for cache tree
  for(Quantity& i: fQuantConf) {

    fCacheTree->Branch(i.nam.c_str(), &i.cache_val, (i.nam+"/D").c_str());

  }//quantity loop for cache tree

  //solution inputs for cache tree
  for(size_t i=0; i<fSolCache.size(); i++) {

    //name for the solution branch
    stringstream snam;
    snam << "sol_" << i;

    fCacheTree->Branch(snam.str().c_str(), &fSolCache[i], (snam.str()+"/D").c_str());
  }

  //fCacheFile->ls();
  //fCacheTree->Print();

  //layers configuration from optional parameters
  Int_t nlayers = opt_map.at(fNam+".nlayers").as<int>(); // number of all layers
  Int_t layer0_bins = opt_map.at(fNam+".layer0_bins").as<int>(); // segmentation of the first layer #0
  Int_t layer_dec = opt_map.at(fNam+".layer_dec").as<int>(); // decrement in bins to each next layer

  //segmentation in powers of 2 when requested by non-zero value
  Int_t layer0_bins_pow2 = opt_map.at(fNam+".layer0_bins_pow2").as<int>();

  //create the layers
  for(Int_t i=0; i<nlayers; i++) {

    //bins for the layer
    Int_t nbins = layer0_bins-i*layer_dec;

    //bins in powers of 2 when requested
    if( layer0_bins_pow2 > 0 ) { nbins = TMath::Power(2, layer0_bins_pow2-i*layer_dec); }

    cout << fNam << ", constructing layer #" << i << ", bins: " << nbins << endl;
    fLayers.push_back( Layer(i, nbins) ); // layer number and segmentation
    //fLayers.back().sol_cache = make_shared<vector<double>>(fSolCache);
    fLayers.back().sol_cache = &fSolCache;
  }

  //layer loop to initialize the quantities
  for(Layer& i: fLayers) {

    //set quantities on the layer from quantity configuration
    i.lay_quant = fQuantConf;

    //quantity loop
    for(size_t iq=0; iq<fQuantConf.size(); iq++) {

      //set the corresponding configuration quantity
      i.lay_quant[iq].quant_conf = &fQuantConf[iq];
    }//quantity loop

  }//layer loop

}//Initialize

//_____________________________________________________________________________
void LookupProblemSolver::AddInput(const vector<double>& quant, const vector<double>& sol) {

  //cout << "LPS: AddInput: ";
  //cout << quant[0] << " " << quant[1] << " " << quant[2] << " " << quant[3] << " ";
  //cout << sol[0] << " " << sol[1] << " " << sol[2] << " " << endl;

  //quantity loop for quantity input
  for(size_t iquant=0; iquant<fQuantConf.size(); iquant++) {
    Quantity& i = fQuantConf[iquant];

    //set cache value from the input
    i.cache_val = quant[iquant];

    //range in the quantity
    if( i.minmax_uninit ) {

      //starting values at the first call
      i.min = i.cache_val;
      i.max = i.cache_val;
      i.minmax_uninit = false;
    } else {
      //move the range values
      if( i.cache_val < i.min ) i.min = i.cache_val;
      if( i.cache_val > i.max ) i.max = i.cache_val;
    }

    //cout << i.min << " " << i.max << endl;
    //cout << fQuantConf[iquant].min << " " << fQuantConf[iquant].max << endl;
  }//quantity loop

  //input solution
  for(size_t isol=0; isol<fSolCache.size(); isol++) {

    fSolCache[isol] = sol[isol];
  }//input solution

  //fill the cache tree
  fCacheTree->Fill();

}//AddInput

//_____________________________________________________________________________
void LookupProblemSolver::Finalize() {

  //add inputs to layers from cached values

  //layer loop to set quantity distributions
  for(Layer& i: fLayers) {

    //quantity loop
    for(Quantity& iq: i.lay_quant) {

      //binning for the layer and range for the quantity
      iq.hx_nam = fNam+"_hx_"+to_string(i.lay_idx)+"_"+iq.nam;
      Double_t hmin = iq.quant_conf->min;
      Double_t hmax = iq.quant_conf->max;
      //quantity histogram, extend by one more bin
      Double_t binsiz = (hmax-hmin)/i.lay_nbins;

      //cout << iq.hx_nam << " " << i.lay_nbins << " " << hmin << " " << hmax << endl;

      if( fCellFinder != CellFinder::kNo ) {
        iq.hx_init = make_shared<TH1D>((iq.hx_nam+"_init").c_str(), (iq.hx_nam+"_init").c_str(), i.lay_nbins, hmin-binsiz, hmax+binsiz);
      } else {
        iq.hx = make_shared<TH1D>(iq.hx_nam.c_str(), iq.hx_nam.c_str(), i.lay_nbins, hmin-binsiz, hmax+binsiz);
      }
    }//quantity loop
  }//layer loop

  //variable binning for quantity distributions when requested
  if( fCellFinder != CellFinder::kNo ) {

    //cache loop for initial distribution
    Long64_t ncache = fCacheTree->GetEntries();
    for(Long64_t ic=0; ic<ncache; ic++) {
      fCacheTree->GetEntry(ic);

      for(Layer& i: fLayers) {
        for(Quantity& iq: i.lay_quant) {
          iq.hx_init->Fill( iq.quant_conf->cache_val );
        }
      }
    }//cache loop for initial distribution

    //make the final quantity distributions
    for(Layer& i: fLayers) {
      for(Quantity& iq: i.lay_quant) {

        //select the specific finder
        switch(fCellFinder) {
          case CellFinder::kFoam: {
            //variable binning based on TFoam
            FoamCellFinder finder(*iq.hx_init, 2*iq.hx_init->GetNbinsX());
            iq.hx = make_shared<TH1D>( finder.MakeH1D(iq.hx_nam.c_str()) );
          } break;
          case CellFinder::kMed: {
            //variable binning based on median cell split
            MedCellFinder finder(iq.hx_init.get(), 2*iq.hx_init->GetNbinsX());
            iq.hx = make_shared<TH1D>( finder.MakeH1D(iq.hx_nam.c_str()) );
          } break;
          case CellFinder::kSum: {
            //variable binning based on even sum per bin
            iq.hx = make_shared<TH1D>( SumCellFinder(iq.hx_init.get(), iq.hx_nam.c_str()) );
          } break;
        }//switch
      }//quantity loop
    }//layer loop

  }//variable binning condition

  //cache loop
  Long64_t ncache = fCacheTree->GetEntries();
  for(Long64_t ic=0; ic<ncache; ic++) {

    fCacheTree->GetEntry(ic);

    //input from the cache to all layers
    for(Layer& ilay: fLayers) {
      ilay.AddInput();
    }

  }//cache loop

}//Finalize

//_____________________________________________________________________________
void LookupProblemSolver::Export() {

  //layer loop
  for(Layer& i: fLayers) {

    //quantities for the layer
    TList qlist, qlist_init;

    //quantity loop
    for(Quantity& iq: i.lay_quant) {

      qlist.AddLast(iq.hx.get());
      if(iq.hx_init) { qlist_init.AddLast(iq.hx_init.get()); }
    }//quantity loop

    //export the quantities
    string qnam = fNam+"_quantities_"+to_string(i.lay_idx);
    qlist.Write(qnam.c_str(), TObject::kSingleKey);

    //initial quantity distribution
    if( qlist_init.GetEntries() > 0 ) {

      qnam = fNam+"_quantities_init_"+to_string(i.lay_idx);
      qlist_init.Write(qnam.c_str(), TObject::kSingleKey);
    }

    //tree on links for the layer
    string tnam = fNam+"_links_"+to_string(i.lay_idx);
    TTree link_tree(tnam.c_str(), tnam.c_str());
    ULong64_t idx;
    Int_t ninp;
    vector<double> sol_tree, sol_err_tree;

    //index branches
    link_tree.Branch("idx", &idx, "idx/l");
    link_tree.Branch("ninp", &ninp, "ninp/I");

    //solution branches
    sol_tree.resize(fSolCache.size());
    sol_err_tree.resize(fSolCache.size());
    for(size_t i=0; i<fSolCache.size(); i++) {
      string snam = "sol_"+to_string(i);
      string snam_err = "sol_err_"+to_string(i);
 
      link_tree.Branch(snam.c_str(), &sol_tree[i], (snam+"/D").c_str());
      link_tree.Branch(snam_err.c_str(), &sol_err_tree[i], (snam_err+"/D").c_str());
    }

    //link loop for the layer
    for(map<ULong64_t, Link>::iterator ilnk = i.links.begin(); ilnk != i.links.end(); ilnk++) {

      //link index
      idx = (*ilnk).first;

      //get the link and evaluate its input particles
      Link& lnk = (*ilnk).second;
      lnk.Evaluate();

      //link parameters to the tree
      ninp = lnk.GetNinp();
      sol_tree = lnk.solution;
      sol_err_tree = lnk.solution_err;

      link_tree.Fill();

    }//link loop

    link_tree.Write(0, TObject::kOverwrite);

  }//layer loop

  //close and remove the cache file
  fCacheFile->Close();

  string rm_command("rm -f ");
  rm_command += fCacheFile->GetName();
  int stat = system( rm_command.c_str() );
  if( stat != 0 ) {
    cout << fNam << ", can't remove cache file: " << fCacheFile->GetName() << endl;
    cout << "Status code: " << stat << endl;
  }

}//Export

//_____________________________________________________________________________
void LookupProblemSolver::Layer::AddInput() {

  //add input for the layer

  //cout << "LPS::Layer::AddInput" << endl;

  //quantity values for the input
  vector<Double_t> quant;
  quant.reserve( lay_quant.size() );

  //quantity loop
  for(size_t i=0; i<lay_quant.size(); i++) {

    //current cached value
    Double_t v = lay_quant[i].quant_conf->cache_val;

    //set the value for distribution and for link index
    lay_quant[i].hx->Fill( v );
    quant[i] = v;

  }//quantity loop

  //index from quantity values for the link on the layer
  ULong64_t idx = GetLinkIdx( quant );

  //test for the range
  if( idx == 0 ) return;

  //link to the solution
  map<ULong64_t, Link>::iterator ilnk = links.find(idx);

  //add the link if not present
  if( ilnk == links.end() ) {
    ilnk = links.insert( make_pair(idx, Link(sol_cache->size())) ).first;
  }

  //add solution input to the link
  Link& lnk = (*ilnk).second;
  lnk.AddInput(*sol_cache);

}//Layer::AddInput

//_____________________________________________________________________________
ULong64_t LookupProblemSolver::Layer::GetLinkIdx(const vector<Double_t>& quant) {

  //linear index for link corresponding to provided quantities

  ULong64_t idx = 0;

  //quantity loop
  for(unsigned int iq=0; iq<lay_quant.size(); iq++) {

    //offset by bin index
    Int_t ibin = lay_quant[iq].hx->FindBin( quant[iq] );

    //test for the range
    if( ibin < 1 or ibin > lay_quant[iq].hx->GetNbinsX() ) return 0;

    //segment base by quantity position
    ULong64_t seg = 1;
    for(unsigned int is=0; is<iq; is++) {
      seg *= ULong64_t( lay_quant[is].hx->GetNbinsX() );
    }

    idx += ULong64_t(ibin)*seg;

  }//quantity loop

  return idx;

}//Layer::GetLinkIdx

//_____________________________________________________________________________
void LookupProblemSolver::Link::Evaluate() {

  //cout << "LPS::Link::Evaluate" << endl;

  //assign link solution from the inputs provided to the link

  //more inputs for the link
  if( input_solutions.size() > 1 ) {

    //loop over individual solutions
    for(size_t isol=0; isol<solution.size(); isol++) {

      //inputs for the given solution
      vector<double> inp;

      //inputs loop
      for(const vector<double>& i:input_solutions) {

        inp.push_back( i[isol] );
      }//inputs loop

      //set the solution and its error from the inputs
      solution[isol] = GetMean(inp);
      solution_err[isol] = GetErr(inp, solution[isol]);

    }//loop over individual solutions

  } else {

    //just one input solution in the link

    //solutions loop
    for(size_t isol=0; isol<solution.size(); isol++) {

      solution[isol] = input_solutions[0][isol];

    }//solutions loop

  }


}//Link::Evaluate

//_____________________________________________________________________________
Double_t LookupProblemSolver::Link::GetMean(std::vector<double>& v) {

  //mean for a set of values

  Double_t nx = Double_t(v.size());

  Double_t sum = 0.;
  for(unsigned int i=0; i<v.size(); i++) {
    sum += v[i];
  }

  return sum/nx;

}//GetMean

//_____________________________________________________________________________
Double_t LookupProblemSolver::Link::GetErr(std::vector<double>& v, Double_t m) {

  //error for the mean

  if( v.size() <= 1 ) return -1;

  Double_t nx = Double_t(v.size());

  Double_t sum = 0.;
  for(unsigned int i=0; i<v.size(); i++) {
    sum += (v[i]-m)*(v[i]-m);
  }

  return TMath::Sqrt( sum/(nx*(nx-1)) );

}//GetMean



























