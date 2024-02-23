
//_____________________________________________________________________________
//
// Particle energy and polar theta and azimuthal phi angles
// are reconstructed out of a set of measured quantities,
// version 2
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
#include "TMath.h"
#include "TUnixSystem.h"

//local classes
#include "EThetaPhiRecoV2.h"
#include "FoamCellFinder.h"
#include "MedCellFinder.h"
#include "SumCellFinder.h"

using namespace std;
using namespace boost;

//_____________________________________________________________________________
EThetaPhiRecoV2::EThetaPhiRecoV2(string nam, program_options::options_description *opt):
    fNam(nam), fOpt(opt) {

  if(!opt) return;

  //program options for layers configuration with default values
  fOpt->add_options()
    ((fNam+".nlayers").c_str(), program_options::value<int>()->default_value(10), "Number of layers")
    ((fNam+".layer0_bins").c_str(), program_options::value<int>()->default_value(1000), "First layer #0 segmentation")
    ((fNam+".layer_dec").c_str(), program_options::value<int>()->default_value(100), "Decrement in bins to each next layer")
    ((fNam+".layer0_bins_pow2").c_str(), program_options::value<int>()->default_value(0), "Segmentation in powers of 2")
  ;

}//EThetaPhiRecoV2

//_____________________________________________________________________________
void EThetaPhiRecoV2::MakeQuantity(std::string qnam) {

  fQuantConf.push_back( Quantity(qnam) );

}//MakeQuantity

//_____________________________________________________________________________
void EThetaPhiRecoV2::Initialize(program_options::variables_map *opt_map) {

  cout << "EThetaPhiRecoV2::Initialize, " << fNam << endl;

  //input cache
  fCacheFile = new TFile((fNam+"_tmp.root").c_str(), "recreate");
  fCacheTree = new TTree((fNam+"_cache").c_str(), (fNam+"_cache").c_str());

  //quantity loop for cache tree
  for(Quantity& i: fQuantConf) {

    fCacheTree->Branch(i.nam.c_str(), &i.cache_val, (i.nam+"/D").c_str());

  }//quantity loop for cache tree

  //kinematics input for cache tree
  fCacheTree->Branch("inp_en", &fCacheEn, "inp_en/D");
  fCacheTree->Branch("inp_theta", &fCacheTheta, "inp_theta/D");
  fCacheTree->Branch("inp_phi", &fCachePhi, "inp_phi/D");

  //fCacheFile->ls();
  //fCacheTree->Print();

  //layers configuration from optional parameters
  Int_t nlayers = opt_map->at(fNam+".nlayers").as<int>(); // number of all layers
  Int_t layer0_bins = opt_map->at(fNam+".layer0_bins").as<int>(); // segmentation of the first layer #0
  Int_t layer_dec = opt_map->at(fNam+".layer_dec").as<int>(); // decrement in bins to each next layer

  //segmentation in powers of 2 when requested by non-zero value
  Int_t layer0_bins_pow2 = opt_map->at(fNam+".layer0_bins_pow2").as<int>();

  //create the layers
  for(Int_t i=0; i<nlayers; i++) {

    //bins for the layer
    Int_t nbins = layer0_bins-i*layer_dec;

    //bins in powers of 2 when requested
    if( layer0_bins_pow2 > 0 ) { nbins = TMath::Power(2, layer0_bins_pow2-i*layer_dec); }

    cout << fNam << ", constructing layer #" << i << ", bins: " << nbins << endl;
    fLayers.push_back( Layer(i, nbins) ); // layer number and segmentation
  }

  //layer loop to initialize the quantities
  for(Layer& i: fLayers) {
    i.lay_quant = fQuantConf;

    //quantity loop
    for(size_t iq=0; iq<fQuantConf.size(); iq++) {

      //set the pointer to corresponding configuration quantity
      i.lay_quant[iq].quant_conf = &fQuantConf[iq];
    }//quantity loop
  }//layer loop

}//Initialize

//_____________________________________________________________________________
void EThetaPhiRecoV2::AddInput(Double_t *quant, Double_t en, Double_t theta, Double_t phi) {

  //loop over input quantites, same order as they were added
  Int_t iquant=0;
  for(Quantity& i: fQuantConf) {

    i.cache_val = quant[iquant];
    iquant++;

    //range in the quantity
    if( i.minmax_def ) continue;

    //starting values at the first call
    if( i.minmax_uninit ) {

      i.min = i.cache_val;
      i.max = i.cache_val;
      i.minmax_uninit = false;
    } else {
      //move the range values
      if( i.cache_val < i.min ) i.min = i.cache_val;
      if( i.cache_val > i.max ) i.max = i.cache_val;
    }

  }//loop over input quantites

  //input kinematics
  fCacheEn = en;
  fCacheTheta = theta;
  fCachePhi = phi;

  //fill the cache tree
  fCacheTree->Fill();

}//AddInput

//_____________________________________________________________________________
void EThetaPhiRecoV2::Finalize() {

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
      if( fCellFinder != CellFinder::kNo ) {
        iq.hx_init = new TH1D((iq.hx_nam+"_init").c_str(), (iq.hx_nam+"_init").c_str(), i.lay_nbins, hmin-binsiz, hmax+binsiz);
      } else {
        iq.hx = new TH1D(iq.hx_nam.c_str(), iq.hx_nam.c_str(), i.lay_nbins, hmin-binsiz, hmax+binsiz);
      }

    }//quantity loop
  }//layer loop

  //variable binning for quantity distributions when requested
  if( fCellFinder != CellFinder::kNo ) {

    //turn off kinematics to read only the quantities
    fCacheTree->SetBranchStatus("inp_en", false);
    fCacheTree->SetBranchStatus("inp_theta", false);
    fCacheTree->SetBranchStatus("inp_phi", false);

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

    //enable the kinematics
    fCacheTree->SetBranchStatus("inp_en", true);
    fCacheTree->SetBranchStatus("inp_theta", true);
    fCacheTree->SetBranchStatus("inp_phi", true);

    //make the final quantity distributions
    for(Layer& i: fLayers) {
      for(Quantity& iq: i.lay_quant) {

        //select the specific finder
        switch(fCellFinder) {
          case CellFinder::kFoam: {
            //variable binning based on TFoam
            FoamCellFinder finder(*iq.hx_init, 2*iq.hx_init->GetNbinsX());
            iq.hx = new TH1D( finder.MakeH1D(iq.hx_nam.c_str()) );
          } break;
          case CellFinder::kMed: {
            //variable binning based on median cell split
            MedCellFinder finder(iq.hx_init, 2*iq.hx_init->GetNbinsX());
            iq.hx = new TH1D( finder.MakeH1D(iq.hx_nam.c_str()) );
          } break;
          case CellFinder::kSum: {
            //variable binning based on even sum per bin
            iq.hx = new SumCellFinder(iq.hx_init, iq.hx_nam.c_str());
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
      ilay.AddInput(fCacheEn, fCacheTheta, fCachePhi);
    }

  }//cache loop

}//Finalize

//_____________________________________________________________________________
void EThetaPhiRecoV2::Export() {

  //layer loop
  for(Layer& i: fLayers) {

    //quantities for the layer
    TList qlist, qlist_init;

    //quantity loop
    for(Quantity& iq: i.lay_quant) {

      qlist.AddLast(iq.hx);
      if(iq.hx_init) { qlist_init.AddLast(iq.hx_init); }
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
    Double_t en, theta, phi, en_err, theta_err, phi_err;
    Int_t ninp;
    link_tree.Branch("idx", &idx, "idx/l");
    link_tree.Branch("en", &en, "en/D");
    link_tree.Branch("theta", &theta, "theta/D");
    link_tree.Branch("phi", &phi, "phi/D");
    link_tree.Branch("ninp", &ninp, "ninp/I");
    link_tree.Branch("en_err", &en_err, "en_err/D");
    link_tree.Branch("theta_err", &theta_err, "theta_err/D");
    link_tree.Branch("phi_err", &phi_err, "phi_err/D");

    //link loop for the layer
    for(map<ULong64_t, Link>::iterator ilnk = i.links.begin(); ilnk != i.links.end(); ilnk++) {

      //link index
      idx = (*ilnk).first;

      //get the link and evaluate its input particles
      Link& lnk = (*ilnk).second;
      lnk.Evaluate();

      //link parameters to the tree
      en = lnk.en;
      theta = lnk.theta;
      phi = lnk.phi;
      ninp = lnk.GetNinp();
      en_err = lnk.en_err;
      theta_err = lnk.theta_err;
      phi_err = lnk.phi_err;

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
void EThetaPhiRecoV2::Import(TFile *in) {

  //import from file

  cout << "EThetaPhiRecoV2::Import, " << fNam << endl;

  //count how many layers are there
  Int_t nlay=0;
  while( in->Get( (fNam+"_quantities_"+to_string(nlay)).c_str() ) ) {
    nlay++;
  }

  cout << fNam << ", layers to import: " << nlay << endl;

  //layers loop
  for(Int_t ilay=0; ilay<nlay; ilay++) {

    //make the layer
    fLayers.push_back( Layer(ilay) );
    Layer& lay = fLayers.back();

    //set the quantities for the layer
    TList *qlist = dynamic_cast<TList*>(in->Get( (fNam+"_quantities_"+to_string(ilay)).c_str() ));

    //quantity loop to attach them to the layer
    for(Int_t iq=0; iq<qlist->GetEntries(); iq++) {

      lay.lay_quant.push_back( Quantity(qlist->At(iq)->GetName()) );
      lay.lay_quant.back().hx = dynamic_cast<TH1D*>( qlist->At(iq) );
    }//quantity loop

    //import the links for the layer
    TTree *link_tree = dynamic_cast<TTree*>(in->Get( (fNam+"_links_"+to_string(ilay)).c_str() ));
    ULong64_t idx;
    Double_t en, theta, phi, en_err, theta_err, phi_err;
    Int_t ninp;
    link_tree->SetBranchAddress("idx", &idx);
    link_tree->SetBranchAddress("en", &en);
    link_tree->SetBranchAddress("theta", &theta);
    link_tree->SetBranchAddress("phi", &phi);
    link_tree->SetBranchAddress("ninp", &ninp);
    link_tree->SetBranchAddress("en_err", &en_err);
    link_tree->SetBranchAddress("theta_err", &theta_err);
    link_tree->SetBranchAddress("phi_err", &phi_err);

    //links tree loop
    for(Long64_t i=0; i<link_tree->GetEntries(); i++) {
      link_tree->GetEntry(i);

      //make the link for the layer from the tree entry
      map<ULong64_t, Link>::iterator ilnk = lay.links.insert( make_pair(idx, Link()) ).first;

      //set the link
      Link& lnk = (*ilnk).second;
      lnk.en = en;
      lnk.theta = theta;
      lnk.phi = phi;
      lnk.ninp = ninp;
      lnk.ilay = ilay;
      lnk.en_err = en_err;
      lnk.theta_err = theta_err;
      lnk.phi_err = phi_err;

    }//links tree loop

  }//layers loop

}//Import

//_____________________________________________________________________________
Bool_t EThetaPhiRecoV2::Reconstruct(Double_t *quant,
    Double_t& el_en, Double_t& el_theta, Double_t& el_phi, Int_t *ipar, Double_t *dpar) {

  //run reconstruction for the measured quantities

  //layer loop to find the first one meeting selection criteria
  for(Layer& lay: fLayers) {

    //link on the given layer
    map<ULong64_t, Link>::iterator ilnk = lay.links.find( lay.GetLinkIdx(quant) );

    //link is present
    if( ilnk == lay.links.end() ) continue;

    //get the link
    Link& lnk = (*ilnk).second;

    //selection criteria

    //minimal number of inputs
    if( lnk.ninp < fMinNinp ) continue;

    //maximal relative error in energy
    if( fMaxRelEn > 0 and lnk.en != 0 and ( lnk.en_err < 0 or lnk.en_err/lnk.en > fMaxRelEn ) ) continue;

    //maximal relative error in phi
    if( fMaxRelPhi > 0 and lnk.phi != 0 and ( lnk.phi_err < 0 or lnk.phi_err/lnk.phi > fMaxRelPhi ) ) continue;

    //maximal error in theta
    if( fMaxErrTheta > 0 and ( lnk.theta_err < 0 or lnk.theta_err > fMaxErrTheta ) ) continue;

    //return the result of reconstruction
    el_en = lnk.en;
    el_theta = lnk.theta;
    el_phi = lnk.phi;
    //if(ninp) *ninp = lnk.ninp;
    if(ipar) {
      ipar[0] = lnk.ninp;
      ipar[1] = lnk.ilay;
    }
    if(dpar) {
      dpar[0] = lnk.en_err;
      dpar[1] = lnk.theta_err;
      dpar[2] = lnk.phi_err;
    }

    return kTRUE;

  }//layer loop

  //no link satisfying the criteria found in any layer
  return kFALSE;

}//Reconstruct

//_____________________________________________________________________________
void EThetaPhiRecoV2::Layer::AddInput(Double_t en, Double_t theta, Double_t phi) {

  //quantity values
  vector<Double_t> quant;
  quant.reserve( lay_quant.size() );

  //quantity loop
  for(size_t i=0; i<lay_quant.size(); i++) {

    //current cached value
    Double_t v = lay_quant[i].quant_conf->cache_val;

    //set the value for distribution and link index
    lay_quant[i].hx->Fill( v );
    quant[i] = v;

  }//quantity loop

  //index for the link on the layer
  ULong64_t idx = GetLinkIdx( quant.data() );

  //test for the range
  if( idx == 0 ) return;

  //link to the particle
  map<ULong64_t, Link>::iterator ilnk = links.find(idx);

  //add the link if not present
  if( ilnk == links.end() ) {
    ilnk = links.insert( make_pair(idx, Link()) ).first;
  }

  //add particle input to the link
  Link& lnk = (*ilnk).second;
  lnk.AddParticle(en, theta, phi);

}//Layer::AddInput

//_____________________________________________________________________________
ULong64_t EThetaPhiRecoV2::Layer::GetLinkIdx(Double_t *quant) {

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
void EThetaPhiRecoV2::Link::AddParticle(double e, double t, double p) {

  //add particle for the link

  inp_en.push_back(e);
  inp_theta.push_back(t);
  inp_phi.push_back(p);

}//Link::AddParticle

//_____________________________________________________________________________
void EThetaPhiRecoV2::Link::Evaluate() {

  //assign particle energy and angles from inputs for a given link

  //more inputs for the link
  if( inp_en.size() > 1 ) {

    en = GetMean(inp_en);
    theta = GetMean(inp_theta);
    phi = GetMean(inp_phi);
    en_err = GetErr(inp_en, en);
    theta_err = GetErr(inp_theta, theta);
    phi_err = GetErr(inp_phi, phi);

  } else {

    //just one particle input in the link
    en = inp_en[0];
    theta = inp_theta[0];
    phi = inp_phi[0];
    en_err = -1;
    theta_err = -1;
    phi_err = -1;
  }

}//Link::Evaluate

//_____________________________________________________________________________
Double_t EThetaPhiRecoV2::Link::GetMean(std::vector<double>& v) {

  //mean for a set of values

  Double_t nx = Double_t(v.size());

  Double_t sum = 0.;
  for(unsigned int i=0; i<v.size(); i++) {
    sum += v[i];
  }

  return sum/nx;

}//GetMean

//_____________________________________________________________________________
Double_t EThetaPhiRecoV2::Link::GetErr(std::vector<double>& v, Double_t m) {

  //error for the mean

  if( v.size() <= 1 ) return -1;

  Double_t nx = Double_t(v.size());

  Double_t sum = 0.;
  for(unsigned int i=0; i<v.size(); i++) {
    sum += (v[i]-m)*(v[i]-m);
  }

  return TMath::Sqrt( sum/(nx*(nx-1)) );

}//GetMean

















