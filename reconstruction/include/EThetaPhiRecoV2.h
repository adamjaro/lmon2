
#ifndef EThetaPhiRecoV2_h
#define EThetaPhiRecoV2_h

// Particle energy and angular reconstruction from measured quantities

class TH1D;

class EThetaPhiRecoV2 {

  public:

    EThetaPhiRecoV2(std::string nam, boost::program_options::options_description *opt=0x0);

    void MakeQuantity(std::string qnam);
    void Initialize(boost::program_options::variables_map *opt_map);
    void AddInput(Double_t *quant, Double_t en, Double_t theta, Double_t phi);
    void Finalize();
    void Export();

    void Import(TFile *in);
    void SetMinNinp(Int_t n) { fMinNinp = n; }
    void SetMaxRelEn(Double_t e) { fMaxRelEn = e; }
    void SetMaxRelPhi(Double_t e) { fMaxRelPhi = e; }
    void SetMaxErrTheta(Double_t e) { fMaxErrTheta = e; }
    Bool_t Reconstruct(Double_t *quant, Double_t& el_en, Double_t& el_theta, Double_t& el_phi, Int_t *ipar=0, Double_t *dpar=0);

  private:

    std::string fNam; // tool name
    boost::program_options::options_description *fOpt; // program options

    //cache for inputs
    TFile *fCacheFile=0x0; // cache file
    TTree *fCacheTree=0x0; // cache tree

    //inputs energy, theta and phi
    Double_t fCacheEn; // GeV
    Double_t fCacheTheta; // rad
    Double_t fCachePhi; //rad

    Int_t fMinNinp=0; // minimal number of inputs used for reconstruction

    //error limits on reconstructed kinematics, negative values for inactive criteria
    Double_t fMaxRelEn=-1; // maximal relative error in reconstructed energy
    Double_t fMaxRelPhi=-1; // maximal relative error in reconstructed phi
    Double_t fMaxErrTheta=-1; // maximal (absolute) error in reconstructed theta (rad)

    //measured quantity
    class Quantity {
      public:
      Quantity(std::string qnam): nam(qnam) {}

      std::string nam; // quantity name
      TH1D *hx=0x0; // quantity distribution from all inputs

      Double_t cache_val; // value for input cache

      //range in the quantity, uninitialized by max<min
      Double_t min=0; // minimal value
      Double_t max=0; // maximal value
      bool minmax_def=false; // range defined from configuration if true
      bool minmax_uninit=true; // range unitialized before first input call

      Quantity *quant_conf=0x0; // pointer to configuration quantity when in layer

    };//Quantity

    std::vector<Quantity> fQuantConf; // quantities configuration

    //link from detector quantities to particle energy and angles
    class Link {
      public:

      Link(): en(0), theta(0), phi(0), en_err(0), theta_err(0), phi_err(0) {}

      Double_t en=0; // particle energy, GeV
      Double_t theta=0; // particle polar angle, rad
      Double_t phi=0; // particle azimuthal angle, rad
      Double_t en_err=0; // corresponding errors, same units
      Double_t theta_err=0;
      Double_t phi_err=0;
      Int_t ninp=-1; // number of inputs set from import
      Int_t ilay=-1; // layer index for the link, set from import

      void AddParticle(double e, double t, double p);

      void Evaluate();

      unsigned int GetNinp() { return inp_en.size(); }

      private:

      Double_t GetMean(std::vector<double>& v);
      Double_t GetErr(std::vector<double>& v, Double_t m);

      std::vector<double> inp_en; // inputs in energy, GeV
      std::vector<double> inp_theta; // inputs in theta, rad
      std::vector<double> inp_phi; // inputs in phi, rad

    };//Link

    //Layer of particular segmentation holding links and quantities
    class Layer {
      public:

      Layer(Int_t idx, Int_t nbins=0): lay_idx(idx), lay_nbins(nbins) {}

      void AddInput(Double_t en, Double_t theta, Double_t phi); // add input to the layer

      ULong64_t GetLinkIdx(Double_t *quant); // link index

      Int_t lay_idx; // layer index, 0 is the top with the finest segmentation
      Int_t lay_nbins; // segmentation for the layer

      std::vector<Quantity> lay_quant; // quantities on the layer

      std::map<ULong64_t, Link> links; // links from detector quantities to particles

    };//Layer

    std::vector<Layer> fLayers;

};

#endif


















