
#ifndef EThetaPhiReco_h
#define EThetaPhiReco_h

// Particle energy and angular reconstruction from measured quantities

class TH1D;

class EThetaPhiReco {

  public:

    EThetaPhiReco(std::string nam, boost::program_options::options_description *opt=0x0);

    void MakeQuantity(std::string qnam, double conv=1);
    void Initialize(boost::program_options::variables_map *opt_map);

    void AddInput(Double_t *quant, Double_t en, Double_t theta, Double_t phi);
    void Export();

    void Import(TFile *in);

    void SetMinNinp(Int_t n) { fMinNinp = n; }

    Bool_t Reconstruct(Double_t *quant, Double_t& el_en, Double_t& el_theta, Double_t& el_phi, Int_t *ninp=0);

  private:

    std::string fNam; // detector name
    boost::program_options::options_description *fOpt; // program options

    //measured quantity
    class Quantity {
      public:
      Quantity(std::string qnam, double c=1): nam(qnam), conv(c), hist(0x0) {}

      std::string nam; // quantity name
      double conv; // conversion units for the range
      TH1D *hist; // quantity range and distribution

    };//Quantity

    std::vector<Quantity> fQuant; // measured quantities for reconstruction

    //link from detector quantities to particle energy and angles
    class Link {
      public:

      Link(): en(0), theta(0), phi(0), en_err(0), theta_err(0), phi_err(0) {}

      Double_t en; // particle energy, GeV
      Double_t theta; // particle polar angle, rad
      Double_t phi; // particle azimuthal angle, rad
      Double_t en_err; // corresponding errors, same units
      Double_t theta_err;
      Double_t phi_err;
      Int_t ninp=0; // number of inputs set from import

      void AddParticle(double e, double t, double p);

      void Evaluate();

      unsigned int GetNinp() { return fEn.size(); }

      private:

      Double_t GetMean(std::vector<double>& v);
      Double_t GetErr(std::vector<double>& v, Double_t m);

      std::vector<double> fEn; // inputs in energy, GeV
      std::vector<double> fTheta; // inputs in theta, rad
      std::vector<double> fPhi; // inputs in phi, rad

    };//Link

    std::map<ULong64_t, Link> fLinks; // links from detector to particles

    ULong64_t GetIdx(Double_t *quant);

    template<typename T> void AddOpt(std::string onam);
    template<typename T> T GetOpt(std::string onam, boost::program_options::variables_map *opt_map);

    Int_t fMinNinp=0; // minimal number of inputs used for reconstruction

};

#endif

