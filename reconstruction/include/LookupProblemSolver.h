
#ifndef LookupProblemSolver_h
#define LookupProblemSolver_h

class TH1D;

class LookupProblemSolver {

  public:

    //constructor for both training and solving
    LookupProblemSolver(size_t nsol, std::string nam, boost::program_options::options_description *opt=nullptr);

    //training functions
    void MakeQuantity(std::string qnam); // add a new independent quantity
    void Initialize(boost::program_options::variables_map& opt_map); // initialize adding inputs
    void AddInput(const std::vector<double>& quant, const std::vector<double>& sol); // add input known solution
    void Finalize(); // process all added inputs before exporting
    void Export(); // export trained links to a file (opened outside before calling)

    //cell finders for processing the added inputs
    void SetUseFoamCellFinder(bool u=true) { if(u) { fCellFinder = CellFinder::kFoam; } }
    void SetUseMedCellFinder(bool u=true) { if(u) { fCellFinder = CellFinder::kMed; } }
    void SetUseSumCellFinder(bool u=true) { if(u) { fCellFinder = CellFinder::kSum; } }

    //solving functions
    void Import(TFile& in); // import trained links from input file

    //selection criteria for solutions
    void SetMinNinp(Int_t n) { fMinNinp = n; } // minimal number of inputs for a given solution
    void SetMaxErr(size_t isol, double v) { fMaxErr[isol].first = true; fMaxErr[isol].second = v; } // maximal allowed error
    void SetMaxRel(size_t isol, double v) { fMaxRel[isol].first = true; fMaxRel[isol].second = v; } // maximal relative error

    //get solution for a given quantities, optionally provide error on the solution and integer parameters
    bool Solve(std::vector<double>& quant, std::vector<double>& sol, std::vector<double> *err=nullptr, Int_t *ipar=nullptr);

  private:

    std::string fNam; // tool name

    //cache for inputs
    std::unique_ptr<TFile> fCacheFile;
    std::unique_ptr<TTree> fCacheTree;

    //input cache for solutions, holding also its dimensionality
    std::vector<double> fSolCache;

    //minimal number of inputs for desired solution
    Int_t fMinNinp=0;

    //criteria for absolute and relative errors on individual solutions
    std::vector<std::pair<bool, double>> fMaxErr;
    std::vector<std::pair<bool, double>> fMaxRel;

    //cell finder for variable binning in independent quantities
    enum class CellFinder { kNo=0, kFoam, kMed, kSum };
    CellFinder fCellFinder = CellFinder::kNo;

    //independent quantity
    class Quantity {
      public:
      Quantity(std::string qnam): nam(qnam) {}

      std::string nam; // quantity name
      std::string hx_nam; // name for quantity distribution hx
      std::shared_ptr<TH1D> hx_init; // initial distribution when hx is set for variable binning
      std::shared_ptr<TH1D> hx; // quantity distribution from all inputs

      Double_t cache_val; // value for input cache

      //range in the quantity, uninitialized by max<min
      Double_t min=0; // minimal value
      Double_t max=0; // maximal value
      bool minmax_uninit=true; // range unitialized before first input call

      Quantity *quant_conf = nullptr;

    };//Quantity

    //quantities configuration
    std::vector<Quantity> fQuantConf;

    //link from indepentent quantities to solutions
    class Link {
      public:

      Link(size_t nsol) {
        solution.assign(nsol, 0);
        solution_err.assign(nsol, 0);
      }

      //solution vector and its error
      std::vector<double> solution, solution_err;

      Int_t ninp=-1; // number of inputs set from import
      Int_t ilay=-1; // layer index for the link, set from import

      void AddInput(const std::vector<double>& s) { input_solutions.push_back(s); }

      void Evaluate();

      unsigned int GetNinp() { return input_solutions.size(); }

      private:

      Double_t GetMean(std::vector<double>& v);
      Double_t GetErr(std::vector<double>& v, Double_t m);

      //structure holding input solutions
      std::vector<std::vector<double>> input_solutions;

    };//Link

    //layer of particular segmentation holding links and quantities
    class Layer {
      public:

      Layer(Int_t idx, Int_t nbins=0): lay_idx(idx), lay_nbins(nbins) {}

      void AddInput(); // add input solution to the layer

      ULong64_t GetLinkIdx(const std::vector<Double_t>& quant); // link index

      Int_t lay_idx; // layer index, 0 is the top with the finest segmentation
      Int_t lay_nbins; // segmentation for the layer

      std::vector<Quantity> lay_quant; // quantities on the layer
      std::vector<double> *sol_cache = nullptr;

      std::map<ULong64_t, Link> links; // links from detector quantities to particles

    };//Layer

    std::vector<Layer> fLayers;


};

#endif



















