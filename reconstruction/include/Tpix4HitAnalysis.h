
#ifndef Tpix4HitAnalysis_h
#define Tpix4HitAnalysis_h

//hit analysis for Timepix4

class GeoParser;
namespace TrkPlaneBasicHits {
  class Coll;
}

#include "TH2D.h"

class Tpix4HitAnalysis {

  public:

    Tpix4HitAnalysis() {}

    void Run(const char *conf);

  private:

    std::string GetStr(boost::program_options::variables_map& opt_map, std::string par);
    void ShowProgress(Double_t xi, Double_t xall);

    class plane_hits {

      public: 

        plane_hits(std::string nam, GeoParser *geo, TTree *in_tree);

        void run_event();

        Int_t ipix, irow;
        Double_t en;

        std::unique_ptr<TrkPlaneBasicHits::Coll> hits;

        TH2D h_counts;
        TTree h_tree;

    }; // plane hits

};

#endif










