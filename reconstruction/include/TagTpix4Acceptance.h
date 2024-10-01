
#ifndef TagTpix4Acceptance_h
#define TagTpix4Acceptance_h

namespace TrkPlaneBasicHits {
  class Coll;
}

class TagTpix4Acceptance {

  public:

    void Run(const std::vector<std::string>& argvv);

  private:

  Double_t fTrueQ2;
  Double_t fTrueX;

    class plane_hits {

      public:

        plane_hits(std::string nam, TTree *in_tree);

        unsigned long GetNhits();

        std::unique_ptr<TrkPlaneBasicHits::Coll> hits;

    }; // plane hits

};

#endif

