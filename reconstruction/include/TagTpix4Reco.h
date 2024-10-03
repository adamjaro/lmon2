
#ifndef TagTpix4Reco_h
#define TagTpix4Reco_h

namespace TrkPlaneBasicHits {
  class Coll;
}

class TagTpix4Reco {

  public:

    void Run(const std::vector<std::string>& argvv);

  private:

  Double_t fTrueQ2;
  Double_t fTrueX;

};

#endif

