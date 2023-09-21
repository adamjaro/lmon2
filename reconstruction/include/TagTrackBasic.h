
#ifndef TagTrackBasic_h
#define TagTrackBasic_h

// Reconstruction task to run TagTrackFindBasic tracker finder

class EThetaPhiReco;
class EThetaPhiRecoV2;
class TagTrackFindBasic;

namespace MCParticles {
  class Coll;
}

class TagTrackBasic {

  public:

    TagTrackBasic() {}

    void Run(const char *conf);

  private:

    //void ElectronRec(TagTrackFindBasic *tag, EThetaPhiReco *rec);
    void ElectronRec(TagTrackFindBasic *tag, EThetaPhiRecoV2 *rec);

    std::string GetStr(boost::program_options::variables_map& opt_map, std::string par);

    Double_t fBeamEn=0; // beam energy, GeV

    //input true kinematics
    Double_t fTrueEn=0;
    Double_t fTrueTheta=0;
    Double_t fTruePhi=0;
    Double_t fTrueQ2=0;
    Double_t fTrueX=0;

    MCParticles::Coll *fMC=0x0;

};

#endif

