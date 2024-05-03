
#ifndef TagTrackRecLPS_h
#define TagTrackRecLPS_h

// Reconstruction task to run TagTrackFindBasic tracker finder

class LookupProblemSolver;
class TagTrackFindBasic;
class CalPWOClusterWavg;

namespace MCParticles {
  class Coll;
}

class TagTrackRecLPS {

  public:

    TagTrackRecLPS() {}

    void Run(const char *conf);

  private:

    void ElectronRec(TagTrackFindBasic& tag, std::unique_ptr<LookupProblemSolver>& rec);
    //void TrackCalMatch(TagTrackFindBasic *tag, CalPWOClusterWavg *cal);

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

