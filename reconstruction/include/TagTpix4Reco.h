
#ifndef TagTpix4Reco_h
#define TagTpix4Reco_h

// Reconstruction for tagger detectors with Timepix4

namespace MCParticles { class Coll;}
namespace TagTracks { class Track;}

class TagTpix4Reco {

  public:

    void Run(const std::vector<std::string>& argvv);

  private:

    void ElectronRec(std::vector<TagTracks::Track>& trk);

    //input true kinematics
    Double_t fTrueEn=0;
    Double_t fTrueTheta=0;
    Double_t fTruePhi=0;
    Double_t fTrueQ2=0;
    Double_t fTrueX=0;
    Double_t fNumInteractions=0;

    std::shared_ptr<MCParticles::Coll> fMC;

};

#endif

