
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

    Double_t fTrueQ2;
    Double_t fTrueX;
    Double_t fNumInteractions;

    std::shared_ptr<MCParticles::Coll> fMC;

};

#endif

