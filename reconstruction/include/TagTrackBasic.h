
#ifndef TagTrackBasic_h
#define TagTrackBasic_h

// Reconstruction task to run TagTrackFindBasic tracker finder

#include <TMVA/MethodBase.h>
#include <TMVA/Reader.h>
#include <TH2F.h>

class EThetaPhiReco;
class EThetaPhiRecoV2;
class TagTrackFindBasic;

namespace MCParticles {
  class Coll;
}

enum LowQ2NNIndexIn{PosY,PosZ,DirX,DirY};
enum LowQ2NNIndexOut{MomX,MomY,MomZ};

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

    float nnInput[4];
    TMVA::Reader          m_reader{"!Color:!Silent"};
    TMVA::MethodBase*     m_method{nullptr};
    std::string m_method_name{"DNN_CPU"};
    std::string m_file_path{"/home/simon/geant4/lmon2/reconstruction/python/LowQ2_DNN_CPU.weights.xml"};
    bool useTMVA = true;
    float m_electron{0.000510998928}; //TODO: Link to constant elsewhere?

    TFile* oFile;
    TH2F*  fPosHist;
    TH2F*  fVecHist;

    MCParticles::Coll *fMC=0x0;

};

#endif

