
#ifndef TagTrackBasicMakeResp_h
#define TagTrackBasicMakeResp_h

// Task to create reconstruction response by EThetaPhiReco

class EThetaPhiReco;
class EThetaPhiRecoV2;
class TagTrackFindBasic;

class TagTrackBasicMakeResp {

  public:

    TagTrackBasicMakeResp() {}

    void Run(const char *conf);

  private:

    //void AddInput(TagTrackFindBasic *tag, EThetaPhiReco *rec);
    void AddInput(TagTrackFindBasic *tag, EThetaPhiRecoV2 *rec);

    std::string GetStr(boost::program_options::variables_map& opt_map, std::string par);

    //input true kinematics
    Double_t fTrueEn=0;
    Double_t fTrueTheta=0;
    Double_t fTruePhi=0;

};

#endif

