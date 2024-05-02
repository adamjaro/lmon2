
#ifndef TagTrackMakeLPS_h
#define TagTrackMakeLPS_h

class TagTrackFindBasic;
class LookupProblemSolver;

class TagTrackMakeLPS {

  public:

    TagTrackMakeLPS() {}

    void Run(const char *conf);

  private:

    void AddInput(TagTrackFindBasic *tag, LookupProblemSolver *rec);

    std::string GetStr(boost::program_options::variables_map& opt_map, std::string par);

    //input true kinematics
    Double_t fTrueEn=0;
    Double_t fTrueTheta=0;
    Double_t fTruePhi=0;

};

#endif

