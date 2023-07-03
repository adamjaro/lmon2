
#ifndef TagTrackBasicVis_h
#define TagTrackBasicVis_h

#include <vector>

class TagTrackBasicVis {

  public:

    TagTrackBasicVis(const char *conf);

    std::string GetDetName() { return tag->GetName(); }
    void SetDet(int i);

    int ProcessEvent(bool *stat=0x0);
    int NextEvent(int di=1);
    int PreviousEvent();
    void SetEvent(int i) { iev = i; }

    void SetMaxChi2ndf(double chi2);
    void SetClsLimMdist(double d);

    void SetMinNtrk(int n) { min_ntrk = n; }
    void SetMinNcls(int n) { min_ncls = n; }
    void SetMinNcnt(int n) { min_ncnt = n; }
    void SetMinEtrk(int n) { min_etrk = n; }
    void SetMinSigTrk(int n) { min_sig_trk = n; }

    int GetNumberOfClusters();
    int GetNumberOfClusters(int iplane);
    void GetCluster(int iplane, int icls, double& x, double& y, double& z, double& md);

    int GetNumberOfTracks();
    int GetSigTracks();
    void GetTrack(int i, double& x0, double& y0, double& slope_x, double& slope_y, double& chi2, int& itrk);

    int GetNumberOfRefTracks();

    Double_t GetMaxChi2ndf() { return tag->GetMaxChi2ndf(); }
    Double_t GetClsLimMdist() { return tag->GetClsLimMdist(); }

  private:

    std::string GetStr(boost::program_options::variables_map& opt_map, std::string par);

    TChain *tree; // input tree
    TFile *out; // output file
    TTree *otree; // output tree

    TagTrackFindBasic *s1; // Tagger 1
    TagTrackFindBasic *s2; // Tagger 2

    //RefCounter *cnt_s1; // Reference counter 1
    //RefCounter *cnt_s2; // Reference counter 2

    Long64_t iev; // event number

    TagTrackFindBasic *tag; // active tagger
    //RefCounter *cnt; // active reference counter

    int min_ntrk; // minimal number of tracks in event
    int min_ncls; // minimal number of clusters on all planes in event
    int min_ncnt; // minimal number of reference tracks in event
    int min_etrk; // minimal number of excess tracks as reconstructed tracks - reference tracks
    int min_sig_trk; // minimal number of signal tracks (itrk == 1)

    //input MC particles
    //std::vector<Int_t> *fMCItrk=0x0; 
    //std::vector<Double_t> *fMCEn=0x0, *fMCTheta=0x0, *fMCPhi=0x0;

};



#endif

