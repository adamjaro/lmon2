
#ifndef TagTrackFinder_h
#define TagTrackFinder_h

// Track finder for timepix tracker

#include "TagTracks.h"

class GeoParser;

class TagTrackFinder {

  public:

    TagTrackFinder(std::string nam, std::vector<std::shared_ptr<TrkClusterFinder>>& cls, size_t iofs=0);

    void LoadV6p3_rev3(GeoParser& geo, std::string p1, std::string p2, std::string par);

    void CreateOutput(TTree *tree) { fTracks.CreateOutput(fNam, tree); }

    void SetMaxChi2Ndf(Double_t max_chi2ndf) { fChi2ndfMax = max_chi2ndf; }
    Double_t GetMaxChi2ndf() { return fChi2ndfMax; }

    void ProcessEvent();
    void FinishEvent() { fTracks.FinishEvent(); }
    void PrintCounts();

    std::vector<TagTracks::Track>& GetTracks() { return fTracks.GetStore(); }

  private:

    void MakeTrack(Double_t *x, Double_t& pos, Double_t& slope, Double_t& theta, Double_t& chi2);
    Double_t TrackChi2(Double_t *x, Double_t *y, TagTracks::Track& trk);

    template<std::size_t N> void ClusterAnalysis(TrkPlaneClusters::Cluster* (&cls)[N], TagTracks::Track& trk);

    std::string fNam; // finder name, used for output

    std::vector<std::shared_ptr<TrkClusterFinder>> fPlanes; // planes for the station

    Double_t fL=0; // plane spacing, mm
    Double_t fZ[4] = {0,0,0,0}; // local z positions for planes, mm

    TagTracks::Coll fTracks; // tracks collection

    Double_t fChi2ndfMax=0; // maximal reduced chi2 for tracks

    //counter for tracks
    unsigned long fNall=0;

};

#endif

