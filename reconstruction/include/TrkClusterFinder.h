
#ifndef TrkClusterFinder_h
#define TrkClusterFinder_h

#include "TrkPlaneClusters.h"
#include "TrkPlaneBasicHits.h"

class TrkClusterFinder {

  public:

    TrkClusterFinder(std::string nam, TrkPlaneBasicHits::Coll& hits);

    void CreateOutput(TTree *tree);

    void ProcessEvent();
    void FinishEvent() { fClusters.FinishEvent(); }
    void PrintCounts();

    void SetLimMdist(Double_t d) { fClsMinLimMdist = d; }
    Double_t GetLimMdist() { return fClsMinLimMdist; }

    const std::string& GetName() { return fNam; }

    TrkPlaneClusters::Coll& GetClusters() { return fClusters; }

  private:

    int GetHitsCount();
    unsigned long FindHitEmax();
    int FindAdjHits(unsigned long ih, std::list<unsigned long>& adj_hits);
    void GradientTest(TrkPlaneClusters::Cluster& cls);
    Int_t SignI(Int_t i) {return (0<i)-(i<0);}

    TrkPlaneBasicHits::Coll& fHits; // input hits

    TrkPlaneClusters::Coll fClusters; // output clusters

    std::string fNam; // detector name

    //selection criteria for clusters
    Double_t fClsMinLimMdist=0; // limit on minimal distance to another cluster, mm

    //counters for processed clusters
    unsigned long fNall=0;
    unsigned long fNsel=0;

};

#endif

