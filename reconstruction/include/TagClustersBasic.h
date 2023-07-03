
#ifndef TagClustersBasic_h
#define TagClustersBasic_h

class GeoParser;

#include<unordered_map>
#include <list>

#include "DetectorData.h"

#include "TrkPlaneBasicHits.h"

class TagClustersBasic {

  public:

    TagClustersBasic(std::string nam);

    void SetGeometry(GeoParser *geo);
    void ConnectHitsInput(TTree *tree) { fHits.ConnectInput("lowQ2_"+fNam, tree); }

    void CreateOutput(TTree *tree) { fClusters.CreateOutput(fNam+"_clusters", tree); }

    void ConnectClustersInput(TTree *) {}

    void ProcessEvent();
    void FinishEvent();

    void SetLimMdist(Double_t d) { fClsMinLimMdist = d; }
    Double_t GetLimMdist() { return fClsMinLimMdist; }

    //cluster representation
    struct Cluster {

      Double_t GetSigma(Double_t swx2, Double_t pos);

      Double_t x=0; // cluster x position, mm
      Double_t y=0; // cluster y position, mm
      Double_t en=0; // cluster energy, keV
      Int_t nhits=0; // number of hits for the cluster
      Bool_t is_prim=1; // flag for primary particle
      Double_t sigma_x=0; // uncertainty in cluster x position, mm
      Double_t sigma_y=0; // uncertainty in cluster y position, mm
      Int_t itrk=-1; // MC track index associated with the cluster
      Int_t pdg=0; // PDG code for the MC track
      Int_t prim_id=0; // ID of primary particle associated with the cluster
      Int_t ntrk=0; // number of tracks for which the cluster was used
      Double_t min_dist=-1; // minimal distance to another cluster, mm
      Int_t id=0; // cluster ID on the plane
      Int_t iplane=0; // plane ID
      Bool_t stat=kTRUE; // cluster status

      std::list<unsigned long> hits; // indices for hits contributing to cluster

    };//Cluster

    //clusters collection
    class Coll : public DetectorData<Cluster> {

      public:

        Coll();

        std::vector<Cluster>& GetStore() { return fStorage; }

    };//Coll

    Coll& GetClusters() { return fClusters; } 

  private:

    int GetHitsCount();
    unsigned long FindHitEmax();
    int FindAdjHits(unsigned long ih, std::list<unsigned long>& adj_hits);
    void GradientTest(Cluster& cls);
    Int_t SignI(Int_t i) {return (0<i)-(i<0);}

    std::string fNam; // plane name

    TrkPlaneBasicHits::Coll fHits; // input hits

    Coll fClusters; // output clusters

    std::unordered_map<unsigned long, bool> fHitStat; // status flags for hits in plane

    //selection criteria
    Double_t fEmin; // keV, threshold in energy
    Double_t fClsMinLimMdist=0; // limit on minimal distance to another cluster, mm

    G4double fXpos=0; // plane position in x, mm
    G4double fYpos=0; // plane position in y, mm
    G4double fZpos=0; // plane position in z, mm
    G4double fThetaX=0; // plane rotation along x, rad
    G4double fThetaY=0; // plane rotation along y, rad

};

#endif

