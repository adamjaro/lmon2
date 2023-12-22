
#ifndef CalPWOClusterWavg_h
#define CalPWOClusterWavg_h

//clusterizer for CalPWO based on weighted average over all hits

class GeoParser;

#include "CalPWOHits.h"
#include "CaloCluster.h"

class CalPWOClusterWavg {

  public:

    CalPWOClusterWavg(std::string nam);

    void ConnectInput(TTree *tree);
    void SetGeometry(std::string geo_nam, std::string det_nam, GeoParser *geo);

    void CreateOutput(TTree *tree);

    void ProcessEvent();

    //access to created clusters during reconstruction
    const std::vector<CaloCluster::Cls>& GetClusters() { return fCls.GetStore(); }

  private:

    std::string fNam; // CalPWO name

    CalPWOHits::Coll fHits; // input hits

    G4double fXpos=0; // calorimeter position in x, mm
    G4double fYpos=0; // calorimeter position in y, mm
    G4double fZpos=0; // calorimeter position in z, mm
    G4double fThetaX=0; // calorimeter rotation along x, rad
    G4double fThetaY=0; // calorimeter rotation along y, rad

    CaloCluster::Coll fCls; // reconstructed clusters

};

#endif

