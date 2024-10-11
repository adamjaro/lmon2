
#ifndef TrkPlaneClusters_h
#define TrkPlaneClusters_h

// Cluster collection for a tracking plane

#include <list>

#include "DetectorData.h"
#include "TMath.h"

namespace TrkPlaneClusters {

  //cluster representation
  struct Cluster {

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

    Double_t GetSigma(Double_t swx2, Double_t pos) {

      //Bevington, Robinson, eq. 4.22, p. 58

      if( hits.size() < 2 ) return 0;

      Double_t sig2 = ((swx2/en) - pos*pos)/(hits.size() - 1);

      if(sig2 > 0) return TMath::Sqrt(sig2);

      return 0;
    }

  };//Cluster

  //clusters collection
  class Coll : public DetectorData<Cluster> {

    public:

      Coll() {

        DATA_ADD_UNIT_ATTR( x )
        DATA_ADD_UNIT_ATTR( y )
        DATA_ADD_UNIT_ATTR( en )
        DATA_ADD_UNIT_ATTR( nhits )
        DATA_ADD_UNIT_ATTR( is_prim )
        DATA_ADD_UNIT_ATTR( sigma_x )
        DATA_ADD_UNIT_ATTR( sigma_y )
        DATA_ADD_UNIT_ATTR( ntrk )
        DATA_ADD_UNIT_ATTR( min_dist )
        DATA_ADD_UNIT_ATTR( stat )

      }

      std::vector<Cluster>& GetStore() { return fStorage; }

  };//Coll

}

#endif

