
#ifndef TagTracks_h
#define TagTracks_h

// Tracks collection

#include "DetectorData.h"

namespace TagTracks {

    //track representation
  struct Track {

    Double_t x=0; // track position in x, mm
    Double_t y=0; // track position in y, mm
    Double_t slope_x=0; // track slope in x
    Double_t slope_y=0; // track slope in y
    Double_t theta_x=0; // track angle along x, rad
    Double_t theta_y=0; // track angle along y, rad
    Double_t chi2_x=0; // track chi^2 in x
    Double_t chi2_y=0; // track chi^2 in y
    Double_t chi2_xy=0; // track chi^2 in xy plane
    Bool_t is_prim=0; // track for primary particle
    Int_t itrk=-1; // index for MC particle
    Int_t pdg=0; // pdg for MC particle
    Int_t prim_id=-1; // ID of primary particle associated with the track
    Int_t num_diff_itrk=0; // number of unique MC track indices in the clusters
    Bool_t is_rec = 0; // reconstruction flag, 1 = track is reconstructed
    Double_t rec_en = 0; // reconstructed electron energy, GeV
    Double_t rec_theta = 0; // electron polar angle, rad
    Double_t rec_phi = 0; // electron azimuthal angle, rad
    Double_t rec_en_err = 0; // error in reconstructed electron energy, GeV
    Double_t rec_theta_err = 0; // error in electron polar angle, rad
    Double_t rec_phi_err = 0; // error in electron azimuthal angle, rad
    Double_t rec_Q2 = 0; // reconstructed electron Q^2, GeV^2
    Int_t ninp = 0; // number of inputs used for reconstruction
    Int_t ilay = -1; // layer index from reconstruction
    Bool_t has_cal=0; // flag for matched calorimeter cluster
    Double_t cal_x=0; // matched calorimeter cluster, x position (mm)
    Double_t cal_y=0; // matched calorimeter cluster, y position (mm)
    Double_t cal_en=0; // matched calorimeter cluster, energy (GeV)
    Double_t cal_ext_x=0; // extrapolated x track position to calorimeter (mm)
    Double_t cal_ext_y=0; // extrapolated x track position to calorimeter (mm)
    Bool_t has_mcp = 0; // track is paired with MC particle
    Double_t mcp_en = 0; // MC particle energy, GeV
    Double_t mcp_theta = 0; // MC particle polar angle, rad
    Double_t mcp_phi = 0; // MC particle azimuthal angle, rad
    Double_t mcp_Q2 = 0; // MC particle Q^2, GeV^2

    //std::vector<TagClustersBasic::Cluster*> cls; // track clusters

    //void ExtrapolateZ(Double_t z, Double_t& xe, Double_t& ye);

  };//Track

    //tracks collection
    class Coll : public DetectorData<Track> {

      public:

        Coll();

        std::vector<Track>& GetStore() { return fStorage; }

    };//Coll

}

#endif

