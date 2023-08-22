
#ifndef TagTrackFindBasic_h
#define TagTrackFindBasic_h

// Track finder for clusters by TagClustersBasic

class TagClustersBasic;
class GeoParser;

#include "DetectorData.h"

#include "TagClustersBasic.h"
#include "Math/Point3D.h" 
#include "Math/Vector3D.h" 
#include "Math/RotationZYX.h" 

class TagTrackFindBasic {

  public:

    TagTrackFindBasic(std::string nam);

    void SetGeometry(GeoParser *geo);
    void ConnectHitsInput(TTree *tree) { for(TagClustersBasic *i: fPlanes) i->ConnectHitsInput(tree); }

    void CreateOutput(TTree *tree, bool create_planes=true, bool create_tracks=true);

    void ConnecTracksInput(TTree *) {}

    void ProcessEvent();
    void FinishEvent();

    void SetMaxChi2Ndf(Double_t max_chi2ndf) { fChi2ndfMax = max_chi2ndf; }
    void SetClsLimMdist(Double_t d);

    std::string GetName() { return fNam; }

    Double_t GetMaxChi2ndf() { return fChi2ndfMax; }
    Double_t GetClsLimMdist() { return fPlanes[0]->GetLimMdist(); }

    ROOT::Math::XYZVector getOffset(){return detectorOffset;}
    void setOffset(double x, double y, double z){ detectorOffset = ROOT::Math::XYZPoint(x,y,z);}

    double getAngle()        { return detectorAngle;}
    void   setAngle(double angle){ detectorAngle = angle;}

    TagClustersBasic* GetPlane(int i) { return fPlanes[i]; }
    Double_t GetPlaneZ(int i) { return fZ[i]; }

    int GetNumberOfClusters();

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
      Bool_t is_associate=0; // track association to a reference MC particle
      Double_t ref_x=0; // reference position in x, mm
      Double_t ref_y=0; // reference position in y, mm
      Double_t ref_theta_x=0; // reference angle along x, rad
      Double_t ref_theta_y=0; // reference angle along y, rad
      Int_t evt_ntrk=0; // number of all tracks in event for a given track
      Int_t num_shared_cls=0; // number of track clusters shared with another track
      Int_t num_diff_itrk=0; // number of unique MC track indices in the clusters
      Bool_t is_rec = 0; // reconstruction flag, 1 = track is reconstructed
      Double_t rec_en = 0; // reconstructed electron energy, GeV
      Double_t rec_theta = 0; // electron polar angle, rad
      Double_t rec_phi = 0; // electron azimuthal angle, rad
      Double_t rec_Q2 = 0; // reconstructed electron Q^2, GeV^2
      Double_t rec2_en = 0; // reconstructed electron energy, GeV
      Double_t rec2_theta = 0; // electron polar angle, rad
      Double_t rec2_phi = 0; // electron azimuthal angle, rad
      Double_t rec2_Q2 = 0; // reconstructed electron Q^2, GeV^2
      Int_t ninp = 0; // number of inputs used for reconstruction
      Int_t ilay = -1; // layer index from reconstruction
      Bool_t has_mcp = 0; // track is paired with MC particle
      Double_t mcp_en = 0; // MC particle energy, GeV
      Double_t mcp_theta = 0; // MC particle polar angle, rad
      Double_t mcp_phi = 0; // MC particle azimuthal angle, rad
      Double_t mcp_Q2 = 0; // MC particle Q^2, GeV^2
      Double_t true_en=0; // true generated electron energy, GeV
      Double_t true_theta=0; // true generated electron polar angle, rad
      Double_t true_phi=0; // true generated electron azimuthal angle, rad
      Double_t true_Q2=0; // true generated event Q^2, GeV^2

      std::vector<TagClustersBasic::Cluster*> cls; // track clusters

    };//Track

    //tracks collection
    class Coll : public DetectorData<Track> {

      public:

        Coll();

        std::vector<Track>& GetStore() { return fStorage; }

    };//Coll

    std::vector<Track>& GetTracks() { return fTracks.GetStore(); }

  private:

    void MakeTrack(Double_t *x, Double_t& pos, Double_t& slope, Double_t& theta, Double_t& chi2);
    Double_t TrackChi2(Double_t *x, Double_t *y, Track& trk);

    template<std::size_t N> void ClusterAnalysis(TagClustersBasic::Cluster* (&cls)[N], Track& trk);

    std::vector<TagClustersBasic*> fPlanes; // planes for the station

    Coll fTracks; // tracks collection

    Double_t fChi2ndfMax=0; // maximal reduced chi2 for tracks

    Double_t fL=0; // plane spacing, mm
    Double_t fZ[4] = {0,0,0,0}; // local z positions for planes, mm
    ROOT::Math::XYZVector  detectorOffset{0,0,0};
    double detectorAngle{0};

    std::string fNam;

};

#endif

















