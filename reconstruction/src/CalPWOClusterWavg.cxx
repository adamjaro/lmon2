
//_____________________________________________________________________________
//
// Clusterizer for CalPWO based on weighted average over all hits
//
// Works correct only for single-track events since
// all hits are used for one cluster per event, aimed for development
//
//_____________________________________________________________________________

//C++
#include <iostream>

//ROOT
#include "TTree.h"
#include "TMath.h"

//Geant
#include "G4String.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//local classes
#include "CalPWOClusterWavg.h"
#include "GeoParser.h"

using namespace std;

//_____________________________________________________________________________
CalPWOClusterWavg::CalPWOClusterWavg(string nam): fNam(nam) {

}//CalPWOClusterWavg

//_____________________________________________________________________________
void CalPWOClusterWavg::ConnectInput(TTree *tree) {

  fHits.ConnectInput(fNam, tree);

}//ConnectInput

//_____________________________________________________________________________
void CalPWOClusterWavg::SetGeometry(string geo_nam, std::string det_nam, GeoParser *geo) {

  //position of volume where the CalPWO is placed
  geo->GetOptD(geo_nam, "xpos", fXpos, GeoParser::Unit(mm));
  geo->GetOptD(geo_nam, "ypos", fYpos, GeoParser::Unit(mm));
  geo->GetOptD(geo_nam, "zpos", fZpos, GeoParser::Unit(mm));

  //rotation of volume where the CalPWO is placed
  G4double theta = 0;
  geo->GetOptD(geo_nam, "theta", theta, GeoParser::Unit(rad));

  G4bool rotate_x = false;
  geo->GetOptB(geo_nam, "rotate_x", rotate_x);

  if( rotate_x ) {
    fThetaX = theta;
  } else {
    fThetaY = theta;
  }

  //shift in x position inside the volume where the CalPWO is placed
  G4double shift_x=0;
  geo->GetOptD(det_nam, "xpos", shift_x, GeoParser::Unit(mm));
  fXpos += shift_x;

  cout << "xyz: " << fXpos << " " << fYpos << " " << fZpos << " " << fThetaX << " " << fThetaY << endl;

}//SetGeometry

//_____________________________________________________________________________
void CalPWOClusterWavg::CreateOutput(TTree *tree) {

  fCls.CreateOutput(fNam+"_clusters", tree);

}//CreateOutput

//_____________________________________________________________________________
void CalPWOClusterWavg::ProcessEvent() {

  fHits.LoadInput();
  fCls.ClearEvent();

  if(fHits.GetN() <= 0) return;

  //cout << "nhits: " << fHits.GetN() << endl;

  //cluster position (mm), energy (GeV) and number of contributors
  Double_t cls_x=0, cls_y=0, cls_en=0;
  Int_t cls_nhit=0;

  for(CalPWOHits::Hit& i: fHits.GetReadData()) {

    //transformation to local coordinates
    i.Translate(-fXpos, -fYpos, -fZpos);
    i.RotateXY(-fThetaX, -fThetaY);

    //position
    cls_x += i.en*i.x;
    cls_y += i.en*i.y;

    //energy and contributors
    cls_en += i.en;
    cls_nhit++;

    //cout << "hit: " << i.cell_id << " " << i.x << " " << i.y << " " << i.z << " " << i.en << " " << i.prim_id << endl;

  }

  //finalize the weighted average
  cls_x = cls_x/cls_en;
  cls_y = cls_y/cls_en;

  //add the cluster
  CaloCluster::Cls& cls = fCls.Add();
  cls.x = cls_x;
  cls.y = cls_y;
  cls.en = cls_en;
  cls.nhit = cls_nhit;

  //cout << "cluster: " << cls_x << " " << cls_y << " " << cls_en << " " << cls_nhit << endl;

  fCls.FinishEvent();

}//ProcessEvent



























