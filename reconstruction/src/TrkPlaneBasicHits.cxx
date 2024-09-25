
//_____________________________________________________________________________
//
// Acceptance analysis for tagger detectors with Timepix4
//
//_____________________________________________________________________________

//C++
#include <iostream>
#include <array>

//ROOT
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

//local classes
#include "TrkPlaneBasicHits.h"
#include "TagTpix4Acceptance.h"

using namespace std;

//_____________________________________________________________________________
void TagTpix4Acceptance::Run() {

  cout << "Tpix4HitAnalysis::Run" << endl;

  //input tree
  TChain tree("DetectorTree");
  tree.Add("lmon.root");

  //true event kinematics
  tree.SetBranchAddress("true_Q2", &fTrueQ2);

  TFile out("acc.root", "recreate");

  TTree otree("event", "event");
  otree.Branch("true_Q2", &fTrueQ2, "true_Q2/D");

  Bool_t s1_allhit;
  otree.Branch("s1_allhit", &s1_allhit, "s1_allhit/O");

  //plane hits
  std::array<plane_hits, 4> ahits{
    plane_hits("lowQ2_s1_1", &tree),
    plane_hits("lowQ2_s1_2", &tree),
    plane_hits("lowQ2_s1_3", &tree),
    plane_hits("lowQ2_s1_4", &tree)
  };

  Long64_t nev = tree.GetEntries();
  nev = 12;

  cout << "Number of events: " << nev << endl;

  //event loop
  for(Long64_t iev=0; iev<nev; iev++) {

    cout << "iev: " << iev << endl;

    tree.GetEntry(iev);

    unsigned long ns1, ns2, ns3, ns4;
    ns1 = ahits[0].GetNhits();
    ns2 = ahits[1].GetNhits();
    ns3 = ahits[2].GetNhits();
    ns4 = ahits[3].GetNhits();

    s1_allhit = kFALSE;
    if( ns1!=0 and ns2!=0 and ns3!=0 and ns4!=0 ) { s1_allhit = kTRUE; }

    otree.Fill();

  }//event loop

  otree.Write(0, TObject::kOverwrite);

  out.Close();

  cout << "All done" << endl;

}//Run

//_____________________________________________________________________________
TagTpix4Acceptance::plane_hits::plane_hits(std::string nam, TTree *in_tree) {

  //input hits
  hits = make_unique<TrkPlaneBasicHits::Coll>();
  hits->ConnectInput(nam, in_tree);

}//plane_hits

//_____________________________________________________________________________
unsigned long TagTpix4Acceptance::plane_hits::GetNhits() {

  hits->LoadInput();

  return hits->GetN();

}//GetNhits























