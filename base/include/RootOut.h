#ifndef RootOut_h
#define RootOut_h

//helper class for ROOT TTree output

class TFile;
class TTree;

class G4GenericMessenger;

class RootOut {

  public:

    RootOut();

    void Open(); // open the output and create the tree
    void Close(); // write the tree and close outpu file

    TTree *GetTree() {return fDetTree;} // provide the tree to the detectors
    void FillTree(); // call fill for the tree

  private:

    G4String fOutputName; // name of output file
    G4GenericMessenger *fMsg; // messenger for name of output file

    TFile *fOut; // output file
    TTree *fDetTree; // output tree

};

#endif

