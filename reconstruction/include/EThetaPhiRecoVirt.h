
#ifndef EThetaPhiRecoVirt_h
#define EThetaPhiRecoVirt_h

// abstract base for EThetaPhiReco versions, remains left unused
// because virtual functions don't work with CDLL in ctypes
// nor with gSystem in ROOT

class EThetaPhiRecoVirt {

  public:

    virtual ~EThetaPhiRecoVirt() {};

    virtual void Initialize(boost::program_options::variables_map*);

    virtual void AddInput(Double_t*, Double_t, Double_t, Double_t);
    virtual void Finalize();
    virtual void Export();

    virtual Bool_t Reconstruct(Double_t*, Double_t&, Double_t&, Double_t&, Int_t*) = 0;
    virtual void Import(TFile*);
    virtual void SetMinNinp(Int_t);

};

#endif

