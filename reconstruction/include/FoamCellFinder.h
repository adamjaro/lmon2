
#ifndef FoamCellFinder_h
#define FoamCellFinder_h

// Variable binning based on TFoam

#include "TFoam.h"
#include "TFoamIntegrand.h"
#include "TRandom3.h"

class FoamCellFinder: public TFoam {


  public:

    FoamCellFinder(const Char_t* Name): TFoam(Name) {}
    FoamCellFinder(const TH1D& hx, Int_t ncells);

    Int_t GetLastCe() { return fLastCe; }
    const std::vector<Long_t>& GetActiveIdx() { return fCellsAct; }
    std::vector<Long_t> GetActiveIdxSort();

    TFoamCell* GetCell(std::size_t i) { return fCells[i]; }

    void GetCellPos1D(std::size_t i, Double_t& pos, Double_t& siz);

    TH1D MakeH1D(const char *name_title);

    TH1D GetInpH1() { return fInpHist; }

  private:

    Double_t fInpMin;
    Double_t fInpMax;
    TH1D fInpHist;

    TRandom3 fRand;

    class FuncInpH1: public TFoamIntegrand {
      TH1D *hx;

      public:

      void SetH1(TH1D *x) { hx = x; }

      Double_t Density(Int_t n, Double_t *x) {
        return hx->Interpolate(x[0]);
      }
    };

    FuncInpH1 fFuncInpH1;

};


#endif

