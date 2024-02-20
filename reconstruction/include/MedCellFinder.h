
#ifndef MedCellFinder_h
#define MedCellFinder_h

// variable binning based on median cell split

#include "TF1.h"

class MedCellFinder {

  public:

    MedCellFinder(TH1D *hx, std::size_t ncells);
    TH1D MakeH1D(const char *name_title);

    TF1& GetFuncInpH1() { return fFuncInpH1; }
    void PrintCells();

    class Cell {
      public:

      Double_t xmin=0;
      Double_t xmax=0;

      bool is_active=false;
    };

  private:

    Double_t EvalInpH1(Double_t *x, Double_t*);

    TH1D *fInpH1;
    TF1 fFuncInpH1;

    std::vector<Cell> fCells;
    std::size_t fLastCell=0;
    std::vector<size_t> fActiveCellIdx;



};

#endif

