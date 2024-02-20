
#ifndef SumCellFinder_h
#define SumCellFinder_h

// cell finder based on even sum per bin from interpolated input histogram

#include "TH1D.h"

class SumCellFinder: public TH1D {

  public:

    SumCellFinder(TH1D *hx, const char *name_title, int ndiv=1e6);

};

#endif

