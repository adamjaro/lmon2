
//Helper class for optic fiber bent in y-z plane

#ifndef FiberYZ_h
#define FiberYZ_h

class FiberYZ {

  public:

    FiberYZ(Double_t L, Double_t yL, Double_t r);

    Double_t f(Double_t z);
    Double_t fP(Double_t z);

    class slice {
      public:

      slice(Double_t zz, Int_t nphi, Double_t dphi, FiberYZ *fi);

      Double_t z; // slice center in y-z
      Double_t y;
      Double_t r; // radius
      Double_t a; // tangent unit vector
      Double_t b;

      //points on the slice
      std::vector<std::array<Double_t, 3>> points;

      private:
      FiberYZ *fib;
    };

    size_t GetNSlice() { return fSlices.size(); }
    slice& GetSlice(size_t i) { return fSlices[i]; }

    Double_t GetR() const { return fR; }

    void MakeNormVect(Double_t z, Double_t& a, Double_t& b);

  private:

    Double_t fL; // length along z, mm
    Double_t fyL; // elevation at z = L, mm
    Double_t fR; // fiber radius, mm
    Double_t fA, fB; // polynomial coeffitients

    std::vector<slice> fSlices;

};

#endif

