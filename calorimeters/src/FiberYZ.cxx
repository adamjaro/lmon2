
//_____________________________________________________________________________
//
// Helper class for optic fiber bent in y-z plane
//
//_____________________________________________________________________________

//C++
#include <vector>
#include <iostream>

//ROOT
#include "TMath.h"

//local classes
#include "FiberYZ.h"

using namespace std;

//_____________________________________________________________________________
FiberYZ::FiberYZ(Double_t L, Double_t yL, Double_t r): fL(L), fyL(yL), fR(r) {

  //axis polynomial
  fA = -2.*yL/(L*L*L);
  fB = 3.*yL/(L*L);

  //initial slices along z
  Int_t nz = 12;
  Double_t dz = L/(nz-1);

  Double_t rmin = 0.3;
  Double_t rminsq = rmin*rmin;

  //finder loop
  bool stat = true;
  while(stat) {

    //z-loop
    for(Int_t iz=0; iz<nz; iz+=2) {

      Double_t z1 = iz*dz;
      Double_t z2 = (iz+1)*dz;

      Double_t y1 = f(z1);
      Double_t y2 = f(z2);

      Double_t a1, b1, a2, b2;
      MakeNormVect(z1, a1, b1);
      MakeNormVect(z2, a2, b2);

      Double_t yz1[2], yz2[2];
      yz1[0] = z1 + a1*r;
      yz1[1] = y1 + b1*r;
      yz2[0] = z2 + a2*r;
      yz2[1] = y2 + b2*r;

      Double_t rsq1 = (yz2[1]-yz1[1])*(yz2[1]-yz1[1]) + (yz2[0]-yz1[0])*(yz2[0]-yz1[0]);

      if(rsq1 < rminsq) {
        stat = false;
        break;
      }

      yz1[0] = z1 - a1*r;
      yz1[1] = y1 - b1*r;
      yz2[0] = z2 - a2*r;
      yz2[1] = y2 - b2*r;

      Double_t rsq2 = (yz2[1]-yz1[1])*(yz2[1]-yz1[1]) + (yz2[0]-yz1[0])*(yz2[0]-yz1[0]);

      if(rsq2 < rminsq) {
        stat = false;
        break;
      }

    //cout << nz << " " << z1 << " " << z2 << " " << rsq1 << " " << rsq2 << endl;
    }//z-loop

    nz += 1;
    dz = L/(nz-1);

    //cout << nz << " " << dz << endl;

  }//finder loop

  cout << "nz: " << nz << " dz: " << dz << endl;

  Int_t nphi = 12;
  Double_t dphi = TMath::Pi()/nphi;

  //slice loop
  for(Int_t i=0; i<nz; i++) {
    Double_t z = i*dz;

    fSlices.push_back(slice(z, nphi, dphi, this));
  }//slice loop

}//FiberYZ

//_____________________________________________________________________________
Double_t FiberYZ::f(Double_t z) {

  //fiber axis f(z)

  return fA*z*z*z + fB*z*z;

}//f

//_____________________________________________________________________________
Double_t FiberYZ::fP(Double_t z) {

  //derivative f'(z)

  return 3*fA*z*z + 2*fB*z;

}//f

//_____________________________________________________________________________
void FiberYZ::MakeNormVect(Double_t z, Double_t& a, Double_t& b) {

  //perpendicular unit vector to the axis

  Double_t fPz = fP(z);
  Double_t a1 = TMath::Sqrt(1./(1+fPz*fPz));
  Double_t b1 = a1*fPz;

  a = -b1;
  b = a1;

}//MakeNormVect

//_____________________________________________________________________________
FiberYZ::slice::slice(Double_t zz, Int_t nphi, Double_t dphi, FiberYZ *fi): z(zz), fib(fi) {

  //slice center
  y = fib->f(z);
  r = fib->fR;

  //normal unit vector
  fib->MakeNormVect(z, a, b);

  //cout << a << " " << b << endl;

  //points on the slice
  //for(Int_t ip=0; ip<nphi; ip++) {

  //}
  points.push_back({0, r, 0});
  points.push_back({0, -r, 0});

  //rotate the points in y-z
  for(array<Double_t, 3>& i: points) {
    Double_t zp = i[2];
    Double_t yp = i[1];

    i[2] = zp*b + yp*a;
    i[1] = -zp*a + yp*b;
  }

  //translate the points to the given y and z
  for(array<Double_t, 3>& i: points) {

    i[2] += z;
    i[1] += y;
  }

  //cout << "slice " << zz << " " << y << " " << r << endl;

}//slice















