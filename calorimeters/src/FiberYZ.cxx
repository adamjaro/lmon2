
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
FiberYZ::FiberYZ(Double_t L, Double_t yL, Double_t r, Double_t rmin): fL(L), fyL(yL), fR(r) {

  //axis polynomial
  fA = -2.*yL/(L*L*L);
  fB = 3.*yL/(L*L);

  //initial slices along z
  Int_t nz = 12;
  Double_t dz = L/(nz-1);

  Double_t rminsq = rmin*rmin;

  //finder loop
  bool stat = true;
  while(true) {

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

    if( !stat ) break;

    nz += 1;
    dz = L/(nz-1);

    //cout << nz << " " << dz << endl;

  }//finder loop

  //segmentation in y-x plane
  Int_t nphi = TMath::Pi()/TMath::ASin( rmin/(2*r) );
  Double_t dphi = 2*TMath::Pi()/nphi;

  cout << "nz: " << nz << " dz: " << dz << " nphi: " << nphi << endl;

  //slice loop
  for(Int_t i=0; i<nz; i++) {
    Double_t z = i*dz;

    fSlices.push_back(slice(z, nphi, dphi, this));
  }//slice loop

  //output facets
  MakeFacets();

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
void FiberYZ::InvertZ() {

  //invert the shape along z

  //slice loop
  for(slice& slc: fSlices) {

    //center in z
    slc.z *= -1;
    slc.z += fL;

    //point loop
    for(array<Double_t, 3>& p: slc.points) {

      p[2] *= -1;
      p[2] += fL;
    }//point loop
  }//slice loop

  //recalculate the facets
  MakeFacets();

}//InvertZ

//_____________________________________________________________________________
void FiberYZ::MakeFacets() {

  //construct triangular facets representing the shape

  fFacets.clear();

  //slice loop
  for(size_t is=0; is<(fSlices.size()-1); is++) {

  //size_t is = 1;
  //size_t ip = 2;

    const slice& s0 = fSlices[is];
    const slice& s1 = fSlices[is+1];

    //points loop
    for(size_t ip=0; ip<(s0.points.size()-1); ip++) {

      const array<Double_t, 3>& p0 = s0.points[ip];
      const array<Double_t, 3>& p1 = s0.points[ip+1];
      const array<Double_t, 3>& p2 = s1.points[ip+1];
      const array<Double_t, 3>& p3 = s1.points[ip];

      //cout << "p0: " << p0[0] << " " << p0[1] << " " << p0[2] << endl;
      //cout << "p1: " << p1[0] << " " << p1[1] << " " << p1[2] << endl;
      //cout << "p2: " << p2[0] << " " << p2[1] << " " << p2[2] << endl;
      //cout << "p3: " << p3[0] << " " << p3[1] << " " << p3[2] << endl;

      //two triangles, first
      fFacets.push_back( facet() );
      facet& fct = fFacets.back();
      fct.p0 = p0;
      fct.p1 = p1;
      fct.p2 = p2;

      //two triangles, second
      fFacets.push_back( facet() );
      facet& fct1 = fFacets.back();
      fct1.p0 = p2;
      fct1.p1 = p3;
      fct1.p2 = p0;

    }//points loop

    //closing triangles
    const array<Double_t, 3>& p0 = s0.points.back();
    const array<Double_t, 3>& p1 = s0.points[0];
    const array<Double_t, 3>& p2 = s1.points[0];
    const array<Double_t, 3>& p3 = s1.points.back();

      //two closing triangles, first
      fFacets.push_back( facet() );
      facet& fct = fFacets.back();
      fct.p0 = p0;
      fct.p1 = p1;
      fct.p2 = p2;

      //two closing triangles, second
      fFacets.push_back( facet() );
      facet& fct1 = fFacets.back();
      fct1.p0 = p2;
      fct1.p1 = p3;
      fct1.p2 = p0;

  }//slice loop

  //closing facets at first slice
  const slice& s0 = fSlices[0];
  //points loop
  for(size_t ip=0; ip<(s0.points.size()-1); ip++) {

    const array<Double_t, 3>& p0 = s0.points[ip];
    const array<Double_t, 3>& p1 = s0.points[ip+1];

    fFacets.push_back( facet() );
    facet& fct = fFacets.back();
    fct.p0 = {0,s0.y,s0.z};
    fct.p1 = p1;
    fct.p2 = p0;
  }//points loop

  //final facet, first slice
  const array<Double_t, 3>& p00 = s0.points.back();
  const array<Double_t, 3>& p01 = s0.points[0];
  fFacets.push_back( facet() );
  facet& fct0 = fFacets.back();
  fct0.p0 = {0,s0.y,s0.z};
  fct0.p1 = p01;
  fct0.p2 = p00;

  //closing facets at last slice
  const slice& s1 = fSlices.back();
  //points loop
  for(size_t ip=0; ip<(s1.points.size()-1); ip++) {

    const array<Double_t, 3>& p0 = s1.points[ip];
    const array<Double_t, 3>& p1 = s1.points[ip+1];

    fFacets.push_back( facet() );
    facet& fct = fFacets.back();
    fct.p0 = {0,s1.y,s1.z};
    fct.p1 = p0;
    fct.p2 = p1;
  }//points loop

  //final facet, last slice
  const array<Double_t, 3>& p0 = s1.points.back();
  const array<Double_t, 3>& p1 = s1.points[0];
  fFacets.push_back( facet() );
  facet& fct = fFacets.back();
  fct.p0 = {0,s1.y,s1.z};
  fct.p1 = p0;
  fct.p2 = p1;

}//MakeFacets

//_____________________________________________________________________________
FiberYZ::slice::slice(Double_t zz, Int_t nphi, Double_t dphi, FiberYZ *fi): z(zz), fib(fi) {

  //slice center
  y = fib->f(z);
  r = fib->fR;

  //normal unit vector
  fib->MakeNormVect(z, a, b);

  //cout << a << " " << b << endl;
  //cout << "slice:" << endl;

  //points on the slice
  for(Int_t ip=0; ip<nphi; ip++) {

    Double_t phi = ip*dphi;
    Double_t xs = r*TMath::Cos(phi);
    Double_t ys = r*TMath::Sin(phi);

    //cout << ip << " " << phi << " " << xs << " " << ys << endl;

    points.push_back({xs, ys, 0});
  }

  //points.push_back({0, r, 0});
  //points.push_back({0, -r, 0});

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















