#include "pSphere.h"
#include "math.h"
#include<iostream>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor, copy and destructor---------------------------------------------
pSphere::pSphere(unsigned int dim, double* center, double radius) {
  p = dim;
  r = radius;
  c = new double[p];
  for (unsigned int i = 0; i < p; i++) {
    c[i] = center[i];
  }
}

pSphere::pSphere(const pSphere &sphere) {
  p = sphere.p;
  r = sphere.r;
  c = new double[p];
  for (unsigned int i = 0; i < p; i++) {
    c[i] = sphere.c[i];
  }
}

pSphere::~pSphere() {
  delete [] c;
  c = NULL;
}

//accessory---------------------------------------------------------------------
unsigned int  pSphere::get_p() const { return p; }
double pSphere::get_r() const { return r; }
double* pSphere::get_c() const { return c; }

//tools-------------------------------------------------------------------------
double pSphere::get_dist(double* pnt1, double* pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < p; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);
}
void pSphere::idpSphere(unsigned int dim, double* center, double radius) {
  p = dim;
  r = radius;
  for (unsigned int i = 0; i < p; i++) {
    c[i] = center[i];
  }
}

void pSphere::createSphere(unsigned int dim, unsigned int i, unsigned int t, double** &csdY, double** &csdY2, double* &lCosts) {
  pCost cost = pCost(dim);
  cost.idpCost( p, i, t, csdY, csdY2, lCosts);
  double r2 = (lCosts[t+1] - lCosts[i] - cost.get_kVYit()) / cost.get_k();
  if (r2 > 0) {
    idpSphere( p, cost.get_EYit(), sqrt(r2));
  } else {
    idpSphere(p,cost.get_EYit(), 0);
  }
}

bool pSphere ::isIntersection(const pSphere  &sphere) {
  if (get_dist(c, sphere.get_c()) >= (r + sphere.get_r())) {
    return false;
  } else {
    return true;
  }
}

bool pSphere ::isnotIntersection(const pSphere  &sphere) {
  if (get_dist(c, sphere.get_c()) < (r + sphere.get_r())) {
    return false;
  } else {
    return true;
  }
}

bool pSphere ::isInclusion(const pSphere  &sphere) {
  if (get_dist(c, sphere.get_c()) <= (sphere.get_r() - r)) {
    return true;
  } else {
    return false;
  }
}
//----------------------------------------------------------------------------//
