#include "pCost.h"
#include <iostream>
#include "math.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor, constructor copy, destructor-------------------------------------
pCost::pCost(unsigned int dim) {
  p = dim;
  k = 0;
  mi1beta = 0;
  kVYit = 0;
  EYit = new double[p];
}

pCost::pCost (unsigned int dim, unsigned int i, unsigned int t, double** &csdY, double** &csdY2, double* &lCosts) {
  p = dim;
  k = t - i + 1;
  mi1beta = lCosts[i];
  EYit = new double[p];

  double sum_EYit2 = 0;
  double sum_dif_Yit2 = 0;

  for (unsigned int j = 0; j < p; j++) {
    EYit[j] = (csdY[t+1][j] - csdY[i][j])/k;

    sum_EYit2 = sum_EYit2 + EYit[k] * EYit[k];
    sum_dif_Yit2 = sum_dif_Yit2 + (csdY2[t+1][k] - csdY2[i][k]);
  }
  kVYit = sum_dif_Yit2 - k * sum_EYit2;
}

pCost::pCost(const pCost &cost) {
  p = cost.p;
  k = cost.k;
  mi1beta = cost.mi1beta;
  kVYit =cost.kVYit;
  EYit = new double[p];
  for (unsigned int j = 0; j < p; j++) {EYit[j] = cost.EYit[j];}
}

pCost::~pCost() {
  delete [] EYit;
  EYit = NULL;
}

//accessory---------------------------------------------------------------------
unsigned int pCost::get_p() const { return p; }
unsigned int pCost::get_k() const { return k; }
double pCost::get_kVYit() const { return kVYit; }
double pCost::get_mi1beta() const { return mi1beta; }
double* pCost::get_EYit() { return EYit; }

//tools-------------------------------------------------------------------------
void pCost::idpCost(unsigned int dim, unsigned int i, unsigned int t, double** &csdY, double** &csdY2, double* &lCosts) {
  p = dim;
  k = t - i + 1;
  mi1beta = lCosts[i];

  double sum_EYit2 = 0;
  double sum_dif_Yit2 = 0;

  for (unsigned int j = 0; j < p; j++) {
    EYit[j] = (csdY[t+1][j] - csdY[i][j])/k;

    sum_EYit2 = sum_EYit2 + EYit[j] * EYit[j];
    sum_dif_Yit2 = sum_dif_Yit2 + (csdY2[t+1][j] - csdY2[i][j]);
  }
  kVYit = sum_dif_Yit2 - k * sum_EYit2;
}

double pCost::get_min() { return (kVYit + mi1beta); }
//----------------------------------------------------------------------------//
