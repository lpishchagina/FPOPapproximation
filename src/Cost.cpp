#include "Cost.h"

#include <iostream>
#include "math.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
//constructor*******************************************************************
Cost::Cost(unsigned int dim){
  p = dim;
  coef = 0;
  mi_1_p = 0;
  coef_Var = 0;
  mu = new double[p];
}

Cost::Cost(unsigned int dim, unsigned int i, unsigned int t, double* si_1, double* st, double mi_1pen){
  p = dim;
  coef = t - i + 1;
  mi_1_p = mi_1pen;
  mu = new double[p];

  double sum_mu2 = 0;
  double sum_dif_x2 = 0;

  for (unsigned int k = 0; k < p; k++){
    mu[k] = (st[k] - si_1[k])/coef;

    sum_mu2 = sum_mu2 + mu[k] * mu[k];
    sum_dif_x2 = sum_dif_x2 + (st[p + k] - si_1[p + k]);
  }
  coef_Var = sum_dif_x2 - coef * sum_mu2;
}
//constructor copy***************r**********************************************
Cost::Cost(const Cost &cost){
  p = cost.p;
  coef = cost.coef;
  mi_1_p = cost.mi_1_p;
  coef_Var =cost.coef_Var;
  mu = new double[p];
  for (unsigned int k = 0; k < p; k++){mu[k] = cost.mu[k];}
}
//destructor********************************************************************
Cost::~Cost(){delete [] mu; mu = NULL;}

//InitialCost***********************************************************
void Cost::InitialCost(unsigned int dim, unsigned int i, unsigned int t, double* &si_1, double* &st, double mi_1pen){
  p = dim;
  coef = t - i + 1;
  mi_1_p = mi_1pen;

  double sum_mu2 = 0;
  double sum_dif_x2 = 0;

  for (unsigned int k = 0; k < p; k++){
    mu[k] = (st[k] - si_1[k])/coef;

    sum_mu2 = sum_mu2 + mu[k] * mu[k];
    sum_dif_x2 = sum_dif_x2 + (st[p + k] - si_1[p + k]);
  }
  coef_Var = sum_dif_x2 - coef * sum_mu2;
}

//accessory*********************************************************************

unsigned int Cost::get_p() const{return p;}

unsigned int Cost::get_coef() const{return coef;}

double Cost::get_coef_Var() const{return coef_Var;}

double Cost::get_mi_1_p() const{return mi_1_p;}

double* Cost::get_mu(){return mu;}

double Cost::get_min(){ return(coef_Var + mi_1_p);}
//******************************************************************************

