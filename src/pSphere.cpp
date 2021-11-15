#include "pSphere.h"
#include "math.h"

#include<iostream>
#include <Rcpp.h>

using namespace Rcpp;

using namespace std;
//constructor, copy and destructor**********************************************
pSphere::pSphere(unsigned int dim){
  p = dim;
  radius = 0;
  center = new double[p];
}

pSphere::pSphere(unsigned int dim, double* c, double r){
  p = dim;
  radius = r;
  center = new double[p];
  for (unsigned int i = 0; i < p; i++){center[i] = c[i];}
}

pSphere::pSphere(const pSphere &sphere){
  p = sphere.p;
  radius = sphere.radius;
  center = new double[p];
  for (unsigned int i = 0; i < p; i++){center[i] = sphere.center[i];}
}

pSphere::~pSphere(){delete [] center; center = NULL;}

//accessory*********************************************************************
unsigned int  pSphere::get_p()const{return p;}

double pSphere::get_radius() const{return radius;}

double* pSphere::get_center()const{return center;}
//InitialpSphere*****************************************************************
void pSphere::InitialpSphere(unsigned int dim, double* c, double r){
  p = dim;
  radius = r;
  for (unsigned int i = 0; i < p; i++){center[i] = c[i];}
}
//******************************************************************************


