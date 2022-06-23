#include "pRectangle.h"
#include "pSphere.h"
#include <Rcpp.h>

using namespace Rcpp;

using namespace std;
//constructor, copy, destructor-------------------------------------------------
pRectangle::pRectangle(unsigned int dim) {
  p = dim;
  borders = new double*[p];
  for (unsigned int i = 0; i < p; i++) {
    borders[i] = new double[2];
    borders[i][0] = -INFINITY;
    borders[i][1] = INFINITY;
  }
}

pRectangle::pRectangle(unsigned int dim, double** coords) {
  p = dim;
  borders = new double*[p];
  for (unsigned int i = 0; i < p; i++) {
    borders[i] = new double[2];
    borders[i][0] = coords[i][0];
    borders[i][1] = coords[i][1];
  }
}

pRectangle::pRectangle(const pRectangle &rect) {
  p = rect.p;
  borders = new double*[p];
  for(unsigned int i = 0; i < p; i++) {
    borders[i] = new double[2];
    borders[i][0] = rect.borders[i][0];
    borders[i][1] = rect.borders[i][1];
  }
}

pRectangle::~pRectangle() {
  for(unsigned int i = 0; i < p; i++) {
    delete[]borders[i];
  }
  delete[]borders;
  borders = NULL;
}

//accessory---------------------------------------------------------------------
double** pRectangle::get_borders() const { return borders; }
unsigned int pRectangle::get_p() const { return p; }

//tools-------------------------------------------------------------------------
double pRectangle::min_ab(double a, double b) { if (a <= b) { return a; } else { return b; } }
double pRectangle::max_ab(double a, double b){ if (a >= b) { return a; } else { return b; } }

double pRectangle::get_dist(double* pnt1, double* pnt2) {
  double res = 0;
  for (unsigned int k = 0; k < p; k++) {
    res = res + (pnt2[k] - pnt1[k]) * (pnt2[k] - pnt1[k]);
  }
  return sqrt(res);
}

void pRectangle::DoEmptyRect() { borders[1][0] = borders[1][1]; }

bool pRectangle::IsEmptyRect() {
  for (unsigned int k = 0; k < p; k++) {
    if (borders[k][0] >= borders[k][1]) {
      return true;
    }
  }
  return false;
}

bool pRectangle::EmptyIntersection(const pSphere &disk) {
  double* c = disk.get_c();
  double* pnt_min = new double[p];//closest point
  for (unsigned int k = 0; k < p; k++){
    pnt_min[k] = c[k];
    if (c[k] <= borders[k][0]) { pnt_min[k] = borders[k][0]; }
    if (c[k] >= borders[k][1]) { pnt_min[k] = borders[k][1]; }
  }
  if (get_dist(pnt_min, c) >= disk.get_r()) {
    delete [] pnt_min;
    pnt_min = NULL;
    return true;
    } else {
    delete [] pnt_min;
    pnt_min = NULL;
    return false;
  }
}

void pRectangle::SphereApproximation(const pSphere &disk) {
  double* c = disk.get_c();
  for (unsigned int k = 0; k < p; k++) {
    borders[k][0] = c[k] - disk.get_r();
    borders[k][1] = c[k] + disk.get_r();
  }
}

void pRectangle::IntersectionSphere(const pSphere &disk) {
  double* c = disk.get_c();
  double* pnt_min = new double[p];//closest point
  for (unsigned int k = 0; k < p; k++) {
    pnt_min[k] = c[k];
    if (c[k] <= borders[k][0]){ pnt_min[k] = borders[k][0];}
    if (c[k] >= borders[k][1]){ pnt_min[k] = borders[k][1];}
  }
  double* dx2 = new double[p];//discriminant
  double dx2i = 1;
  unsigned int i = 0;
  while ((dx2i > 0) && (i < p)) {
    dx2i = 0;
    for (unsigned int j = 0; j < p; j++) {
      if (j != i) {
        dx2i = dx2i + (pnt_min[j] - c[j]) * (pnt_min[j] - c[j]);
      }
    }
    dx2i = disk.get_r() * disk.get_r() - dx2i;
    dx2[i] = dx2i;
    ++i;
  }
  if (i != p) {
    borders[0][0] =  borders[0][1];
  } else {
    for (unsigned int k = 0; k < p; k++) {
      borders[k][0] = max_ab(borders[k][0], c[k] - sqrt(dx2[k]));
      borders[k][1] = min_ab(borders[k][1], c[k] + sqrt(dx2[k]));
    }
  }
  delete [] pnt_min;
  delete [] dx2;
  pnt_min = NULL;
  dx2 = NULL;
}

void pRectangle::ExclusionSphere(const pSphere &disk) {
  double* c = disk.get_c();
  double dx2;
  double* pnt_max = new double[p];//farthest point
  for (unsigned int k = 0; k < p; k++){
    if (abs(c[k] - borders[k][1]) >= abs(c[k] - borders[k][0])) { pnt_max[k] = borders[k][1]; }
    else { pnt_max[k] = borders[k][0]; }
  }
  for (unsigned int k = 0; k < p; k++) {//discriminant
    dx2 = 0;
    for (unsigned int j = 0; j < p; j++) {
      if (j != k) {
        dx2 = dx2 + (pnt_max[j] - c[j]) * (pnt_max[j] - c[j]);
      }
    }
    dx2 = disk.get_r() * disk.get_r() - dx2;
    if (dx2 > 0) {
      if ((pnt_max[k] == borders[k][0]) && (borders[k][1] <=  c[k] + sqrt(dx2))) {
        borders[k][1] = min_ab(borders[k][1], c[k] - sqrt(dx2));
      }
      if ((pnt_max[k] == borders[k][1]) && (borders[k][0] >=  c[k] - sqrt(dx2))) {
        borders[k][0] = max_ab(borders[k][0], c[k] + sqrt(dx2));
      }
    }
  }
  delete [] pnt_max;
  pnt_max = NULL;
}

//----------------------------------------------------------------------------//
