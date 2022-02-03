#include "pRectangle.h"
#include "pSphere.h"
//#include "math.h"
//#include<iostream>
#include <Rcpp.h>

using namespace Rcpp;

using namespace std;
//constructor*******************************************************************
pRectangle::pRectangle(unsigned int dim){
  p = dim;
  coordinates = new double*[p];
  for (unsigned int i = 0; i < p; i++) {
    coordinates[i] = new double[2];
    coordinates[i][0] = -INFINITY;
    coordinates[i][1] = INFINITY;
  }
}

pRectangle::pRectangle(unsigned int dim, double** coords){
  p = dim;
  coordinates = new double*[p];
  for (unsigned int i = 0; i < p; i++) {
    coordinates[i] = new double[2];
    coordinates[i][0] = coords[i][0];
    coordinates[i][1] = coords[i][1];
  }
}
//constructor copy**************************************************************
pRectangle::pRectangle(const pRectangle &rect){
  p = rect.p;
  coordinates = new double*[p];
  for(unsigned int i = 0; i < p; i++) {
    coordinates[i] = new double[2];
    coordinates[i][0] = rect.coordinates[i][0];
    coordinates[i][1] = rect.coordinates[i][1];
  }
}
//constructor, copy and destructor**********************************************
pRectangle::~pRectangle(){
  for(unsigned int i = 0; i < p; i++) {delete[]coordinates[i];}
  delete[]coordinates;
  coordinates = NULL;
}

//accessory*********************************************************************
double** pRectangle::get_coordinates()const{return coordinates;}
unsigned int pRectangle::get_p()const{return p;}

double pRectangle::min_ab(double a, double b){if (a < b) {return a;} else {return b;}}
double pRectangle::max_ab(double a, double b){if (a > b) {return a;} else {return b;}}

//IsEmpty_rect******************************************************************
bool pRectangle::IsEmpty_rect(){
  for (unsigned int k = 0; k < p; k++){if (coordinates[k][0] >= coordinates[k][1]){return true;}}
  return false;
}
void pRectangle::DoEmpty_rect(){coordinates[1][0] = coordinates[1][1];}

//EmptyIntersection*************************************************************
bool pRectangle::EmptyIntersection(const pSphere &disk){
  double* c = disk.get_center();
  //point_min-------------------------------------------------------------------
  double* pnt_min = new double[p];
  for (unsigned int k = 0; k < p; k++){
    pnt_min[k] = c[k];
    if (c[k] <= coordinates[k][0]){ pnt_min[k] = coordinates[k][0];}
    if (c[k] >= coordinates[k][1]){ pnt_min[k] = coordinates[k][1];}
  }
  //distance--------------------------------------------------------------------
  double  d = 0;
  for (unsigned int k = 0; k < p; k++){ d = d + (pnt_min[k] - c[k]) * (pnt_min[k] - c[k]);}
  //memory----------------------------------------------------------------------
  delete [] pnt_min;
  pnt_min = NULL;
  //check-----------------------------------------------------------------------
  if (sqrt(d) >= disk.get_radius()) {return true;}
  else {return false;}
}


//clean_rect------------------------------------------------------------------
void pRectangle::Clean_rect(){
  for (unsigned int i = 0; i < p; i++) {
    coordinates[i][0] = -INFINITY;
    coordinates[i][1] = INFINITY;
  }
}
//Intersection_disk*************************************************************
void pRectangle::Intersection_disk(const pSphere &disk){
  double r = disk.get_radius();
  double* c = disk.get_center();
  //point_min-------------------------------------------------------------------
  double* pnt_min = new double[p];
  for (unsigned int k = 0; k < p; k++){
    pnt_min[k] = c[k];
    if (c[k] <= coordinates[k][0]){ pnt_min[k] = coordinates[k][0];}
    if (c[k] >= coordinates[k][1]){ pnt_min[k] = coordinates[k][1];}
  }
  //discriminant----------------------------------------------------------------
  double* dx2 = new double[p];
  double dx2i = 1;
  unsigned int i = 0;
  while ((dx2i > 0)&&(i < p)){
    dx2i = 0;
    for (unsigned int j = 0; j < p; j++){if (j != i){dx2i = dx2i + (pnt_min[j] - c[j]) * (pnt_min[j] - c[j]);}}
    dx2i = r * r - dx2i;
    dx2[i] = dx2i;
    ++i;
  }
  //----------------------------------------------------------------------------
  if (i != p) {coordinates[0][0] =  coordinates[0][1];}
  else{
    for (unsigned int k = 0; k < p; k++){
      coordinates[k][0] = max_ab(coordinates[k][0], c[k] - sqrt(dx2[k]));
      coordinates[k][1] = min_ab(coordinates[k][1], c[k] + sqrt(dx2[k]));
    }
  }
  //memory----------------------------------------------------------------------
  delete [] pnt_min;
  delete [] dx2;
  pnt_min = NULL;
  dx2 = NULL;
}

//Exclusion_disk****************************************************************
void pRectangle::Exclusion_disk(const pSphere &disk){
  double r = disk.get_radius();
  double* c = disk.get_center();
  double dx2;
  //-point_max------------------------------------------------------------------
  double* pnt_max = new double[p];
  for (unsigned int k = 0; k < p; k++){
    if (abs(c[k] - coordinates[k][1]) >= abs(c[k] - coordinates[k][0])) {pnt_max[k] = coordinates[k][1];}
    else{pnt_max[k] = coordinates[k][0];}
  }
  //discriminant----------------------------------------------------------------
  for (unsigned int k = 0; k < p; k++){
    dx2 = 0;
    for (unsigned int j = 0; j < p; j++){if (j != k){dx2 = dx2 + (pnt_max[j] - c[j]) * (pnt_max[j] - c[j]);}}
    dx2 = r * r - dx2;
    if (dx2 > 0){
      if ((pnt_max[k] == coordinates[k][0]) && (coordinates[k][1] <=  c[k] + sqrt(dx2))) {coordinates[k][1] = min_ab(coordinates[k][1], c[k] - sqrt(dx2));}
      if ((pnt_max[k] == coordinates[k][1]) && (coordinates[k][0] >=  c[k] - sqrt(dx2))) {coordinates[k][0] = max_ab(coordinates[k][0], c[k] + sqrt(dx2));}
    }
  }
  //memory----------------------------------------------------------------------
  delete [] pnt_max;
  pnt_max = NULL;
}
//CubeApproximation*************************************************************
void pRectangle::CubeApproximation(const pSphere &disk){
  double r = disk.get_radius();
  double* c = disk.get_center();
  for (unsigned int k = 0; k < p; k++){
    coordinates[k][0] = c[k] - r;
    coordinates[k][1] = c[k] + r;
  }
}
//******************************************************************************
