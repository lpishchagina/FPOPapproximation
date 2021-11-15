#include "Geom_sphere_sphere_1.h"

#include <iostream>
#include <iterator>
#include <list>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

//constructor copy**************************************************************
Geom_sphere_sphere_1::Geom_sphere_sphere_1(const Geom_sphere_sphere_1 & geom3){
  p = geom3.p;
  label_t = geom3.label_t;
  disks_t_1.clear();
  disks_t_1 = geom3.disks_t_1;
  fl_empty = geom3.fl_empty;
}

//accessory*********************************************************************
unsigned int Geom_sphere_sphere_1::get_p()const{return p;}

unsigned int Geom_sphere_sphere_1::get_label_t()const{return label_t;}

bool Geom_sphere_sphere_1::get_fl_empty() const {return fl_empty;}

std::list<pSphere> Geom_sphere_sphere_1::get_disks_t_1() const{return disks_t_1;}

//Dist**************************************************************************
double Geom_sphere_sphere_1::Dist(double* a, double*b){
  double dist = 0;
  for (unsigned int k = 0; k < p; k++){dist = dist + (a[k] - b[k])*(a[k] - b[k]);}
  return sqrt(dist);
}

//CleanGeometry*****************************************************************
void Geom_sphere_sphere_1::CleanGeometry(){disks_t_1.clear();}

//EmptyGeometry*****************************************************************
bool Geom_sphere_sphere_1::EmptyGeometry(){return fl_empty;}

//InitialGeometry***************************************************************
void Geom_sphere_sphere_1::InitialGeometry(unsigned  int dim, unsigned  int t, const std::list<pSphere> &disks){
  label_t = t;
  fl_empty = false;
  disks_t_1.clear();
  disks_t_1 = disks;
}

//UpdateGeometry****************************************************************
void Geom_sphere_sphere_1::UpdateGeometry(const pSphere &disk_t){
  double dist;
  std::list<pSphere>::iterator iter = disks_t_1.begin();
  while( iter != disks_t_1.end()){
    dist = Dist(disk_t.get_center(),(*iter).get_center());
    if (dist < ((*iter).get_radius() + disk_t.get_radius())){
      if (dist <= ((*iter).get_radius() - disk_t.get_radius())){ //Exclusion is empty
        fl_empty = true;
        return;
      }
      else {++iter;}
    }
    else { iter = disks_t_1.erase(iter);}//Remove disks
  }
}
//******************************************************************************

