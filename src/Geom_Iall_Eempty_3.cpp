#include "Geom_Iall_Eempty_3.h"

#include <math.h>
#include <iostream>
#include <list>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;


//constructor copy**************************************************************
Geom_Iall_Eempty_3::Geom_Iall_Eempty_3(const Geom_Iall_Eempty_3 & geom1){
  p = geom1.p;
  label_t = geom1.label_t;
  rect_t = new pRectangle(p);
}

//destructor********************************************************************
Geom_Iall_Eempty_3::~Geom_Iall_Eempty_3(){delete rect_t;}

//accessory*********************************************************************
unsigned int Geom_Iall_Eempty_3::get_p()const{return p;}

unsigned int Geom_Iall_Eempty_3::get_label_t()const{return label_t;}

std::list<pSphere> Geom_Iall_Eempty_3::get_disks_t_1()const{
  std::list<pSphere> list_NULL;
  list_NULL.clear();
  return list_NULL;
}

//EmptyGeometry*****************************************************************
bool Geom_Iall_Eempty_3::EmptyGeometry(){return rect_t->IsEmpty_rect();}

//InitialGeometry***************************************************************
void Geom_Iall_Eempty_3::InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks){label_t = t;}

//UpdateGeometry****************************************************************
void Geom_Iall_Eempty_3::UpdateGeometry(const pSphere &disk_t){rect_t->Intersection_disk(disk_t);}

//******************************************************************************

