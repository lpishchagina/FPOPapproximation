#include "pSphere.h"
#include "pRectangle.h"
#include "Geom_Iempty_Eall_4.h"

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
//#include <list>
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;

//constructor copy**************************************************************
Geom_Iempty_Eall_4::Geom_Iempty_Eall_4(const Geom_Iempty_Eall_4 & geom2){
  p = geom2.p;
  label_t = geom2.label_t;
  rect_t = new pRectangle(p);
  disks_t_1.clear();
  disks_t_1 = geom2.disks_t_1;
}
//destructor********************************************************************
Geom_Iempty_Eall_4::~Geom_Iempty_Eall_4(){delete rect_t;}

//accessory*********************************************************************
unsigned int Geom_Iempty_Eall_4::get_p()const{return p;}

unsigned int Geom_Iempty_Eall_4::get_label_t()const{return label_t;}

std::list<pSphere> Geom_Iempty_Eall_4::get_disks_t_1()const{return disks_t_1;}

//CleanGeometry*****************************************************************
void Geom_Iempty_Eall_4::CleanGeometry(){disks_t_1.clear();}

//EmptyGeometry*****************************************************************
bool Geom_Iempty_Eall_4::EmptyGeometry(){return rect_t->IsEmpty_rect();}

//InitialGeometry***************************************************************
void Geom_Iempty_Eall_4::InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks){
  label_t = t;
  disks_t_1.clear();
  disks_t_1 = disks;
}

//UpdateGeometry****************************************************************
void Geom_Iempty_Eall_4::UpdateGeometry(const pSphere &disk_t){
  //Intersection
  rect_t->Intersection_disk(disk_t);
  // Exclusions
  std::list<pSphere>::iterator iter = disks_t_1.begin();
  while(iter != disks_t_1.end() && (!rect_t->IsEmpty_rect())){
    if (rect_t->EmptyIntersection(*iter)) {iter = disks_t_1.erase(iter);}//isn't intersection => Remove disks
    else {
      rect_t->Exclusion_disk(*iter);
      ++iter;
    }
  }
}


