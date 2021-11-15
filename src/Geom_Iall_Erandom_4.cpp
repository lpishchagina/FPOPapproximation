#include "pSphere.h"
#include "pRectangle.h"
#include "Geom_Iall_Erandom_4.h"


#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
//#include <list>
#include <Rcpp.h>

#include <cstdlib>
#include <random>
#include <ctime>


using namespace Rcpp;
using namespace std;

//constructor copy**************************************************************
Geom_Iall_Erandom_4::Geom_Iall_Erandom_4(const Geom_Iall_Erandom_4 & geom2){
  p = geom2.p;
  label_t = geom2.label_t;
  rect_t = new pRectangle(p);
  disks_t_1.clear();
  disks_t_1 = geom2.disks_t_1;
}
//destructor********************************************************************
Geom_Iall_Erandom_4::~Geom_Iall_Erandom_4(){delete rect_t;}

//accessory*********************************************************************
unsigned int Geom_Iall_Erandom_4::get_p()const{return p;}

unsigned int Geom_Iall_Erandom_4::get_label_t()const{return label_t;}

std::list<pSphere> Geom_Iall_Erandom_4::get_disks_t_1()const{return disks_t_1;}

//CleanGeometry*****************************************************************
void Geom_Iall_Erandom_4::CleanGeometry(){disks_t_1.clear();}

//EmptyGeometry*****************************************************************
bool Geom_Iall_Erandom_4::EmptyGeometry(){return rect_t->IsEmpty_rect();}

//InitialGeometry***************************************************************
void Geom_Iall_Erandom_4::InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks){
  label_t = t;
  disks_t_1.clear();
  disks_t_1 = disks;
}

//UpdateGeometry****************************************************************
void Geom_Iall_Erandom_4::UpdateGeometry(const pSphere &disk_t){
  //Intersection
  rect_t->Intersection_disk(disk_t);
  // Exclusion random
  if (!rect_t->IsEmpty_rect()){
    int nb_disks = disks_t_1.size();
    if (nb_disks > 0){
      int nb_rand = get_Number(nb_disks);
      std::list<pSphere>::iterator iter = disks_t_1.begin();
      for (int i = 0; i < (nb_rand-2); i ++){  ++iter;}
      rect_t->Exclusion_disk(*iter);
    }
  }
}
//*****************************************************************************
int Geom_Iall_Erandom_4::get_Number(int N){
  system("sleep 0.1");
  srand(time(NULL));
  int res = rand()% N + 1;
  return res;
}
