#include "pSphere.h"
#include "pRectangle.h"
<<<<<<< HEAD:src/Geom_Ilast1_Erandom_6.cpp
#include "Geom_Ilast1_Erandom_6.h"
=======
#include "Geom_Iall_Erandom_7.h"

>>>>>>> 1197feb576c482346f0434fbeea9d18b0b9a9950:src/Geom_Iall_Erandom_7.cpp

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
<<<<<<< HEAD:src/Geom_Ilast1_Erandom_6.cpp
Geom_Ilast1_Erandom_6::Geom_Ilast1_Erandom_6(const Geom_Ilast1_Erandom_6 & geom2){
=======
Geom_Iall_Erandom_7::Geom_Iall_Erandom_7(const Geom_Iall_Erandom_7 & geom2){
>>>>>>> 1197feb576c482346f0434fbeea9d18b0b9a9950:src/Geom_Iall_Erandom_7.cpp
  p = geom2.p;
  label_t = geom2.label_t;
  rect_t = new pRectangle(p);
  disks_t_1.clear();
  disks_t_1 = geom2.disks_t_1;
}
//destructor********************************************************************
<<<<<<< HEAD:src/Geom_Ilast1_Erandom_6.cpp
Geom_Ilast1_Erandom_6::~Geom_Ilast1_Erandom_6(){delete rect_t;}

//accessory*********************************************************************
unsigned int Geom_Ilast1_Erandom_6::get_p()const{return p;}

unsigned int Geom_Ilast1_Erandom_6::get_label_t()const{return label_t;}

std::list<pSphere> Geom_Ilast1_Erandom_6::get_disks_t_1()const{return disks_t_1;}

//CleanGeometry*****************************************************************
void Geom_Ilast1_Erandom_6::CleanGeometry(){disks_t_1.clear();}

//EmptyGeometry*****************************************************************
bool Geom_Ilast1_Erandom_6::EmptyGeometry(){return rect_t->IsEmpty_rect();}

//InitialGeometry***************************************************************
void Geom_Ilast1_Erandom_6::InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks){
=======
Geom_Iall_Erandom_7::~Geom_Iall_Erandom_7(){delete rect_t;}

//accessory*********************************************************************
unsigned int Geom_Iall_Erandom_7::get_p()const{return p;}

unsigned int Geom_Iall_Erandom_7::get_label_t()const{return label_t;}

std::list<pSphere> Geom_Iall_Erandom_7::get_disks_t_1()const{return disks_t_1;}

//CleanGeometry*****************************************************************
void Geom_Iall_Erandom_7::CleanGeometry(){disks_t_1.clear();}

//EmptyGeometry*****************************************************************
bool Geom_Iall_Erandom_7::EmptyGeometry(){return rect_t->IsEmpty_rect();}

//InitialGeometry***************************************************************
void Geom_Iall_Erandom_7::InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks){
>>>>>>> 1197feb576c482346f0434fbeea9d18b0b9a9950:src/Geom_Iall_Erandom_7.cpp
  label_t = t;
  disks_t_1.clear();
  disks_t_1 = disks;
}

//UpdateGeometry****************************************************************
<<<<<<< HEAD:src/Geom_Ilast1_Erandom_6.cpp
void Geom_Ilast1_Erandom_6::UpdateGeometry(const pSphere &disk_t){
  //last1
  pRectangle* pcube = new pRectangle(p);
  pcube->Intersection_disk(disk_t);
  rect_t = pcube;
  // Exclusion random
  if (!rect_t->IsEmpty_rect()){
    int nb_disks = disks_t_1.size();
    //Rcpp::Rcout<<"nb disks"<<nb_disks<<endl;
=======
void Geom_Iall_Erandom_7::UpdateGeometry(const pSphere &disk_t){
  //Intersection
  rect_t->Intersection_disk(disk_t);
  // Exclusion random
  if (!rect_t->IsEmpty_rect()){
    int nb_disks = disks_t_1.size();
//    Rcpp::Rcout<<"nb disks"<<nb_disks<<endl;
>>>>>>> 1197feb576c482346f0434fbeea9d18b0b9a9950:src/Geom_Iall_Erandom_7.cpp
    if (nb_disks > 0){
      int nb_rand = get_Number(nb_disks);
      std::list<pSphere>::iterator iter = disks_t_1.begin();
      for (int i = 0; i < (nb_rand-1); i ++){  ++iter;}
      if (rect_t->EmptyIntersection(*iter)) {
        iter = disks_t_1.erase(iter);
      }//isn't intersection => Remove disks
      else {rect_t->Exclusion_disk(*iter);}
    }
  }
}
<<<<<<< HEAD:src/Geom_Ilast1_Erandom_6.cpp
//******************************************************************************
int Geom_Ilast1_Erandom_6::get_Number(int N){
=======
//*****************************************************************************
int Geom_Iall_Erandom_7::get_Number(int N){
>>>>>>> 1197feb576c482346f0434fbeea9d18b0b9a9950:src/Geom_Iall_Erandom_7.cpp
  system("sleep 0.1");
  srand(time(NULL));
  int res = rand()% N + 1;
  return res;
}

