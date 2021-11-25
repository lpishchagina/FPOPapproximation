#include "pSphere.h"
#include "pRectangle.h"
#include "Geom_Irandom_Erandom_8.h"


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
Geom_Irandom_Erandom_8::Geom_Irandom_Erandom_8(const Geom_Irandom_Erandom_8 & geom2){
  p = geom2.p;
  label_t = geom2.label_t;
  rect_t = new pRectangle(p);
  disks_t_1.clear();
  disks_t_1 = geom2.disks_t_1;
  disks_tu.clear();
  disks_tu = geom2.disks_tu;
}
//destructor********************************************************************
Geom_Irandom_Erandom_8::~Geom_Irandom_Erandom_8(){delete rect_t;}

//accessory*********************************************************************
unsigned int Geom_Irandom_Erandom_8::get_p()const{return p;}

unsigned int Geom_Irandom_Erandom_8::get_label_t()const{return label_t;}

std::list<pSphere> Geom_Irandom_Erandom_8::get_disks_t_1()const{return disks_t_1;}

std::list<pSphere> Geom_Irandom_Erandom_8::get_disks_tu()const{return disks_tu;}

//CleanGeometry*****************************************************************
void Geom_Irandom_Erandom_8::CleanGeometry(){
  disks_t_1.clear();
  disks_tu.clear();
}

//EmptyGeometry*****************************************************************
bool Geom_Irandom_Erandom_8::EmptyGeometry(){return rect_t->IsEmpty_rect();}

//InitialGeometry***************************************************************
void Geom_Irandom_Erandom_8::InitialGeometry(unsigned int dim, unsigned int t, const std::list<pSphere> &disks){
  label_t = t;
  disks_t_1.clear();
  disks_t_1 = disks;
}

//UpdateGeometry****************************************************************
void Geom_Irandom_Erandom_8::UpdateGeometry(const pSphere &disk_t){
  //Intersection random
  pRectangle* pcube = new pRectangle(p);
  std::list<pSphere>::iterator iter = disks_tu.begin();
  int nb_disks = disks_tu.size();
  int nb_rand;
  if (nb_disks != 0){
    nb_rand = get_Number(nb_disks);
    for (int i = 0; i < (nb_rand-1); i ++){  ++iter;}
    pcube->Intersection_disk(*iter);
  }
  pcube->Intersection_disk(disk_t);
  disks_tu.push_back(disk_t);
  rect_t = pcube;
  // Exclusion random
  if (!rect_t->IsEmpty_rect()){
    //Rcpp::Rcout<<"nb disks"<<nb_disks<<endl;
    nb_disks = disks_t_1.size();
    if (nb_disks > 0){
      nb_rand = get_Number(nb_disks);
      iter = disks_t_1.begin();
      for (int i = 0; i < (nb_rand-1); i ++){  ++iter;}
      if (rect_t->EmptyIntersection(*iter)) {
        iter = disks_t_1.erase(iter);
      }//isn't intersection => Remove disks
      else {rect_t->Exclusion_disk(*iter);}
    }
  }

}
//*****************************************************************************
int Geom_Irandom_Erandom_8::get_Number(int N){
  system("sleep 0.1");
  srand(time(NULL));
  int res = rand()% N + 1;
  return res;
}
