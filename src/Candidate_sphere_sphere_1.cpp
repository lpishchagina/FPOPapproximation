#include "pSphere.h"
#include "Candidate_sphere_sphere_1.h"
#include "Cost.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;

Candidate_sphere_sphere_1::Candidate_sphere_sphere_1(const Candidate_sphere_sphere_1 & candidate){
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;//ATTENTION!!!parce que c'est pointeur sur CumSum et vecteur dans FPOP(c'est le lien!!!)
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_sphere_sphere_1::~Candidate_sphere_sphere_1(){CumSumData = NULL;  VectOfCosts = NULL;}

unsigned int Candidate_sphere_sphere_1::GetTau()const{return Tau;}

double Candidate_sphere_sphere_1::Dist(double* a, double*b){
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++){dist = dist + (a[k] - b[k])*(a[k] - b[k]);}
  return sqrt(dist);
}


void Candidate_sphere_sphere_1::CleanOfCandidate(){CumSumData = NULL;  VectOfCosts = NULL;}

bool Candidate_sphere_sphere_1::EmptyOfCandidate() {return fl_empty;}

void Candidate_sphere_sphere_1::InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts){
  Tau = t;
  CumSumData = cumsumdata;//ATTENTION!!!parce que c'est pointeur sur CumSum et vecteur dans FPOP(c'est le lien!!!)
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Candidate_sphere_sphere_1::UpdateOfCandidate(unsigned int t, std::vector<std::list<Candidate_sphere_sphere_1>::reverse_iterator> &vectlinktocands){
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, t, CumSumData[Tau], CumSumData[t + 1], VectOfCosts[Tau]);
  double r2 = (VectOfCosts[t + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (r2 < 0){
    fl_empty = true;
    return;
  }
  else{
    pSphere Disk_Taut = pSphere(Dim);
    pSphere Disk_ut1 = pSphere(Dim);
    Disk_Taut.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
    unsigned int i =0;
    double dist;
    while (vectlinktocands[i] -> GetTau() != t){
      unsigned int u = vectlinktocands[i] -> GetTau();
      cost.InitialCost(Dim, u, t-1, CumSumData[u], CumSumData[t], VectOfCosts[u]);
      r2 = (VectOfCosts[t] - VectOfCosts[u] - cost.get_coef_Var())/cost.get_coef();
      Disk_ut1.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
      dist = Dist(Disk_Taut.get_center(),Disk_ut1.get_center());

      if (dist < (Disk_ut1.get_radius() + Disk_Taut.get_radius())){
        if (dist <= (Disk_ut1.get_radius() - Disk_Taut.get_radius())){ //Exclusion is empty
          fl_empty = true;
          return;
        }
      }
      else { i++;}//Remove disks ou step???
    }
  }
}
