#include "pRectangle.h"
#include "Candidate_Iall_Eall_2.h"
#include "Cost.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;


Candidate_Iall_Eall_2::Candidate_Iall_Eall_2(const Candidate_Iall_Eall_2 & candidate){
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;//ATTENTION!!!parce que c'est pointeur sur CumSum et vecteur dans FPOP(c'est le lien!!!)
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iall_Eall_2::~Candidate_Iall_Eall_2(){delete Rect;  CumSumData = NULL;  VectOfCosts = NULL;}

unsigned int Candidate_Iall_Eall_2::GetTau()const{return Tau;}

void Candidate_Iall_Eall_2::CleanOfCandidate(){CumSumData = NULL;  VectOfCosts = NULL;}

bool Candidate_Iall_Eall_2::EmptyOfCandidate() {return Rect -> IsEmpty_rect();}

void Candidate_Iall_Eall_2::InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts){
  Tau = t;
  CumSumData = cumsumdata;//ATTENTION!!!parce que c'est pointeur sur CumSum et vecteur dans FPOP(c'est le lien!!!)
  VectOfCosts = vectofcosts;
}

void Candidate_Iall_Eall_2::UpdateOfCandidate(unsigned int t, std::vector<std::list<Candidate_Iall_Eall_2>::reverse_iterator> &vectlinktocands){
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, t, CumSumData[Tau], CumSumData[t + 1], VectOfCosts[Tau]);
  double r2 = (VectOfCosts[t + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (r2 < 0){ Rect->DoEmpty_rect();}
  else{
    pSphere Disk = pSphere(Dim);
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
    Rect -> Intersection_disk(Disk);

    unsigned int i =0;
    while ( (vectlinktocands[i] -> GetTau() != t) && (!Rect->IsEmpty_rect()) ){
      unsigned int u = vectlinktocands[i] -> GetTau();
      cost.InitialCost(Dim, u, t-1, CumSumData[u], CumSumData[t], VectOfCosts[u]);
      r2 = (VectOfCosts[t] - VectOfCosts[u] - cost.get_coef_Var())/cost.get_coef();
      Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
      if (Rect->EmptyIntersection(Disk)) {i++;}
      else {
        Rect->Exclusion_disk(Disk);
        i++;
      }
    }
  }
}

