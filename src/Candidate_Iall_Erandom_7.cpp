#include "pRectangle.h"
#include "Candidate_Iall_Erandom_7.h"
#include "Cost.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;


Candidate_Iall_Erandom_7::Candidate_Iall_Erandom_7(const Candidate_Iall_Erandom_7 & candidate){
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;//ATTENTION!!!parce que c'est pointeur sur CumSum et vecteur dans FPOP(c'est le lien!!!)
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iall_Erandom_7::~Candidate_Iall_Erandom_7(){delete Rect;  CumSumData = NULL;  VectOfCosts = NULL;}

unsigned int Candidate_Iall_Erandom_7::GetTau()const{return Tau;}

int Candidate_Iall_Erandom_7::get_Number(int N){
  system("sleep 0.1");
  srand(time(NULL));
  int res = rand()% N + 1;
  return res;
}

void Candidate_Iall_Erandom_7::CleanOfCandidate(){CumSumData = NULL;  VectOfCosts = NULL;}

bool Candidate_Iall_Erandom_7::EmptyOfCandidate() {return Rect -> IsEmpty_rect();}

void Candidate_Iall_Erandom_7::InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts){
  Tau = t;
  CumSumData = cumsumdata;//ATTENTION!!!parce que c'est pointeur sur CumSum et vecteur dans FPOP(c'est le lien!!!)
  VectOfCosts = vectofcosts;
}

void Candidate_Iall_Erandom_7::UpdateOfCandidate(unsigned int t, std::vector<std::list<Candidate_Iall_Erandom_7>::reverse_iterator> &vectlinktocands){
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, t, CumSumData[Tau], CumSumData[t + 1], VectOfCosts[Tau]);
  double r2 = (VectOfCosts[t + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (r2 < 0){ Rect->DoEmpty_rect();}
  else{
    pSphere Disk = pSphere(Dim);
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
    Rect -> Intersection_disk(Disk);

    unsigned int N = vectlinktocands.size()-1;
    if(N >0){
      unsigned int i = get_Number(N)-1;
      unsigned int u = vectlinktocands[i] -> GetTau();
      cost.InitialCost(Dim, u, t-1, CumSumData[u], CumSumData[t], VectOfCosts[u]);
      r2 = (VectOfCosts[t] - VectOfCosts[u] - cost.get_coef_Var())/cost.get_coef();
      Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
      Rect->Exclusion_disk(Disk);
    }
  }
}
