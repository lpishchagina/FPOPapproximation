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


Candidate_Iall_Eall_2::Candidate_Iall_Eall_2(const Candidate_Iall_Eall_2 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;//ATTENTION!!!parce que c'est pointeur sur CumSum et vecteur dans FPOP(c'est le lien!!!)
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iall_Eall_2::~Candidate_Iall_Eall_2() { delete Rect;  CumSumData = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Iall_Eall_2::GetTau()const { return Tau; }

void Candidate_Iall_Eall_2::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_Iall_Eall_2::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Iall_Eall_2::InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts) {
  Tau = t;
  CumSumData = cumsumdata;
  VectOfCosts = vectofcosts;
}

void Candidate_Iall_Eall_2::UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Iall_Eall_2>::iterator> &vectlinktocands) {
  unsigned int N = vectlinktocands.size();
  unsigned int t = vectlinktocands[N-1] -> GetTau();
  unsigned int u;
  double r2;
  Rect -> Clean_rect();
  Cost cost = Cost(Dim);
  pSphere Disk = pSphere(Dim);
  //intersection
  for (unsigned int j = i; j < N; j++) {
    u = vectlinktocands[j] -> GetTau();
    cost.InitialCost(Dim, u, t, CumSumData[u], CumSumData[t + 1], VectOfCosts[u]);
    r2 = (VectOfCosts[t + 1] - VectOfCosts[u] - cost.get_coef_Var())/cost.get_coef();
    if (r2 < 0) {
      Rect -> DoEmpty_rect();
      return;
    }
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
    Rect -> Intersection_disk(Disk);
    if (Rect -> IsEmpty_rect()) {
      return;
    }
  }
  //exclusion
  if (N > 1) {
    for (unsigned int j = 0; j < i; j++) {
      u = vectlinktocands[j] -> GetTau();
      cost.InitialCost(Dim, u, t-1, CumSumData[u], CumSumData[t], VectOfCosts[u]);
      r2 = (VectOfCosts[t] - VectOfCosts[u] - cost.get_coef_Var()) / cost.get_coef();
      Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
      if (!Rect -> EmptyIntersection(Disk)) {
        Rect -> Exclusion_disk(Disk);
        if (Rect -> IsEmpty_rect()) {
          return;
        }
      }
    }
  }
}
