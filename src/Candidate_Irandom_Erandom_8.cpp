#include "pRectangle.h"
#include "Candidate_Irandom_Erandom_8.h"
#include "Cost.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;


Candidate_Irandom_Erandom_8::Candidate_Irandom_Erandom_8(const Candidate_Irandom_Erandom_8 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Irandom_Erandom_8::~Candidate_Irandom_Erandom_8() { delete Rect;  CumSumData = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Irandom_Erandom_8::GetTau()const { return Tau; }

int Candidate_Irandom_Erandom_8::get_Number(int N) {
  system("sleep 0.1");
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

void Candidate_Irandom_Erandom_8::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_Irandom_Erandom_8::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Irandom_Erandom_8::InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts) {
  Tau = t;
  CumSumData = cumsumdata;
  VectOfCosts = vectofcosts;
}

void Candidate_Irandom_Erandom_8::UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Irandom_Erandom_8>::iterator> &vectlinktocands){
  unsigned int N = vectlinktocands.size();
  unsigned int t = vectlinktocands[N-1] -> GetTau();
  unsigned int u;
  double r2;
  unsigned int rand_i;
  Rect -> Clean_rect();
  Cost cost = Cost(Dim);
  pSphere Disk = pSphere(Dim);

  u = vectlinktocands[N-1] -> GetTau();
  cost.InitialCost(Dim, u, t, CumSumData[u], CumSumData[t + 1], VectOfCosts[u]);
  r2 = (VectOfCosts[t + 1] - VectOfCosts[u] - cost.get_coef_Var())/cost.get_coef();
  if (r2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  } else {
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
    Rect -> Intersection_disk(Disk);
    if (Rect -> IsEmpty_rect()) {
      return;
    }
  }
  if ((N > 1) && (i != (N - 1))) {
    rand_i = get_Number(N - 1 - i) - 1;
    u = vectlinktocands[i + rand_i] -> GetTau();
    cost.InitialCost(Dim, u, t, CumSumData[u], CumSumData[t + 1], VectOfCosts[u]);
    r2 = (VectOfCosts[t + 1] - VectOfCosts[u] - cost.get_coef_Var())/cost.get_coef();
    if (r2 < 0) {
      Rect -> DoEmpty_rect();
      return;
    } else {
      Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
      Rect -> Intersection_disk(Disk);
      if (Rect -> IsEmpty_rect()) {
        return;
      }
    }
  }
  //random exclusion
  if (i > 0) {
    rand_i = get_Number(i)-1;
    u = vectlinktocands[rand_i] -> GetTau();
    cost.InitialCost(Dim, u, t-1, CumSumData[u], CumSumData[t], VectOfCosts[u]);
    r2 = (VectOfCosts[t] - VectOfCosts[u] - cost.get_coef_Var()) / cost.get_coef();
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
    if (!Rect -> EmptyIntersection(Disk)) {
      Rect -> Exclusion_disk(Disk);
    }
  }
}
