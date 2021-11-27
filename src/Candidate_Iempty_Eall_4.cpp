#include "pRectangle.h"
#include "Candidate_Iempty_Eall_4.h"
#include "Cost.h"
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;

Candidate_Iempty_Eall_4::Candidate_Iempty_Eall_4(const Candidate_Iempty_Eall_4 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iempty_Eall_4::~Candidate_Iempty_Eall_4() { delete Rect;  CumSumData = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Iempty_Eall_4::GetTau()const { return Tau; }

void Candidate_Iempty_Eall_4::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_Iempty_Eall_4::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Iempty_Eall_4::InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts) {
  Tau = t;
  CumSumData = cumsumdata;
  VectOfCosts = vectofcosts;
}

void Candidate_Iempty_Eall_4::UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Iempty_Eall_4>::iterator> &vectlinktocands) {
  unsigned int N = vectlinktocands.size();
  unsigned int t = vectlinktocands[N-1] -> GetTau();
  unsigned int u;
  double r2;
  Rect -> Clean_rect();
  Cost cost = Cost(Dim);
  pSphere Disk = pSphere(Dim);
  //intersection
  cost.InitialCost(Dim, t, t, CumSumData[t], CumSumData[t + 1], VectOfCosts[t]);
  r2 = (VectOfCosts[t + 1] - VectOfCosts[t] - cost.get_coef_Var())/cost.get_coef();
  if (r2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  } else {
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
    Rect -> Intersection_disk(Disk);
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
          Rcpp::Rcout<<"FPOP"<<endl;
          return;
        }
      }
    }
  }
}
