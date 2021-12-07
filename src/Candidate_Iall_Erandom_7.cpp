#include "Candidate_Iall_Erandom_7.h"

using namespace Rcpp;
using namespace std;

Candidate_Iall_Erandom_7::Candidate_Iall_Erandom_7(const Candidate_Iall_Erandom_7 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect = new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iall_Erandom_7::~Candidate_Iall_Erandom_7() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Candidate_Iall_Erandom_7::GetTau()const { return Tau; }

int Candidate_Iall_Erandom_7::get_Number(int N) {
  system("sleep 0.1");
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

void Candidate_Iall_Erandom_7::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

bool Candidate_Iall_Erandom_7::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Iall_Erandom_7::InitialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = t;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Candidate_Iall_Erandom_7::UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Iall_Erandom_7>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  unsigned int N = vectlinktocands.size();
  unsigned int t = vectlinktocands[N-1] -> GetTau();
  unsigned int u;
  double Radius2;
  RealNbExclus = 0;
  Rect -> Clean_rect();
  Cost cost = Cost(Dim);
  pSphere Disk = pSphere(Dim);
  //intersection
  for (unsigned int j = i; j < N; j++) {
    u = vectlinktocands[j] -> GetTau();
    cost.InitialCost(Dim, u, t, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[t + 1] - VectOfCosts[u] - cost.get_coef_Var())/cost.get_coef();
    if (Radius2 < 0) {
      Rect -> DoEmpty_rect();
      return;
    }
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
    Rect -> Intersection_disk(Disk);
    if (Rect -> IsEmpty_rect()) {
      return;
    }
  }
  //random exclusion
  if (i > 0) {
    unsigned int rand_i = get_Number(i) - 1;
    u = vectlinktocands[rand_i] -> GetTau();
    cost.InitialCost(Dim, u, t-1, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[t] - VectOfCosts[u] - cost.get_coef_Var()) / cost.get_coef();
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
    if (!Rect -> EmptyIntersection(Disk)) {
      Rect -> Exclusion_disk(Disk);
      RealNbExclus = 1;
    }
  }
}
