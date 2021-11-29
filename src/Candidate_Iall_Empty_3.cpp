#include "Candidate_Iall_Eempty_3.h"

using namespace Rcpp;
using namespace std;

Candidate_Iall_Eempty_3::Candidate_Iall_Eempty_3(const Candidate_Iall_Eempty_3 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iall_Eempty_3::~Candidate_Iall_Eempty_3(){ delete Rect;  CumSumData = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Iall_Eempty_3::GetTau()const { return Tau; }

void Candidate_Iall_Eempty_3::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_Iall_Eempty_3::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Iall_Eempty_3::InitialOfCandidate(unsigned int t, double** &cumsumdata, double* &vectofcosts) {
  Tau = t;
  CumSumData = cumsumdata;
  VectOfCosts = vectofcosts;
}

void Candidate_Iall_Eempty_3::UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Iall_Eempty_3>::iterator> &vectlinktocands) {
  unsigned int N = vectlinktocands.size();
  unsigned int t = vectlinktocands[N-1] -> GetTau();
  unsigned int u;
  double r2;
  Rect -> Clean_rect();
  Cost cost = Cost(Dim);
  pSphere Disk = pSphere(Dim);
  //$\R_tau^t = cap_{j=u}^t S_u^t $
  for (unsigned int j = i; j < N; j++) {
    u = vectlinktocands[j] -> GetTau();
    cost.InitialCost(Dim, u, t, CumSumData[u], CumSumData[t + 1], VectOfCosts[u]);
    r2 = (VectOfCosts[t + 1] - VectOfCosts[u] - cost.get_coef_Var()) / cost.get_coef();
    if (r2 < 0) {
      Rect -> DoEmpty_rect();
      return;
    }
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(r2));
    Rect -> Intersection_disk(Disk);
    if (Rect -> IsEmpty_rect()) {
      Rcpp::Rcout<< "FPOP"<<endl;
      return;
    }
  }
}

