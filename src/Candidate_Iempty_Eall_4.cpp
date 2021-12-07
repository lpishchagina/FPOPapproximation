#include "Candidate_Iempty_Eall_4.h"

using namespace Rcpp;
using namespace std;

Candidate_Iempty_Eall_4::Candidate_Iempty_Eall_4(const Candidate_Iempty_Eall_4 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iempty_Eall_4::~Candidate_Iempty_Eall_4() { delete Rect;  CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Candidate_Iempty_Eall_4::GetTau()const { return Tau; }

void Candidate_Iempty_Eall_4::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

bool Candidate_Iempty_Eall_4::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Iempty_Eall_4::InitialOfCandidate(unsigned int t, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = t;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Candidate_Iempty_Eall_4::UpdateOfCandidate(unsigned int i, std::vector<std::list<Candidate_Iempty_Eall_4>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  unsigned int N = vectlinktocands.size();
  unsigned int t = vectlinktocands[N-1] -> GetTau();
  unsigned int u;
  double Radius2;
  RealNbExclus = 0;
  Rect -> Clean_rect();
  Cost cost = Cost(Dim);
  pSphere Disk = pSphere(Dim);
  //intersection
  cost.InitialCost(Dim, t, t, CumSumData, CumSumData2, VectOfCosts);
  Radius2 = (VectOfCosts[t + 1] - VectOfCosts[t] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  }
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  Rect -> Intersection_disk(Disk);
  //exclusion
  if (i > 0) {
    unsigned int NbEmptyInter = 0;
    for (unsigned int j = 0; j < i; j++) {
      u = vectlinktocands[j] -> GetTau();
      cost.InitialCost(Dim, u, t-1, CumSumData, CumSumData2, VectOfCosts);
      Radius2 = (VectOfCosts[t] - VectOfCosts[u] - cost.get_coef_Var()) / cost.get_coef();
      Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
      if (Rect -> EmptyIntersection((Disk))) {
        NbEmptyInter++;
      } else {
        Rect -> Exclusion_disk(Disk);
        if (Rect -> IsEmpty_rect()) {
          RealNbExclus = j + 1- NbEmptyInter;
          return;
        }
      }
    }
    RealNbExclus = i - NbEmptyInter;
  }
}
