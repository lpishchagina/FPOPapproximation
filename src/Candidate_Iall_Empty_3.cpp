// Rect^Tau_t = inter_{j=i^LastT} S^Tau_j = Rect^Tau_(t-1) inter S^Tau_t
#include "Candidate_Iall_Eempty_3.h"

using namespace Rcpp;
using namespace std;

Candidate_Iall_Eempty_3::Candidate_Iall_Eempty_3(const Candidate_Iall_Eempty_3 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iall_Eempty_3::~Candidate_Iall_Eempty_3(){ delete Rect;  CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Candidate_Iall_Eempty_3::GetTau()const { return Tau; }

void Candidate_Iall_Eempty_3::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL; VectOfCosts = NULL; }

bool Candidate_Iall_Eempty_3::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Iall_Eempty_3::InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Candidate_Iall_Eempty_3::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Iall_Eempty_3>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  double Radius2;
  pSphere Disk = pSphere(Dim);
  Cost cost = Cost(Dim);
  unsigned int LastT = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  unsigned int j;
  //intersection : Rect^Tau_t = Rect^Tau_(t-1) inter_{j=i^LastT} S^Tau_j
  for (unsigned int i = IndexToLinkOfUpdCand; i < vectlinktocands.size(); i++) {
    j = vectlinktocands[i] -> GetTau();
    cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
    //pelt
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
}
