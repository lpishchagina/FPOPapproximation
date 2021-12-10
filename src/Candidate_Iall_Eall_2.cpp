// Rect^Tau_t = inter_{j = Tau^t}(S^tau_j) \ (union_{j=1^Tau-1}S^j_(Tau-1))
//= (Rect^Tau_(t-1)inter S^Tau_t) \ (union_{j=1^Tau-1}S^j_(Tau-1))
#include "Candidate_Iall_Eall_2.h"

using namespace Rcpp;
using namespace std;

Candidate_Iall_Eall_2::Candidate_Iall_Eall_2(const Candidate_Iall_Eall_2 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Iall_Eall_2::~Candidate_Iall_Eall_2() { delete Rect;  CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Candidate_Iall_Eall_2::GetTau()const { return Tau; }

void Candidate_Iall_Eall_2::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

bool Candidate_Iall_Eall_2::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Iall_Eall_2::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Candidate_Iall_Eall_2::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Iall_Eall_2>::iterator> &vectlinktocands, unsigned int &RealNbExclus) {
  RealNbExclus = 0;
  //pelt
  Cost cost = Cost(Dim);
  unsigned int LastT = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, LastT, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[LastT + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  }
  //intersection : Rect^Tau_t = Rect^Tau_(t-1) inter S^Tau_t
  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  Rect -> Intersection_disk(Disk);
  //exclusion : Rect^Tau_t= Rect^Tau_t \ (union_{j=1^Tau-1}S^j_(Tau-1))
  if ((IndexToLinkOfUpdCand > 0) && (!Rect -> IsEmpty_rect())) {
    unsigned int j;
    unsigned int NbEmptyInter = 0;
    for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
      j = vectlinktocands[i] -> GetTau();
      cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
      Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
      Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
      if (Rect -> EmptyIntersection((Disk))) {
        NbEmptyInter++;
      } else {
        Rect -> Exclusion_disk(Disk);
        if (Rect -> IsEmpty_rect()) {
          RealNbExclus = i + 1 - NbEmptyInter;
          return;
        }
      }
    }
    RealNbExclus = IndexToLinkOfUpdCand - NbEmptyInter;
  }
}
