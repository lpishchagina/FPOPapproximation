// Rect^Tau_t = inter_{j = Tau^t}(S^tau_j) \ S^randBeforeTau_(Tau-1)
//= (Rect^Tau_(t-1)inter S^Tau_t) \ S^randBeforeTau_(Tau-1)
#include "Candidate_Irandom_Eall_10.h"

using namespace Rcpp;
using namespace std;

Candidate_Irandom_Eall_10::Candidate_Irandom_Eall_10(const  Candidate_Irandom_Eall_10 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect = new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Irandom_Eall_10::~Candidate_Irandom_Eall_10() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Candidate_Irandom_Eall_10::GetTau()const { return Tau; }

int Candidate_Irandom_Eall_10::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

void  Candidate_Irandom_Eall_10::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL; VectOfCosts = NULL; }

bool  Candidate_Irandom_Eall_10::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Irandom_Eall_10::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void  Candidate_Irandom_Eall_10::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Irandom_Eall_10>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  //random after Tau: Rect^tau_t = Cube(S^Tau_RandAfterTau
  unsigned int IndexRandAfterTau = get_Number(vectlinktocands.size() - IndexToLinkOfUpdCand) + IndexToLinkOfUpdCand - 1;
  unsigned int RandCandAfterTau = vectlinktocands[IndexRandAfterTau] -> GetTau();
  //pelt
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, RandCandAfterTau, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[RandCandAfterTau + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  }
  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  Rect -> CubeApproximation(Disk);
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
