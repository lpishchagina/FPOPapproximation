// Rect^tau_t = (Cube(S^Tau_LastT) inter S^Tau_RandAfterTau)\ S^RandBeforeTau_(Tau-1)})
#include "Candidate_Irandom_Erandom_8.h"

using namespace Rcpp;
using namespace std;

Candidate_Irandom_Erandom_8::Candidate_Irandom_Erandom_8(const Candidate_Irandom_Erandom_8 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Irandom_Erandom_8::~Candidate_Irandom_Erandom_8() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL; VectOfCosts = NULL; }

unsigned int Candidate_Irandom_Erandom_8::GetTau()const { return Tau; }

int Candidate_Irandom_Erandom_8::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

void Candidate_Irandom_Erandom_8::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL; VectOfCosts = NULL; }

bool Candidate_Irandom_Erandom_8::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Irandom_Erandom_8::InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Candidate_Irandom_Erandom_8::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Irandom_Erandom_8>::iterator> &vectlinktocands, unsigned int& RealNbExclus){
  RealNbExclus = 0;
  //random after Tau: Rect^tau_t = (Rect^tau_(t-1) inter S^Tau_RandAfterTau)
  unsigned int IndexRandAfterTau = get_Number (vectlinktocands.size() - IndexToLinkOfUpdCand) + IndexToLinkOfUpdCand - 1;
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
  Rect -> Intersection_disk(Disk);
  //random exclusion : Rect^Tau_t = Rect^Tau_t  \ S^RandBeforeTau_(Tau-1)
  if ((IndexToLinkOfUpdCand > 0) && (!Rect -> IsEmpty_rect())) {
    unsigned int IndexRandBeforeTau = get_Number(IndexToLinkOfUpdCand) - 1;
    unsigned int RandCandBeforeTau = vectlinktocands[IndexRandBeforeTau] -> GetTau();
    cost.InitialCost(Dim, RandCandBeforeTau, Tau-1, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[Tau] - VectOfCosts[RandCandBeforeTau] - cost.get_coef_Var()) / cost.get_coef();
    Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
    if (!Rect -> EmptyIntersection(Disk)) {
      Rect -> Exclusion_disk(Disk);
      RealNbExclus = 1;
    }
  }
}
