// Rect^tau_t = Cube(S^tau_t)\ S^RandBeforeTau_(Tau-1)})
#include "Candidate_Ilast_Erandom_6.h"

using namespace Rcpp;
using namespace std;

Candidate_Ilast_Erandom_6::Candidate_Ilast_Erandom_6(const Candidate_Ilast_Erandom_6 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect = new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Ilast_Erandom_6::~Candidate_Ilast_Erandom_6() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Ilast_Erandom_6::GetTau()const { return Tau; }

int Candidate_Ilast_Erandom_6::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

void Candidate_Ilast_Erandom_6::CleanOfCandidate() { CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

bool Candidate_Ilast_Erandom_6::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Ilast_Erandom_6::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
}

void Candidate_Ilast_Erandom_6::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Ilast_Erandom_6>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  //pelt
  Cost cost = Cost(Dim);
  unsigned int tLast = vectlinktocands[vectlinktocands.size() -1] -> GetTau();
  cost.InitialCost(Dim, Tau, tLast, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[tLast + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  }
  //last : Rect^Tau_t = Cube(S^Tau_LastT)
  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  Rect -> CubeApproximation(Disk);
  //exclusion : Rect^Tau_t= Rect^Tau_t \ S^RandBeforeTau_(Tau-1)})
  if (IndexToLinkOfUpdCand > 0) {
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
