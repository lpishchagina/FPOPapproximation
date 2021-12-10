// Z^tau_t = S^Tau_LastT\ S^RandBeforeTau_(Tau-1)
#include "Candidate_Ilast_Erandom_sphere_13.h"

using namespace Rcpp;
using namespace std;

Candidate_Ilast_Erandom_sphere_13::Candidate_Ilast_Erandom_sphere_13(const Candidate_Ilast_Erandom_sphere_13 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Ilast_Erandom_sphere_13::~Candidate_Ilast_Erandom_sphere_13(){ CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Ilast_Erandom_sphere_13::GetTau()const { return Tau; }

int Candidate_Ilast_Erandom_sphere_13::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}


double Candidate_Ilast_Erandom_sphere_13::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Candidate_Ilast_Erandom_sphere_13::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_Ilast_Erandom_sphere_13::EmptyOfCandidate() { return fl_empty; }

void Candidate_Ilast_Erandom_sphere_13::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Candidate_Ilast_Erandom_sphere_13::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Ilast_Erandom_sphere_13>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  RealNbExclus = 0;
  //random after Tau: Rect^tau_t = Cube(S^Tau_RandAfterTau)
  unsigned int IndexRandAfterTau = get_Number (vectlinktocands.size() - IndexToLinkOfUpdCand) + IndexToLinkOfUpdCand - 1;
  unsigned int RandCandAfterTau = vectlinktocands[IndexRandAfterTau] -> GetTau();
  //pelt
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, RandCandAfterTau, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[RandCandAfterTau + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    fl_empty = true;
    return;
  }
  //random after tau
  pSphere DiskTauRandAfterTau = pSphere(Dim);
  DiskTauRandAfterTau.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  //random exclusion
  if (IndexToLinkOfUpdCand > 0) {
    double dist;
    unsigned int IndexRandBeforeTau = get_Number(IndexToLinkOfUpdCand) - 1;
    unsigned int RandBeforeTau = vectlinktocands[IndexRandBeforeTau] -> GetTau();
    cost.InitialCost(Dim, RandBeforeTau, Tau-1, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[Tau] - VectOfCosts[RandBeforeTau] - cost.get_coef_Var()) / cost.get_coef();
    pSphere DiskRandBeforeTauTauMinus1 = pSphere(Dim);
    DiskRandBeforeTauTauMinus1.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
    dist = Dist(DiskTauRandAfterTau.get_center(), DiskRandBeforeTauTauMinus1.get_center());
    if (dist < (DiskRandBeforeTauTauMinus1.get_radius() + DiskTauRandAfterTau.get_radius())) {
      if (dist <= (DiskRandBeforeTauTauMinus1.get_radius() - DiskTauRandAfterTau.get_radius())) {
        fl_empty = true;
      }
    }
  }
}
