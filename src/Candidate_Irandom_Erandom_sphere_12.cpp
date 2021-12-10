// Z^tau_t = S^Tau_RandAfterTau\ S^RandBeforeTau_(Tau-1)
#include "Candidate_Irandom_Erandom_sphere_12.h"

using namespace Rcpp;
using namespace std;

Candidate_Irandom_Erandom_sphere_12::Candidate_Irandom_Erandom_sphere_12(const Candidate_Irandom_Erandom_sphere_12 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Irandom_Erandom_sphere_12::~Candidate_Irandom_Erandom_sphere_12(){ CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Irandom_Erandom_sphere_12::GetTau()const { return Tau; }

int Candidate_Irandom_Erandom_sphere_12::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}


double Candidate_Irandom_Erandom_sphere_12::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Candidate_Irandom_Erandom_sphere_12::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_Irandom_Erandom_sphere_12::EmptyOfCandidate() { return fl_empty; }

void Candidate_Irandom_Erandom_sphere_12::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Candidate_Irandom_Erandom_sphere_12::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Irandom_Erandom_sphere_12>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  //pelt
  Cost cost = Cost(Dim);
  unsigned int LastT = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, LastT, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[LastT + 1] - VectOfCosts[Tau] - cost.get_coef_Var()) / cost.get_coef();
  if (Radius2 < 0) {
    fl_empty = true;
    return;
  }
  pSphere Disk_TauLastT = pSphere(Dim);
  Disk_TauLastT.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  //random exclusion
  if (IndexToLinkOfUpdCand > 0) {
    double dist;
    unsigned int IndexRandBeforeTau = get_Number(IndexToLinkOfUpdCand) - 1;
    unsigned int RandBeforeTau = vectlinktocands[IndexRandBeforeTau] -> GetTau();
    cost.InitialCost(Dim, RandBeforeTau, Tau-1, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[Tau] - VectOfCosts[RandBeforeTau] - cost.get_coef_Var()) / cost.get_coef();
    pSphere DiskRandBeforeTauTauMinus1 = pSphere(Dim);
    DiskRandBeforeTauTauMinus1.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
    dist = Dist(Disk_TauLastT.get_center(), DiskRandBeforeTauTauMinus1.get_center());
    if (dist < (DiskRandBeforeTauTauMinus1.get_radius() + Disk_TauLastT.get_radius())) {
      if (dist <= (DiskRandBeforeTauTauMinus1.get_radius() - Disk_TauLastT.get_radius())) {
        fl_empty = true;
      }
    }
  }
}
