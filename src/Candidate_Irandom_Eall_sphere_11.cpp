// Z^tau_t = S^Tau_LastT\ (union_{j=1^Tau-1}S^j_(Tau-1))
#include "Candidate_Irandom_Eall_sphere_11.h"

using namespace Rcpp;
using namespace std;

Candidate_Irandom_Eall_sphere_11::Candidate_Irandom_Eall_sphere_11(const Candidate_Irandom_Eall_sphere_11 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Irandom_Eall_sphere_11::~Candidate_Irandom_Eall_sphere_11(){ CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Irandom_Eall_sphere_11::GetTau()const { return Tau; }

int Candidate_Irandom_Eall_sphere_11::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}


double Candidate_Irandom_Eall_sphere_11::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Candidate_Irandom_Eall_sphere_11::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_Irandom_Eall_sphere_11::EmptyOfCandidate() { return fl_empty; }

void Candidate_Irandom_Eall_sphere_11::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Candidate_Irandom_Eall_sphere_11::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Irandom_Eall_sphere_11>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
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
  //random
  pSphere Disk_TauRandCandAfterTau = pSphere(Dim);
  Disk_TauRandCandAfterTau.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  //exclusion
  if (IndexToLinkOfUpdCand > 0)  {
    unsigned int j;
    pSphere DiskJTauMinus1 = pSphere(Dim);
    double dist;
    for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
      j = vectlinktocands[i] -> GetTau();
      cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
      Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
      DiskJTauMinus1.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
      dist = Dist(Disk_TauRandCandAfterTau.get_center(), DiskJTauMinus1.get_center());
      if (dist < (DiskJTauMinus1.get_radius() + Disk_TauRandCandAfterTau.get_radius())) {
        if (dist <= (DiskJTauMinus1.get_radius() - Disk_TauRandCandAfterTau.get_radius())) {
          fl_empty = true;
          return;
        }
      }
    }
  }
}

