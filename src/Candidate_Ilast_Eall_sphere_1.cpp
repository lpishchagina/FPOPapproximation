// Z^tau_t = S^Tau_LastT\ (union_{j=1^Tau-1}S^j_(Tau-1))
#include "Candidate_Ilast_Eall_sphere_1.h"

using namespace Rcpp;
using namespace std;

Candidate_Ilast_Eall_sphere_1::Candidate_Ilast_Eall_sphere_1(const Candidate_Ilast_Eall_sphere_1 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Candidate_Ilast_Eall_sphere_1::~Candidate_Ilast_Eall_sphere_1(){ CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Ilast_Eall_sphere_1::GetTau()const { return Tau; }

double Candidate_Ilast_Eall_sphere_1::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Candidate_Ilast_Eall_sphere_1::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Candidate_Ilast_Eall_sphere_1::EmptyOfCandidate() { return fl_empty; }

void Candidate_Ilast_Eall_sphere_1::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Candidate_Ilast_Eall_sphere_1::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Ilast_Eall_sphere_1>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
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
      dist = Dist(Disk_TauLastT.get_center(), DiskJTauMinus1.get_center());
      if (dist < (DiskJTauMinus1.get_radius() + Disk_TauLastT.get_radius())) {
        if (dist <= (DiskJTauMinus1.get_radius() - Disk_TauLastT.get_radius())) {
          fl_empty = true;
          return;
        }
      }
    }
  }
}
