#include "Sph_Rand_vAll.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Sph_Rand_vAll::Sph_Rand_vAll(const Sph_Rand_vAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
}

Sph_Rand_vAll::~Sph_Rand_vAll(){ CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Sph_Rand_vAll::GetTau()const { return Tau; }

int Sph_Rand_vAll::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

double Sph_Rand_vAll::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Sph_Rand_vAll::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }

bool Sph_Rand_vAll::EmptyOfCandidate() { return fl_empty; }

void Sph_Rand_vAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Sph_Rand_vAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Rand_vAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  RealNbExclus = 0;
  //random after Tau:
  unsigned int IndexRandAfterTau = get_Number (vectlinktocands.size() - IndexToLinkOfUpdCand) + IndexToLinkOfUpdCand - 1;
  unsigned int j = vectlinktocands[IndexRandAfterTau] -> GetTau();

  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

  if (Radius2 < 0) {fl_empty = true; return; } //pelt

  pSphere Disk_TauRandCandAfterTau = pSphere(Dim);
  Disk_TauRandCandAfterTau.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  //inclusion :
  if (IndexToLinkOfUpdCand > 0)  {
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

