#include "Sph_Last_vRand.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Sph_Last_vRand::Sph_Last_vRand(const Sph_Last_vRand & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  fl_empty = candidate.fl_empty;
}

Sph_Last_vRand::~Sph_Last_vRand(){ CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

int Sph_Last_vRand::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

double Sph_Last_vRand::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

unsigned int Sph_Last_vRand::GetTau()const { return Tau; }
void Sph_Last_vRand::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }
bool Sph_Last_vRand::EmptyOfCandidate() { return fl_empty; }

void Sph_Last_vRand::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Sph_Last_vRand::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Last_vRand>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  RealNbExclus = 0;
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {fl_empty = true; return; }  //pelt

  pSphere DiskLastTau = pSphere(Dim);
  DiskLastTau.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));

  //random after tau inclusion
  if (IndexToLinkOfUpdCand > 0) {
    unsigned int IndexRandBeforeTau = get_Number(IndexToLinkOfUpdCand) - 1;
    j = vectlinktocands[IndexRandBeforeTau] -> GetTau();
    cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
    pSphere DiskRandBeforeTauTauMinus1 = pSphere(Dim);
    DiskRandBeforeTauTauMinus1.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
    double dist = Dist(DiskLastTau.get_center(), DiskRandBeforeTauTauMinus1.get_center());
    if (dist < (DiskRandBeforeTauTauMinus1.get_radius() + DiskLastTau.get_radius())) {
      if (dist <= (DiskRandBeforeTauTauMinus1.get_radius() - DiskLastTau.get_radius())) {
        fl_empty = true;
      }
    }
  }
}
