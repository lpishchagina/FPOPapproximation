#include "Sph_Last_vAll.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

Sph_Last_vAll::Sph_Last_vAll(const Sph_Last_vAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  fl_empty = candidate.fl_empty;
}

Sph_Last_vAll::~Sph_Last_vAll(){ CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Sph_Last_vAll::GetTau()const { return Tau; }

double Sph_Last_vAll::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) { dist = dist + (a[k] - b[k])*(a[k] - b[k]);}
  return sqrt(dist);
}

void Sph_Last_vAll::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }
bool Sph_Last_vAll::EmptyOfCandidate() { return fl_empty; }

void Sph_Last_vAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
}

void Sph_Last_vAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Last_vAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  RealNbExclus = 0;
  Cost cost = Cost(Dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var()) / cost.get_coef();

  if (Radius2 < 0) { fl_empty = true; return; } //pelt

  pSphere Disk_TauLastT = pSphere(Dim);
  Disk_TauLastT.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  //check inclusion
  if (IndexToLinkOfUpdCand > 0)  {
    pSphere DiskJTauMinus1 = pSphere(Dim);
    double dist;
    for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
      ++RealNbExclus;
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
