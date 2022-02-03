#include "Sph_Rand_lAll.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

Sph_Rand_lAll::Sph_Rand_lAll(const Sph_Rand_lAll & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  fl_empty = candidate.fl_empty;
  DiskListBefore.clear();
  DiskListBefore = candidate.DiskListBefore;
  CreationFl = candidate.CreationFl;
}

Sph_Rand_lAll::~Sph_Rand_lAll() { CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }


unsigned int Sph_Rand_lAll::GetTau()const { return Tau; }
void Sph_Rand_lAll::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }
bool Sph_Rand_lAll::EmptyOfCandidate() { return fl_empty; }

int Sph_Rand_lAll::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

double Sph_Rand_lAll::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) { dist = dist + (a[k] - b[k])*(a[k] - b[k]); }
  return sqrt(dist);
}

void Sph_Rand_lAll::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
  DiskListBefore.clear();
  CreationFl = true;
}

void Sph_Rand_lAll::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Rand_lAll>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  RealNbExclus = 0;
  //random after Tau:
  unsigned int IndexRandAfterTau = get_Number (vectlinktocands.size() - IndexToLinkOfUpdCand) + IndexToLinkOfUpdCand - 1;
  unsigned int j = vectlinktocands[IndexRandAfterTau] -> GetTau();
  Cost cost = Cost(Dim);
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

  if (Radius2 < 0) {  fl_empty = true;  return; }   //pelt
  //random
  pSphere Disk_TauRandCandAfterTau = pSphere(Dim);
  Disk_TauRandCandAfterTau.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));

  double dist;
  //CreationFl = true =>1 iteration : Creation of DiskListBefore
  if (CreationFl) {
    CreationFl = false;
    if (IndexToLinkOfUpdCand > 0) {
      pSphere DiskTau_1 = pSphere(Dim);
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        j = vectlinktocands[i] -> GetTau();
        cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
        Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
        DiskTau_1.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
        /*//check for 1 iteration: inclusion
        dist = Dist(Disk_TauRandCandAfterTau.get_center(), DiskTau_1.get_center());
        if (dist < (DiskTau_1.get_radius() + Disk_TauRandCandAfterTau.get_radius())) {
          if (dist <= (DiskTau_1.get_radius() - Disk_TauRandCandAfterTau.get_radius())) { fl_empty = true; return; }
        }*/
        DiskListBefore.push_back(DiskTau_1);
      }
    }
  }
  //inclusion :
  if ((!fl_empty) && (DiskListBefore.size() > 0)) {
    std::list<pSphere>::iterator iter = DiskListBefore.begin();
    while( iter != DiskListBefore.end()){
      dist = Dist(Disk_TauRandCandAfterTau.get_center(),(*iter).get_center());
      if (dist < ((*iter).get_radius() + Disk_TauRandCandAfterTau.get_radius())){
        if (dist <= ((*iter).get_radius() - Disk_TauRandCandAfterTau.get_radius())){
          fl_empty = true;
          return;
        } else {
          ++iter;
        }
      } else {
        iter = DiskListBefore.erase(iter);
      }
    }
  }
}
