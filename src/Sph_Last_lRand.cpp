#include "Sph_Last_lRand.h"
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

Sph_Last_lRand::Sph_Last_lRand(const Sph_Last_lRand & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  IndexVectBefore.clear();
  IndexVectBefore = candidate.IndexVectBefore;
  CreationFl = candidate.CreationFl;

}

Sph_Last_lRand::~Sph_Last_lRand(){ CumSumData = NULL; CumSumData2 = NULL;  VectOfCosts = NULL; }

int Sph_Last_lRand::get_Number(int N) {
  srand(time(NULL));
  int res = rand() % N + 1;
  return res;
}

double Sph_Last_lRand::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

unsigned int Sph_Last_lRand::GetTau()const { return Tau; }
void Sph_Last_lRand::CleanOfCandidate() { CumSumData = NULL;  VectOfCosts = NULL; }
bool Sph_Last_lRand::EmptyOfCandidate() { return fl_empty; }
std::vector<unsigned int> Sph_Last_lRand::GetIndexVectBefore()const {return IndexVectBefore;}

void Sph_Last_lRand::InitialOfCandidate(unsigned int tau, double** &cumsumdata, double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  fl_empty = false;
  IndexVectBefore.clear();
  CreationFl = true;
}

void Sph_Last_lRand::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Sph_Last_lRand>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  fl_empty = false;
  RealNbExclus = 0;
  Cost cost = Cost(Dim);
  unsigned int j = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, j, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[j + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();

  if (Radius2 < 0) { fl_empty = true; return;}   //pelt

  pSphere Disk_TauLastT = pSphere(Dim);
  Disk_TauLastT.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));

  //CreationFl = true =>1 iteration : Creation of IndexVectBefore
  if (CreationFl) {
    CreationFl = false;
    if (IndexToLinkOfUpdCand > 0) {
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        IndexVectBefore.push_back(vectlinktocands[i] -> GetTau());
      }
    }
  }
  //inclusion :
  if ((!fl_empty) &&(IndexVectBefore.size() > 0)) {
    unsigned int IndexRandBeforeTau = get_Number(IndexVectBefore.size()) - 1;
    j = IndexVectBefore[IndexRandBeforeTau];
    cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
    Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
    pSphere DiskRandBeforeTau = pSphere(Dim);
    DiskRandBeforeTau.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));

    double dist = Dist(Disk_TauLastT.get_center(),DiskRandBeforeTau.get_center());
    if (dist < (DiskRandBeforeTau.get_radius() + Disk_TauLastT.get_radius())){
      if (dist <= (DiskRandBeforeTau.get_radius() - Disk_TauLastT.get_radius())){
        fl_empty = true;
        return;
      }
    } else {
      if (IndexRandBeforeTau < (IndexVectBefore.size() -1 )) {
        IndexVectBefore[IndexRandBeforeTau] = IndexVectBefore.back();
      }
      IndexVectBefore.pop_back();
    }
  }
}
