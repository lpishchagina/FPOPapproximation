// Rect^tau_t = (Rect^tau-(LastT-1) interS^Tau_LastT)\ (union_{j=1^Tau-1}S^j_(Tau-1))
#include "Candidate_Ilast_EallModif2_15.h"

using namespace Rcpp;
using namespace std;

Candidate_Ilast_EallModif2_15::Candidate_Ilast_EallModif2_15(const Candidate_Ilast_EallModif2_15 & candidate) {
  Dim = candidate.Dim;
  Tau = candidate.Tau;
  Rect= new pRectangle(Dim);
  CumSumData = candidate.CumSumData;
  CumSumData2 = candidate.CumSumData2;
  VectOfCosts = candidate.VectOfCosts;
  disks_t_1.clear();
  disks_t_1 = candidate.disks_t_1;
  CreationFl = candidate.CreationFl;
}

Candidate_Ilast_EallModif2_15::~Candidate_Ilast_EallModif2_15() { delete Rect;  CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; }

unsigned int Candidate_Ilast_EallModif2_15::GetTau()const { return Tau; }

std::list<pSphere> Candidate_Ilast_EallModif2_15::get_disks_t_1()const { return disks_t_1; }

void Candidate_Ilast_EallModif2_15::CleanOfCandidate() { CumSumData = NULL;  CumSumData2 = NULL;  VectOfCosts = NULL; disks_t_1.clear();}

bool Candidate_Ilast_EallModif2_15::EmptyOfCandidate() { return Rect -> IsEmpty_rect(); }

void Candidate_Ilast_EallModif2_15::InitialOfCandidate(unsigned int tau, double** &cumsumdata,  double** &cumsumdata2, double* &vectofcosts) {
  Tau = tau;
  CumSumData = cumsumdata;
  CumSumData2 = cumsumdata2;
  VectOfCosts = vectofcosts;
  disks_t_1.clear();
  CreationFl = true;
}

double Candidate_Ilast_EallModif2_15::Dist(double* a, double*b) {
  double dist = 0;
  for (unsigned int k = 0; k < Dim; k++) {
    dist = dist + (a[k] - b[k])*(a[k] - b[k]);
  }
  return sqrt(dist);
}

void Candidate_Ilast_EallModif2_15::UpdateOfCandidate(unsigned int IndexToLinkOfUpdCand, std::vector<std::list<Candidate_Ilast_EallModif2_15>::iterator> &vectlinktocands, unsigned int& RealNbExclus) {
  RealNbExclus = 0;
  //check 1 : pelt
  Cost cost = Cost(Dim);
  unsigned int LastT = vectlinktocands[vectlinktocands.size() - 1] -> GetTau();
  cost.InitialCost(Dim, Tau, LastT, CumSumData, CumSumData2, VectOfCosts);
  double Radius2 = (VectOfCosts[LastT + 1] - VectOfCosts[Tau] - cost.get_coef_Var())/cost.get_coef();
  if (Radius2 < 0) {
    Rect -> DoEmpty_rect();
    return;
  }

  pSphere Disk = pSphere(Dim);
  Disk.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
  unsigned int j;
  std::list<pSphere>::iterator iter;
  //Creation of list of disks(tau-1) for 1 iteration
  if (CreationFl) {
    CreationFl = false;
    if (IndexToLinkOfUpdCand > 0) {
      double dist;
      pSphere DiskTau_1 = pSphere(Dim);
      for (unsigned int i = 0; i < IndexToLinkOfUpdCand; i++) {
        j = vectlinktocands[i] -> GetTau();
        cost.InitialCost(Dim, j, Tau-1, CumSumData, CumSumData2, VectOfCosts);
        Radius2 = (VectOfCosts[Tau] - VectOfCosts[j] - cost.get_coef_Var()) / cost.get_coef();
        DiskTau_1.InitialpSphere(Dim, cost.get_mu(), sqrt(Radius2));
        dist = Dist(Disk.get_center(), DiskTau_1.get_center());
        if (dist < (DiskTau_1.get_radius() + Disk.get_radius())) {
          if (dist <= (DiskTau_1.get_radius() - Disk.get_radius())) {
           Rect -> DoEmpty_rect();
            return;
          }
        }
        disks_t_1.push_back(DiskTau_1);
      }
    }
  }
  //intersection
  Rect -> Intersection_disk(Disk);
  //exclusion
  if ((disks_t_1.size() > 0) && (!Rect -> IsEmpty_rect())) {
    iter = disks_t_1.begin();
    while(iter != disks_t_1.end() && (!Rect -> IsEmpty_rect())){
      if (Rect -> EmptyIntersection(*iter)) {
        iter = disks_t_1.erase(iter);
      }//isn't intersection => Remove disks
      else {
        Rect -> Exclusion_disk(*iter);
        RealNbExclus++;
        ++iter;
      }
    }
  }
}
